from mbc_option_common import prepare_with_previous_params
import datetime
import csv
import io
import os
import re
import time
import threading
import urllib.error
import urllib.parse
import urllib.request
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed


BLAST_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
INITIAL_WAIT = 15
POLL_INTERVAL = 10
MAX_WAIT = 3600
BATCH_SIZE = 50
MIN_VISIBLE_PROGRESS_PAUSE = 0.25
MIAN_RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
NON_INFORMATIVE_TAXA = {
    "uncultured",
    "uncultured bacterium",
    "uncultured archaeon",
    "uncultured fungus",
    "uncultured eukaryote",
    "environmental samples",
    "unclassified sequences",
    "unclassified",
    "unidentified",
    "unknown",
    "metagenome",
}
DEFAULT_MAX_WORKERS = 3
API_KEY_MAX_WORKERS = 10


def _has_any_taxonomy_value(lineage):
    return any((lineage.get(rank) or "").strip() != "" for rank in MIAN_RANKS)


def _has_valid_ncbi_api_key(api_key):
    key = (api_key or "").strip()
    return len(key) >= 20 and " " not in key


def _ncbi_rate_limit(api_key):
    return 10.0 if _has_valid_ncbi_api_key(api_key) else 3.0


def _ncbi_max_workers(api_key):
    return API_KEY_MAX_WORKERS if _has_valid_ncbi_api_key(api_key) else DEFAULT_MAX_WORKERS


def _blast_post(params):
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(BLAST_URL, data=data, method="POST")
    req.add_header("Content-Type", "application/x-www-form-urlencoded")
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            return resp.read().decode("utf-8", errors="replace")
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        raise RuntimeError(f"Network error while contacting NCBI BLAST (submit): {e}")


def _blast_get(params):
    qs = urllib.parse.urlencode(params)
    url = f"{BLAST_URL}?{qs}"
    try:
        with urllib.request.urlopen(url, timeout=120) as resp:
            return resp.read()
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        raise RuntimeError(f"Network error while contacting NCBI BLAST (retrieve): {e}")


def _blast_submit(core, fasta_text, db, evalue, hitlist_size, api_key):
    params = {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "DATABASE": db,
        "QUERY": fasta_text,
        "EXPECT": str(evalue),
        "HITLIST_SIZE": str(hitlist_size),
        "FORMAT_TYPE": "XML",
        "TOOL": "mbctools",
    }
    if api_key:
        params["API_KEY"] = api_key

    print(core.warningStyle + "\nSubmitting query to NCBI BLASTn..." + core.normalStyle)
    try:
        response = _blast_post(params)
    except RuntimeError as e:
        print(core.errorStyle + str(e) + core.normalStyle)
        return None, None

    rid = None
    rtoe = None
    for line in response.splitlines():
        if line.strip().startswith("RID ="):
            rid = line.split("=", 1)[1].strip()
        elif line.strip().startswith("RTOE ="):
            try:
                rtoe = int(line.split("=", 1)[1].strip())
            except ValueError:
                rtoe = INITIAL_WAIT

    if not rid:
        print(core.errorStyle + "Unable to parse RID from NCBI response." + core.normalStyle)
        return None, None

    print(core.successStyle + f"RID: {rid} / RTOE: {rtoe if rtoe is not None else INITIAL_WAIT}s" + core.normalStyle)
    return rid, (rtoe if rtoe is not None else INITIAL_WAIT)


def _blast_wait_for_results(core, rid, rtoe):
    wait = max(rtoe, INITIAL_WAIT)
    print(core.warningStyle + f"Waiting {wait}s before first status check..." + core.normalStyle)
    time.sleep(wait)

    elapsed = wait
    while elapsed < MAX_WAIT:
        try:
            raw = _blast_get(
                {
                    "CMD": "Get",
                    "RID": rid,
                    "FORMAT_TYPE": "XML",
                    "FORMAT_OBJECT": "SearchInfo",
                }
            ).decode("utf-8", errors="replace")
        except RuntimeError as e:
            print(core.errorStyle + str(e) + core.normalStyle)
            return False

        if "Status=WAITING" in raw:
            print(core.warningStyle + f"Status: WAITING (elapsed {elapsed}s)..." + core.normalStyle)
            time.sleep(POLL_INTERVAL)
            elapsed += POLL_INTERVAL
        elif "Status=READY" in raw:
            if "ThereAreHits=yes" in raw:
                print(core.successStyle + "Status: READY - hits found" + core.normalStyle)
            else:
                print(core.warningStyle + "Status: READY - no hits found" + core.normalStyle)
            return True
        elif "Status=FAILED" in raw:
            print(core.errorStyle + "BLAST search failed on NCBI side." + core.normalStyle)
            return False
        elif "Status=UNKNOWN" in raw:
            print(core.errorStyle + "RID is unknown or expired (>24h)." + core.normalStyle)
            return False
        else:
            print(core.warningStyle + f"Unexpected status response (elapsed {elapsed}s), retrying..." + core.normalStyle)
            time.sleep(POLL_INTERVAL)
            elapsed += POLL_INTERVAL

    print(core.errorStyle + f"BLAST timed out after {MAX_WAIT}s." + core.normalStyle)
    return False


def _blast_retrieve(rid, fmt_type, extra=None):
    params = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": fmt_type,
        "HITLIST_SIZE": "500",
    }
    if fmt_type != "Text":
        params["FORMAT_OBJECT"] = "Alignment"
    if extra:
        params.update(extra)
    return _blast_get(params)


def _extract_pre_block(text):
    start_tag = "<PRE>"
    end_tag = "</PRE>"
    start = text.find(start_tag)
    end = text.find(end_tag)
    if start >= 0 and end > start:
        return text[start + len(start_tag):end].strip()
    return text.strip()


def _download_blast_tabular(core, rid):
    # NCBI sometimes returns a temporary HTML wrapper (<PRE>...</PRE>) before tabular content is ready.
    attempt_params = [
        {"ALIGNMENT_VIEW": "Tabular", "FORMAT_OBJECT": "Alignment"},
        {"ALIGNMENT_VIEW": "Tabular"},
    ]
    for i in range(6):
        params = attempt_params[i % len(attempt_params)]
        try:
            raw = _blast_get(
                {
                    "CMD": "Get",
                    "RID": rid,
                    "FORMAT_TYPE": "Text",
                    "HITLIST_SIZE": "500",
                    **params,
                }
            )
        except RuntimeError as e:
            print(core.errorStyle + str(e) + core.normalStyle)
            return None

        text = raw.decode("utf-8", errors="replace").strip()
        text = _extract_pre_block(text)
        if "# Fields:" in text:
            return (text + "\n").encode("utf-8")

        if i < 5:
            print(core.warningStyle + "Tabular BLAST output not ready yet, retrying in 5s..." + core.normalStyle)
            time.sleep(5)

    print(core.errorStyle + "Unable to retrieve a valid BLAST tabular hit table from NCBI." + core.normalStyle)
    return None


def _extract_fasta_records(fasta_text):
    records = []
    current = []
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            if len(current) > 0:
                records.append("\n".join(current))
            current = [line]
        elif len(current) > 0:
            current.append(line)
    if len(current) > 0:
        records.append("\n".join(current))
    return records


def _normalize_field_name(name):
    return re.sub(r"[^a-z0-9]", "", (name or "").lower())


def _parse_hit_table_field_names(hit_table_lines):
    for line in hit_table_lines:
        if line.startswith("#") and "Fields:" in line:
            after_colon = line.split(":", 1)[1]
            return [field.strip() for field in after_colon.split(",")]
    return []


def _find_field_index(field_names, candidates):
    normalized_names = [_normalize_field_name(name) for name in field_names]
    normalized_candidates = {_normalize_field_name(candidate) for candidate in candidates}
    for idx, field_name in enumerate(normalized_names):
        if field_name in normalized_candidates:
            return idx
    return None


def _reverse_complement_dna(seq):
    return seq.translate(str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn"))[::-1]


def _extract_record_accession(record):
    header = record.splitlines()[0][1:].strip()
    first_token = header.split()[0]
    if "|" in first_token:
        split_pipe = [part for part in first_token.split("|") if part != ""]
        for part in reversed(split_pipe):
            if "." in part or re.match(r"^[A-Za-z_]+\d+", part):
                first_token = part
                break
    return first_token


def _fetch_full_fasta_records_from_accessions(accessions, entrez_db, api_key):
    if len(accessions) == 0:
        return {}

    params = {
        "db": entrez_db,
        "id": ",".join(accessions),
        "rettype": "fasta",
        "retmode": "text",
        "tool": "mbctools",
    }
    if api_key:
        params["api_key"] = api_key

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=60) as resp:
        text = resp.read().decode("utf-8", errors="replace")

    records = _extract_fasta_records(text)
    sequence_map = {}
    for record in records:
        accession = _extract_record_accession(record)
        sequence = "".join(record.splitlines()[1:]).replace(" ", "").upper()
        if accession != "" and sequence != "":
            sequence_map[accession] = sequence
            if "." in accession:
                sequence_map[accession.split(".", 1)[0]] = sequence
    return sequence_map


def _download_subject_aligned_regions_fasta(core, hit_table_lines, db_prefix, output_path, api_key):
    field_names = _parse_hit_table_field_names(hit_table_lines)
    if len(field_names) == 0:
        raise ValueError("Unable to parse '# Fields:' line from hit table")

    subject_index = _find_field_index(field_names, ["subject acc.ver", "subject id", "sseqid"])
    sstart_index = _find_field_index(field_names, ["s. start", "subject start", "sstart"])
    send_index = _find_field_index(field_names, ["s. end", "subject end", "send"])

    if None in [subject_index, sstart_index, send_index]:
        raise ValueError(
            "Hit table must contain subject accession and subject start/end coordinates to export aligned intervals"
        )

    distinct_regions = []
    seen = set()
    for line in hit_table_lines:
        if line.startswith("#"):
            continue
        split_line = line.strip().split("\t")
        if len(split_line) <= max(subject_index, sstart_index, send_index):
            continue

        accession = split_line[subject_index].strip()
        if accession == "":
            continue
        try:
            sstart = int(split_line[sstart_index].strip())
            send = int(split_line[send_index].strip())
        except ValueError:
            continue

        region_key = (accession, sstart, send)
        if region_key in seen:
            continue
        seen.add(region_key)
        distinct_regions.append(region_key)

    if len(distinct_regions) == 0:
        raise ValueError("No subject aligned regions found in hit table")

    entrez_db = "nuccore" if db_prefix == "n" else "protein"
    sleep_between_requests = 0.11 if api_key else 0.34
    visible_pause = max(sleep_between_requests, MIN_VISIBLE_PROGRESS_PAUSE)
    total_regions = len(distinct_regions)
    unique_accessions = []
    for accession, _, _ in distinct_regions:
        if accession not in unique_accessions:
            unique_accessions.append(accession)

    accession_to_sequence = {}
    total_accession_batches = max(1, (len(unique_accessions) + BATCH_SIZE - 1) // BATCH_SIZE)
    full_seq_phase_start = time.time()
    for batch_i, accession_batch in enumerate(_chunked(unique_accessions, BATCH_SIZE), start=1):
        time.sleep(visible_pause)
        try:
            accession_to_sequence.update(_fetch_full_fasta_records_from_accessions(accession_batch, entrez_db, api_key))
        except Exception:
            pass

        elapsed = int(time.time() - full_seq_phase_start)
        print(
            core.warningStyle
            + f"Subject full-sequence batch retrieval progress: {batch_i}/{total_accession_batches} (elapsed {elapsed}s)"
            + core.normalStyle,
            flush=True,
        )

    missing_accessions = [accession for accession in unique_accessions if accession_to_sequence.get(accession) in [None, ""]]
    if len(missing_accessions) > 0:
        print(
            core.warningStyle
            + f"Retrying missing subject full sequences individually: {len(missing_accessions)} accession(s)"
            + core.normalStyle,
            flush=True,
        )
    for accession in missing_accessions:
        time.sleep(visible_pause)
        try:
            accession_to_sequence.update(_fetch_full_fasta_records_from_accessions([accession], entrez_db, api_key))
        except Exception:
            pass

    accession_to_taxid = {}
    if db_prefix == "n":
        # Share one limiter across taxid and lineage phases, like OTU taxonomy build.
        rate = _ncbi_rate_limit(api_key)
        limiter = _RateLimiter(rate)

        if len(unique_accessions) > 0:
            def _taxid_progress(completed, total_batches):
                elapsed = int(time.time() - full_seq_phase_start)
                print(
                    core.warningStyle
                    + f"Subject taxid batch retrieval progress: {completed}/{total_batches} (elapsed {elapsed}s)"
                    + core.normalStyle,
                    flush=True,
                )

            accession_to_taxid.update(
                _fetch_taxids_from_accessions(
                    unique_accessions,
                    api_key,
                    limiter=limiter,
                    progress_callback=_taxid_progress,
                )
            )

            missing_taxid_accessions = [accession for accession in unique_accessions if accession_to_taxid.get(accession) in [None, ""]]
            if len(missing_taxid_accessions) > 0:
                print(
                    core.warningStyle
                    + f"Retrying missing subject taxids individually: {len(missing_taxid_accessions)} accession(s)"
                    + core.normalStyle,
                    flush=True,
                )
            for accession in missing_taxid_accessions:
                limiter.acquire()
                accession_to_taxid[accession] = _fetch_taxid_from_accession(accession, api_key)

    unique_taxids = []
    for taxid in accession_to_taxid.values():
        if taxid not in [None, ""] and taxid not in unique_taxids:
            unique_taxids.append(taxid)

    taxid_to_lineage = {}
    if len(unique_taxids) > 0:
        def _lineage_progress(completed, total_batches):
            elapsed = int(time.time() - full_seq_phase_start)
            print(
                core.warningStyle
                + f"Subject lineage batch retrieval progress: {completed}/{total_batches} (elapsed {elapsed}s)"
                + core.normalStyle,
                flush=True,
            )

        taxid_to_lineage.update(
            _fetch_lineages_from_taxids(
                unique_taxids,
                api_key,
                limiter=limiter,
                progress_callback=_lineage_progress,
            )
        )

        missing_lineage_taxids = [
            taxid for taxid in unique_taxids if not _has_any_taxonomy_value(taxid_to_lineage.get(taxid, {rank: "" for rank in MIAN_RANKS}))
        ]
        if len(missing_lineage_taxids) > 0:
            print(
                core.warningStyle
                + f"Retrying missing subject lineages individually: {len(missing_lineage_taxids)} taxid(s)"
                + core.normalStyle,
                flush=True,
            )
        for taxid in missing_lineage_taxids:
            limiter.acquire()
            taxid_to_lineage[taxid] = _fetch_lineage_from_taxid(taxid, api_key)

    accession_candidates = {}
    for region_i, (accession, sstart, send) in enumerate(distinct_regions, start=1):
        full_sequence = accession_to_sequence.get(accession)
        if full_sequence in [None, ""] and "." in accession:
            full_sequence = accession_to_sequence.get(accession.split(".", 1)[0])

        if full_sequence not in [None, ""]:
            seq_start = min(sstart, send)
            seq_stop = max(sstart, send)
            if seq_start >= 1 and seq_stop <= len(full_sequence):
                seq = full_sequence[seq_start - 1:seq_stop]
                if db_prefix == "n" and sstart > send:
                    seq = _reverse_complement_dna(seq)
                seq = seq.upper()
            else:
                seq = ""
            if seq != "":
                taxid = accession_to_taxid.get(accession)
                lineage = taxid_to_lineage.get(taxid, {rank: "" for rank in MIAN_RANKS})

                # Skip sequences with no class-level (or finer) resolution: these
                # are poorly characterised submissions that cannot be separated by
                # marker in downstream analyses.
                if (lineage.get("class") or "").strip() == "":
                    continue

                accession_entries = accession_candidates.setdefault(accession, {})
                if seq not in accession_entries:
                    accession_entries[seq] = {
                        "seq_id": accession,
                        "sequence": seq,
                        "lineage": lineage,
                        "accession": accession,
                        "sstart": sstart,
                        "send": send,
                        "duplicates": 1,
                        "first_region_index": region_i,
                    }
                else:
                    accession_entries[seq]["duplicates"] += 1

    seq_entries = []
    seq_to_index = {}
    representative_entries = []
    for accession_entries in accession_candidates.values():
        representative_entries.append(max(accession_entries.values(), key=_subject_interval_representative_key))

    representative_entries.sort(key=lambda entry: entry["first_region_index"])
    for entry in representative_entries:
        seq = entry["sequence"]
        if seq not in seq_to_index:
            seq_to_index[seq] = len(seq_entries)
            seq_entries.append(entry)
        else:
            seq_entries[seq_to_index[seq]]["duplicates"] += entry["duplicates"]

    taxonomy_output_path = str(Path(output_path).with_name(Path(output_path).stem + "_taxonomy.tsv"))
    with open(output_path, "w", encoding="utf-8") as outfile:
        for entry in seq_entries:
            header = (
                f">{entry['seq_id']} accession={entry['accession']} "
                f"sstart={entry['sstart']} send={entry['send']} duplicates={entry['duplicates']}"
            )
            outfile.write(header + "\n")
            outfile.write(entry["sequence"] + "\n")

    with open(taxonomy_output_path, "w", encoding="utf-8", newline="") as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=["Accession"] + MIAN_RANKS, delimiter="\t")
        writer.writeheader()
        for entry in seq_entries:
            lineage = entry.get("lineage", {rank: "" for rank in MIAN_RANKS})
            row = {"Accession": entry["seq_id"]}
            row.update({rank: lineage.get(rank, "") for rank in MIAN_RANKS})
            writer.writerow(row)

    return len(distinct_regions), len(seq_entries), taxonomy_output_path


def _run_live_blastn_and_save_hit_table(core):
    fasta_path_obj = Path(core.metaXplorFasta)
    if not fasta_path_obj.is_file():
        print(
            core.errorStyle
            + f"Input FASTA file not found: {core.metaXplorFasta}. Please run step 4a first."
            + core.normalStyle
        )
        return None, None
    if fasta_path_obj.stat().st_size == 0:
        print(
            core.errorStyle
            + f"Input FASTA file is empty: {core.metaXplorFasta}. Please run step 4a first."
            + core.normalStyle
        )
        return None, None

    fasta_text = fasta_path_obj.read_text(encoding="utf-8", errors="replace")
    if not fasta_text.strip().startswith(">"):
        print(core.errorStyle + "Input file does not look like FASTA (missing leading '>')." + core.normalStyle)
        return None, None

    db = "core_nt"
    print(core.warningStyle + "Querying BLAST database: core_nt" + core.normalStyle)

    evalue = None
    while evalue is None:
        raw = input(core.promptStyle + "BLAST e-value threshold" + core.normalStyle + " (default = 1e-5): ").strip()
        if raw == "":
            raw = "1e-5"
        try:
            evalue = float(raw)
        except ValueError:
            print(core.errorStyle + f"Invalid e-value: {raw}" + core.normalStyle)

    max_hits = None
    while max_hits is None:
        raw = input(core.promptStyle + "Maximum hits per query" + core.normalStyle + " (default = 10): ").strip()
        if raw == "":
            raw = "10"
        if raw.isnumeric() and int(raw) >= 1:
            max_hits = int(raw)
        else:
            print(core.errorStyle + f"Invalid max hits: {raw}" + core.normalStyle)

    api_key = input(core.promptStyle + "NCBI API key (optional)" + core.normalStyle + " (press Enter to skip): ").strip()
    if api_key == "":
        api_key = os.environ.get("NCBI_API_KEY", "")

    rid, rtoe = _blast_submit(core, fasta_text, db, evalue, max_hits, api_key)
    if rid is None:
        return None, None
    if not _blast_wait_for_results(core, rid, rtoe):
        return None, None

    prefix = fasta_path_obj.stem + "_ncbi_blastn"
    # out_xml = Path(core.current_dir) / "outputs" / f"{prefix}.xml"
    out_tab = Path(core.current_dir) / "outputs" / f"{prefix}_hittable.txt"

    # print(core.warningStyle + "Downloading XML results..." + core.normalStyle)
    # try:
    #     out_xml.write_bytes(_blast_retrieve(rid, "XML"))
    # except RuntimeError as e:
    #     print(core.errorStyle + str(e) + core.normalStyle)
    #     print(core.warningStyle + "You can retry live BLAST or proceed with your own hit-table file." + core.normalStyle)
    #     return None
    # time.sleep(3)
    print(core.warningStyle + "Downloading hit table (tabular) results..." + core.normalStyle)
    tab_bytes = _download_blast_tabular(core, rid)
    if tab_bytes is None:
        print(core.warningStyle + "You can retry live BLAST or proceed with your own hit-table file." + core.normalStyle)
        return None, None
    out_tab.write_bytes(tab_bytes)

    # print(core.successStyle + f"Saved BLAST XML to {out_xml}" + core.normalStyle)
    print(core.successStyle + f"Saved BLAST hit table to {out_tab}" + core.normalStyle)
    return str(out_tab), max_hits


def _rotate_tsv(input_path, output_path):
    with open(input_path, "r", encoding="utf-8") as infile:
        rows = [line.rstrip("\n").split("\t") for line in infile if line.strip() != ""]

    if len(rows) == 0:
        raise ValueError("Input sequence-composition file is empty")

    max_cols = max(len(row) for row in rows)
    padded_rows = [row + [""] * (max_cols - len(row)) for row in rows]
    rotated_rows = list(zip(*padded_rows))

    with open(output_path, "w", encoding="utf-8", newline="") as outfile:
        for row in rotated_rows:
            outfile.write("\t".join(row).rstrip("\t") + "\n")


def _parse_accession_from_sseqid(sseqid):
    return sseqid.split(":", 1)[1] if ":" in sseqid else sseqid


def _load_all_hits(assignments_path):
    """Returns {qseqid: [accession, ...]} keeping only hits that share the best bitscore."""
    # First pass: collect all rows and track the best bitscore per query.
    rows_by_query = {}
    best_bitscore = {}
    with open(assignments_path, "r", encoding="utf-8", newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        for row in reader:
            qseqid = (row.get("qseqid") or "").strip()
            sseqid = (row.get("sseqid") or "").strip()
            if qseqid == "" or sseqid == "":
                continue
            try:
                bitscore = float((row.get("bitscore") or "0").strip())
            except ValueError:
                bitscore = 0.0
            accession = _parse_accession_from_sseqid(sseqid)
            if qseqid not in rows_by_query:
                rows_by_query[qseqid] = []
                best_bitscore[qseqid] = bitscore
            else:
                if bitscore > best_bitscore[qseqid]:
                    best_bitscore[qseqid] = bitscore
            rows_by_query[qseqid].append((accession, bitscore))

    # Second pass: keep only accessions at the best bitscore level.
    hits = {}
    for qseqid, entries in rows_by_query.items():
        top = best_bitscore[qseqid]
        seen = []
        for accession, bitscore in entries:
            if bitscore >= top and accession not in seen:
                seen.append(accession)
        hits[qseqid] = seen
    return hits


def _compute_lca(lineages):
    """Given a list of lineage dicts, return the LCA lineage.

    Walks MIAN_RANKS from broad to specific. Once a rank diverges across hits,
    that rank and all finer ranks are set to empty string.
    """
    result = {}
    diverged = False
    for rank in MIAN_RANKS:
        if diverged:
            result[rank] = ""
            continue
        non_empty = {lg.get(rank, "") for lg in lineages if lg.get(rank, "") != ""}
        if len(non_empty) == 1:
            result[rank] = non_empty.pop()
        elif len(non_empty) == 0:
            result[rank] = ""
        else:
            result[rank] = ""
            diverged = True
    return result


def _is_non_informative_taxon(name):
    value = (name or "").strip().lower()
    if value == "":
        return True
    if value in NON_INFORMATIVE_TAXA:
        return True
    return value.startswith("uncultured ") or value.startswith("unclassified ")


def _is_informative_lineage(lineage):
    for rank in MIAN_RANKS:
        value = (lineage.get(rank) or "").strip()
        if value != "" and not _is_non_informative_taxon(value):
            return True
    return False


def _lineage_score(lineage):
    informative = 0
    non_empty = 0
    for rank in MIAN_RANKS:
        value = (lineage.get(rank) or "").strip()
        if value != "":
            non_empty += 1
            if not _is_non_informative_taxon(value):
                informative += 1
    return (informative, non_empty)


def _select_best_lineage(lineages):
    best = {rank: "" for rank in MIAN_RANKS}
    best_score = (-1, -1)
    for lineage in lineages:
        score = _lineage_score(lineage)
        if score > best_score:
            best = lineage
            best_score = score
    return best


def _chunked(values, size):
    for i in range(0, len(values), size):
        yield values[i:i + size]


def _subject_interval_representative_key(entry):
    return (
        entry["duplicates"],
        len(entry["sequence"]),
        -min(entry["sstart"], entry["send"]),
        -max(entry["sstart"], entry["send"]),
        -entry["first_region_index"],
    )


# --------------------------------------------------------------------------- #
# Token-bucket rate limiter (shared across all worker threads)
# --------------------------------------------------------------------------- #

class _RateLimiter:
    """
    Thread-safe token bucket.  Calling .acquire() blocks until a request
    slot is available, ensuring we never exceed `rate` calls per second
    regardless of how many threads are running.
    """
    def __init__(self, rate: float):
        self._rate      = rate
        self._tokens    = rate          # start full
        self._lock      = threading.Lock()
        self._last_time = time.monotonic()

    def acquire(self):
        with self._lock:
            now             = time.monotonic()
            elapsed         = now - self._last_time
            self._last_time = now
            self._tokens    = min(self._rate, self._tokens + elapsed * self._rate)
            if self._tokens < 1:
                sleep_for    = (1 - self._tokens) / self._rate
                time.sleep(sleep_for)
                self._tokens = 0
            else:
                self._tokens -= 1


# --------------------------------------------------------------------------- #
# Internal eutils helpers
# --------------------------------------------------------------------------- #

def _build_url(endpoint: str, params: dict) -> str:
    return (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        + endpoint
        + "?"
        + urllib.parse.urlencode(params)
    )


def _http_get(url: str, timeout: int = 60) -> bytes | None:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            return resp.read()
    except Exception:
        return None


def _esummary_taxids_batch(accessions: list[str], api_key) -> dict[str, str | None]:
    """
    Fast path: nuccore esummary returns TaxId directly without parsing
    full GenBank XML.  Resolves the vast majority of standard accessions
    at a fraction of the bandwidth cost of efetch rettype=gb.
    """
    result: dict[str, str | None] = {a: None for a in accessions}
    params = {
        "db":      "nuccore",
        "id":      ",".join(accessions),
        "retmode": "xml",
        "tool":    "mbctools",
    }
    if api_key:
        params["api_key"] = api_key

    raw = _http_get(_build_url("esummary.fcgi", params))
    if raw is None:
        return result

    try:
        root = ET.fromstring(raw)
    except ET.ParseError:
        return result

    for docsum in root.findall("DocSum"):
        accession_version = ""
        taxid             = ""
        for item in docsum.findall("Item"):
            name = item.get("Name", "")
            if name == "AccessionVersion" and item.text:
                accession_version = item.text.strip()
            if name == "TaxId" and item.text:
                taxid = item.text.strip()

        if accession_version and taxid:
            result[accession_version] = taxid
            bare = accession_version.split(".")[0]
            if bare in result:
                result[bare] = taxid

    return result


def _efetch_taxids_batch(accessions: list[str], api_key) -> dict[str, str | None]:
    """
    Slow path / fallback: full GenBank XML fetch.  Used only for accessions
    that esummary did not resolve (suppressed, non-standard, or WGS records).
    Logic identical to the original _fetch_taxids_from_accessions inner loop.
    """
    result: dict[str, str | None] = {a: None for a in accessions}
    params = {
        "db":      "nuccore",
        "id":      ",".join(accessions),
        "rettype": "gb",
        "retmode": "xml",
        "tool":    "mbctools",
    }
    if api_key:
        params["api_key"] = api_key

    raw = _http_get(_build_url("efetch.fcgi", params))
    if raw is None:
        return result

    try:
        root = ET.fromstring(raw)
    except ET.ParseError:
        return result

    for gbseq in root.iter("GBSeq"):
        accession         = (gbseq.findtext("GBSeq_primary-accession") or "").strip()
        accession_version = (gbseq.findtext("GBSeq_accession-version")  or "").strip()
        target_accession  = accession_version or accession
        taxid = None

        for feature in gbseq.iter("GBFeature"):
            key = feature.find("GBFeature_key")
            if key is None or key.text != "source":
                continue
            for qualifier in feature.iter("GBQualifier"):
                qname  = qualifier.find("GBQualifier_name")
                qvalue = qualifier.find("GBQualifier_value")
                if (
                    qname  is not None and qname.text  == "db_xref"
                    and qvalue is not None and qvalue.text
                    and "taxon:" in qvalue.text
                ):
                    taxid = qvalue.text.split("taxon:", 1)[1].strip()
                    break
            if taxid is not None:
                break

        for key in (target_accession, accession):
            if key and key in result:
                result[key] = taxid

    return result


def _fetch_taxid_from_accession(accession, api_key):
    r = _esummary_taxids_batch([accession], api_key)
    if r.get(accession):
        return r[accession]
    r = _efetch_taxids_batch([accession], api_key)
    return r.get(accession)


def _fetch_taxids_for_batch(accessions, api_key, limiter):
    limiter.acquire()
    result = _esummary_taxids_batch(accessions, api_key)
    missing = [accession for accession in accessions if not result.get(accession)]
    if len(missing) > 0:
        limiter.acquire()
        result.update(_efetch_taxids_batch(missing, api_key))
    return result


def _fetch_taxids_from_accessions(accessions, api_key, limiter=None, progress_callback=None):
    """
    Concurrent accession -> taxid retrieval across batches.
    Each batch uses esummary first, then efetch fallback for unresolved accessions.
    """
    taxids = {accession: None for accession in accessions}
    if len(accessions) == 0:
        return taxids

    shared_limiter = limiter if limiter is not None else _RateLimiter(_ncbi_rate_limit(api_key))
    batches = list(_chunked(accessions, BATCH_SIZE))
    total_batches = len(batches)

    with ThreadPoolExecutor(max_workers=_ncbi_max_workers(api_key)) as pool:
        futures = {pool.submit(_fetch_taxids_for_batch, batch, api_key, shared_limiter): batch for batch in batches}
        completed = 0
        for future in as_completed(futures):
            taxids.update(future.result())
            completed += 1
            if progress_callback is not None:
                progress_callback(completed, total_batches)

    return taxids


def _fetch_lineage_from_taxid(taxid, api_key):
    result = {rank: "" for rank in MIAN_RANKS}
    params = {
        "db":      "taxonomy",
        "id":      taxid,
        "retmode": "xml",
        "tool":    "mbctools",
    }
    if api_key:
        params["api_key"] = api_key

    raw = _http_get(_build_url("efetch.fcgi", params), timeout=20)
    if raw is None:
        return result

    try:
        root = ET.fromstring(raw)
    except ET.ParseError:
        return result

    top_taxon = root.find("Taxon")
    if top_taxon is not None:
        result = _extract_lineage_from_taxon(top_taxon)

    return result


def _extract_lineage_from_taxon(taxon):
    result        = {rank: "" for rank in MIAN_RANKS}
    current_rank  = (taxon.findtext("Rank")          or "").lower()
    organism_name =  taxon.findtext("ScientificName") or ""

    lineage_ex = taxon.find("LineageEx")
    if lineage_ex is not None:
        for lineage_taxon in lineage_ex.findall("Taxon"):
            rank = (lineage_taxon.findtext("Rank") or "").lower()
            if rank in ("superkingdom", "domain"):
                rank = "kingdom"
            if rank in MIAN_RANKS:
                result[rank] = lineage_taxon.findtext("ScientificName") or ""

    rank = current_rank
    if rank in ("superkingdom", "domain"):
        rank = "kingdom"
    if rank in MIAN_RANKS:
        result[rank] = organism_name

    if result["species"] == "" and current_rank == "species" and organism_name != "":
        parts = organism_name.split()
        result["species"] = parts[1] if len(parts) > 1 else organism_name

    return result


def _fetch_lineages_for_batch(taxids, api_key, limiter):
    limiter.acquire()
    result = {taxid: {rank: "" for rank in MIAN_RANKS} for taxid in taxids}
    params = {
        "db": "taxonomy",
        "id": ",".join(taxids),
        "retmode": "xml",
        "tool": "mbctools",
    }
    if api_key:
        params["api_key"] = api_key

    raw = _http_get(_build_url("efetch.fcgi", params))
    if raw is None:
        return result

    try:
        root = ET.fromstring(raw)
    except ET.ParseError:
        return result

    for taxon in root.findall("Taxon"):
        taxid = (taxon.findtext("TaxId") or "").strip()
        if taxid:
            result[taxid] = _extract_lineage_from_taxon(taxon)

    return result


def _fetch_lineages_from_taxids(taxids, api_key, limiter=None, progress_callback=None):
    """
    Concurrent taxid -> lineage retrieval across batches.
    """
    lineage_map = {taxid: {rank: "" for rank in MIAN_RANKS} for taxid in taxids}
    if len(taxids) == 0:
        return lineage_map

    shared_limiter = limiter if limiter is not None else _RateLimiter(_ncbi_rate_limit(api_key))
    batches = list(_chunked(taxids, BATCH_SIZE))
    total_batches = len(batches)

    with ThreadPoolExecutor(max_workers=_ncbi_max_workers(api_key)) as pool:
        futures = {pool.submit(_fetch_lineages_for_batch, batch, api_key, shared_limiter): batch for batch in batches}
        completed = 0
        for future in as_completed(futures):
            lineage_map.update(future.result())
            completed += 1
            if progress_callback is not None:
                progress_callback(completed, total_batches)

    return lineage_map


def _build_mian_taxonomy(core, assignments_path, output_path, api_key):
    """
    - time.sleep() between requests removed — rate limiting is now handled
      by the token-bucket limiter shared across worker threads.
    - A single _RateLimiter instance is created here and reused across
      both the taxid and lineage fetch phases so the rate budget is
      shared correctly across the whole function.
    - Progress reporting preserved exactly as in the original.
    """
    all_hits = _load_all_hits(assignments_path)
    cache    = {}
    total    = len(all_hits)

    if total == 0:
        raise ValueError("No query hits found in assignments file")

    print(core.warningStyle + f"Building MIAN taxonomy for {total} query sequences (LCA across top-scoring hits)..." + core.normalStyle, flush=True)

    unique_accessions = []
    for accessions in all_hits.values():
        for accession in accessions:
            if accession not in cache and accession not in unique_accessions:
                unique_accessions.append(accession)

    if len(unique_accessions) > 0:
        print(core.warningStyle + "Fetching NCBI taxids for matched subjects..." + core.normalStyle, flush=True)

    # Single limiter shared across all network phases of this function
    rate    = _ncbi_rate_limit(api_key)
    limiter = _RateLimiter(rate)

    # ── accession → taxid ─────────────────────────────────────────────────
    accession_to_taxid = _fetch_taxids_from_accessions(
        unique_accessions,
        api_key,
        limiter=limiter,
        progress_callback=lambda completed, total_batches: print(
            core.warningStyle
            + f"Taxid batch progress: {completed}/{total_batches}"
            + core.normalStyle,
            flush=True,
        ),
    )

    missing_accessions = [a for a in unique_accessions if accession_to_taxid.get(a) in (None, "")]
    if len(missing_accessions) > 0:
        print(
            core.warningStyle
            + f"Retrying taxid retrieval for {len(missing_accessions)} accession(s) individually..."
            + core.normalStyle,
            flush=True,
        )
        for accession in missing_accessions:
            limiter.acquire()
            accession_to_taxid[accession] = _fetch_taxid_from_accession(accession, api_key)

    # ── taxid → lineage ───────────────────────────────────────────────────
    unique_taxids = []
    for accession in unique_accessions:
        taxid = accession_to_taxid.get(accession)
        if taxid is not None and taxid not in unique_taxids:
            unique_taxids.append(taxid)

    taxid_to_lineage = {}

    if len(unique_taxids) > 0:
        print(core.warningStyle + f"Fetching NCBI taxonomy lineages in batches of {BATCH_SIZE}..." + core.normalStyle, flush=True)

    taxid_to_lineage.update(
        _fetch_lineages_from_taxids(
            unique_taxids,
            api_key,
            limiter=limiter,
            progress_callback=lambda completed, total_batches: print(
                core.warningStyle
                + f"Lineage batch progress: {completed}/{total_batches}"
                + core.normalStyle,
                flush=True,
            ),
        )
    )

    missing_taxids = [
        t for t in unique_taxids
        if not _has_any_taxonomy_value(taxid_to_lineage.get(t, {rank: "" for rank in MIAN_RANKS}))
    ]
    if len(missing_taxids) > 0:
        print(
            core.warningStyle
            + f"Retrying lineage retrieval for {len(missing_taxids)} taxid(s) individually..."
            + core.normalStyle,
            flush=True,
        )
        for taxid in missing_taxids:
            limiter.acquire()
            taxid_to_lineage[taxid] = _fetch_lineage_from_taxid(taxid, api_key)

    # ── populate cache ─────────────────────────────────────────────────────
    for accession in unique_accessions:
        taxid = accession_to_taxid.get(accession)
        if taxid is None:
            cache[accession] = {rank: "" for rank in MIAN_RANKS}
        else:
            cache[accession] = taxid_to_lineage.get(taxid, {rank: "" for rank in MIAN_RANKS})

    # ── write TSV ──────────────────────────────────────────────────────────
    with open(output_path, "w", encoding="utf-8", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=["OTU"] + MIAN_RANKS, delimiter="\t")
        writer.writeheader()

        for i, (qseqid, accessions) in enumerate(all_hits.items(), start=1):
            lineages = [cache.get(accession, {rank: "" for rank in MIAN_RANKS}) for accession in accessions]

            informative_lineages = [lineage for lineage in lineages if _is_informative_lineage(lineage)]
            source_lineages      = informative_lineages if len(informative_lineages) > 0 else lineages
            lca                  = _compute_lca(source_lineages)

            if not _has_any_taxonomy_value(lca):
                non_empty_lineages = [lineage for lineage in source_lineages if _has_any_taxonomy_value(lineage)]
                if len(non_empty_lineages) > 0:
                    lca = _select_best_lineage(non_empty_lineages)

            row = {"OTU": qseqid}
            row.update({rank: lca.get(rank, "") for rank in MIAN_RANKS})
            writer.writerow(row)

            if i == 1 or i == total or i % max(1, total // 20) == 0:
                percent = int((i * 100) / total)
                print(core.warningStyle + f"Taxonomy file building progress: {i}/{total} ({percent}%)" + core.normalStyle, flush=True)

    print(core.warningStyle + "Taxonomy progress: done (100%)" + core.normalStyle, flush=True)


def menu4e(core):
    prepare_with_previous_params(core)

    if not os.path.isfile(core.metaXplorSequenceComposition) or os.path.getsize(core.metaXplorSequenceComposition) == 0:
        print(core.errorStyle + "\nFile " + core.metaXplorSequenceComposition + " is missing or empty. Please run step 4a" + core.normalStyle)
        core.rerun(core.main_menu4)
    if not os.path.isfile(core.metaXplorAssignments) or os.path.getsize(core.metaXplorAssignments) == 0:
        print(core.errorStyle + "\nFile " + core.metaXplorAssignments + " is missing or empty. Please run step 4b" + core.normalStyle)
        core.rerun(core.main_menu4)
    if not os.path.isfile(core.metaXplorSamples) or os.path.getsize(core.metaXplorSamples) == 0:
        print(core.errorStyle + "\nFile " + core.metaXplorSamples + " is missing or empty. Please run step 4c" + core.normalStyle)
        core.rerun(core.main_menu4)

    mian_sequences_path = str(Path(core.current_dir) / "tmp_files" / "mian_OTU.tsv")
    mian_taxonomy_path = str(Path(core.current_dir) / "tmp_files" / "mian_taxonomy.tsv")

    try:
        _rotate_tsv(core.metaXplorSequenceComposition, mian_sequences_path)
        print(core.successStyle + "File " + mian_sequences_path + " was successfully written" + core.normalStyle)
    except Exception as e:
        print(core.errorStyle + "Unable to rotate sequence file: " + str(e) + core.normalStyle)
        core.rerun(core.main_menu4)

    api_key = input(core.promptStyle + "NCBI API key (optional)" + core.normalStyle + " (press Enter to skip): ").strip()
    if api_key == "":
        api_key = os.environ.get("NCBI_API_KEY", "")

    try:
        _build_mian_taxonomy(core, core.metaXplorAssignments, mian_taxonomy_path, api_key)
        print(core.successStyle + "File " + mian_taxonomy_path + " was successfully written" + core.normalStyle)
    except Exception as e:
        print(core.errorStyle + "Unable to build taxonomy file: " + str(e) + core.normalStyle)
        core.rerun(core.main_menu4)

    mian_zip_name = "mbctools_mian_export_" + core.date.strftime("%Y%m%d") + ".zip"
    mian_zip_path = str(Path(core.current_dir) / mian_zip_name)
    with zipfile.ZipFile(mian_zip_path, mode="w") as zf:
        zf.write(mian_sequences_path, "mian_OTU.tsv")
        zf.write(mian_taxonomy_path, "mian_taxonomy.tsv")
        zf.write(core.metaXplorSamples, "mian_sample_metadata.tsv")
    print(core.successStyle + "MIAN archive was successfully created as " + mian_zip_path + core.normalStyle)

    core.rerun(core.main_menu4)


def main_menu4(core):
    """Displays submenu 4."""
    core.os.system("cls" if core.winOS else "clear")
    print(
        core.titleStyle
        + "\n--- MENU 4: EXPORTING ANALYSIS RESULTS INTO metaXplor FORMAT ---"
        + core.normalStyle
        + "\n\n"
        "4a -> Generate sequence files\n"
        "\tCompiles all sequences selected for all loci into a single fasta\n"
        "\tOutputs a .tsv file indicating samples weights for each sequence\n\n"
        "4b -> Generate assignment file\n"
        "\tConverts blastn results (obtained from blasting above-mentioned fasta file) from 'Hit table (text)'"
        "\n\t(format #7) into metaXplor format\n\n"
        "4c -> Build metaXplor-format sample metadata file from provided tabulated file\n\n"
        "4d -> Compress all metaXplor files into a final, ready to import, zip archive\n\n"
        "4e -> Generate MIAN data files from metaXplor files\n"
        + core.normalStyle
    )

    core.rmenu = core.promptUser(
        "Please select an option among those listed above",
        None,
        ["4a", "4b", "4c", "4d", "4e", "back", "home", "exit"],
        1,
        core.main,
        "",
    )

    if core.rmenu == "4a":
        menu4a(core)
    elif core.rmenu == "4b":
        menu4b(core)
    elif core.rmenu == "4c":
        menu4c(core)
    elif core.rmenu == "4d":
        menu4d(core, True)
    elif core.rmenu == "4e":
        menu4e(core)


def menu4a(core):
    prepare_with_previous_params(core)
    locusToFastaDict = {}
    for locus in list(set(core.lociPEs + core.lociSEs)):
        files2cat = core.glob.glob("results_by_locus/" + locus + core.fileSep + "*_allseq_select.fasta")
        if len(files2cat) > 0:
            locusToFastaDict[locus] = "results_by_locus/" + locus + core.fileSep + locus + "_allseq_select.fasta"

    print()
    logFile = open(f"{core.current_dir}{core.fileSep}outputs{core.fileSep}res4a.log", "w")
    derepResults = core.derepSeveralFastaFiles(locusToFastaDict, core.metaXplorFasta, core.metaXplorSequenceComposition, logFile)
    logFile.close()

    if derepResults[0] == 0:
        print(core.errorStyle + "\nNo concatenated sequences found for any loci. Please run step 3 before retrying" + core.normalStyle)
        os.remove(core.metaXplorFasta)
        os.remove(core.metaXplorSequenceComposition)
        core.rerun(core.main_menu3)
    else:
        print(core.successStyle + "\n" + str(derepResults[0]) + " distinct sequences were compiled into .fasta and .tsv files")
        if len(derepResults[1]) > 0:
            print(
                core.warningStyle
                + "Warning: not all sequences could be included because concatenation step (#3) was not run on some loci: "
                + ", ".join(derepResults[1])
                + core.successStyle
            )
        print(
            "You may now run blastn on "
            + core.metaXplorFasta
            + ", then launch step 4b to either run BLAST live from within mbctools or provide your own 'Hit table (text)' (format #7)"
            + core.normalStyle
        )
        print(
            core.warningStyle
            + "\nNB: The appropriate format required by step 4b may be obtained by either:"
            + core.normalStyle
            + "\n - if using NCBI online BLAST interface, selecting 'Hit table (text)' from the 'Download All' dropdown list"
            + "\n - if executing command line BLAST, specifying the following argument: -outfmt 7\n"
            + "(Its header must contain a line starting with: "
            + core.citationStyle
            + "'# Fields: query acc.ver, subject acc.ver, '"
            + core.normalStyle
            + ")"
        )
        core.rerun(core.main_menu4)


def menu4b(core):
    prepare_with_previous_params(core)

    print(
        core.warningStyle
        + "\nNB: The appropriate format required by step 4b may be obtained by either:"
        + core.normalStyle
        + "\n - if using NCBI online BLAST interface, selecting 'Hit table (text)' from the 'Download All' dropdown list"
        + "\n - if executing command line BLAST, specifying the following argument: -outfmt 7\n"
        + "(Its header must contain a line starting with: "
        + core.citationStyle
        + "'# Fields: query acc.ver, subject acc.ver, '"
        + core.normalStyle
        + ")"
    )

    generatedHitTable = None
    liveBlastMaxHits = None
    runLiveBlast = core.promptUser(
        "Do you want to run a guided remote NCBI BLASTn now? Enter yes or no",
        None,
        ["yes", "no", "back", "home", "exit"],
        1,
        core.main_menu4,
        "",
    )
    if runLiveBlast == "yes":
        generatedHitTable, liveBlastMaxHits = _run_live_blastn_and_save_hit_table(core)
        if generatedHitTable is None:
            print(core.errorStyle + "Could not generate hit table from live BLAST." + core.normalStyle)
            menu4b(core)
            return

    blastTextHitTable = core.promptUser(
        "Enter path to blastn hit-table (text format #7)",
        generatedHitTable,
        ["back", "home", "exit"],
        3,
        core.main_menu4,
        "",
    )
    with open(blastTextHitTable.strip(), "r") as infile:
        lines = re.sub(r"\s\s+", "\t", infile.read()).splitlines()

    if len(lines) == 0:
        print(core.errorStyle + "Provided hit-table file is empty!" + core.normalStyle)
        menu4b(core)

    database = None
    previousQseqId = None
    blastType = lines[0].split(" ")[1].strip()

    defaultMaxHits = str(liveBlastMaxHits) if liveBlastMaxHits is not None else "5"
    maxHits = None
    print()
    while maxHits is None or not maxHits.isnumeric() or int(maxHits) < 1:
        if maxHits is not None:
            print(core.errorStyle + "\n--> WRONG INPUT: " + maxHits + core.normalStyle)
        maxHits = input(
            core.promptStyle
            + "Enter maximum number of retained hits per query."
            + core.normalStyle
            + " Default is "
            + defaultMaxHits
            + ": "
        )
        if maxHits == "":
            maxHits = defaultMaxHits
    maxHits = int(maxHits)

    print()
    with open(core.metaXplorAssignments, "w") as outfile:
        i = 0
        nHitsForQseqId = 0
        while i < len(lines):
            if lines[i].startswith("#"):
                if database is None and "Database:" in lines[i]:
                    try:
                        database = lines[i].split(" ")[2]
                        nucleotide_databases = {
                            "nt",
                            "core_nt",
                            "refseq_rna",
                            "refseq_genomic",
                            "est",
                            "gss",
                            "pat",
                            "sts",
                            "htgs",
                            "env_nt",
                            "16S_rRNA",
                            "tsa_nr",
                            "wgs",
                            "mitogenomes",
                            "plastid_genomes",
                        }
                        protein_databases = {
                            "nr",
                            "refseq_protein",
                            "swissprot",
                            "pdb",
                            "env_nr",
                            "pat_protein",
                            "uniprotkb",
                            "cdd",
                        }

                        if database in nucleotide_databases:
                            database = "n"
                        elif database in protein_databases:
                            database = "p"
                        else:
                            raise Exception(f"Unsupported Database type: {database}")

                    except Exception as e:
                        print(core.errorStyle + f"Unable to parse accession type prefix in '{lines[i]}': {e}" + core.normalStyle)
                        os.remove(core.metaXplorAssignments)
                        menu4b(core)
                elif previousQseqId is None and "Fields:" in lines[i]:
                    outfile.write(
                        lines[i]
                        .split(":")[1]
                        .strip()
                        .replace(", ", "\t")
                        .replace("query acc.ver", "qseqid")
                        .replace("subject acc.ver", "sseqid")
                        .replace("bit score", "bitscore")
                        .replace("% identity", "pident")
                        .replace("q. start", "qstart")
                        .replace("q. end", "qend")
                        + "\tassignment_method\tbest_hit\n"
                    )
            else:
                if database is None:
                    print(core.errorStyle + "Unable to determine accession type prefix" + core.normalStyle)
                    os.remove(core.metaXplorAssignments)
                    menu4b(core)

                splitLine = lines[i].strip().split("\t")
                newQuery = previousQseqId != splitLine[0]
                if newQuery:
                    nHitsForQseqId = 0

                if nHitsForQseqId < maxHits:
                    j = 0
                    while j < len(splitLine):
                        if j > 0:
                            outfile.write("\t")
                        if j == 1:
                            outfile.write(database + ":")
                        outfile.write(splitLine[j])
                        j += 1

                    outfile.write("\t")
                    outfile.write(blastType)
                    outfile.write("\t")

                    if newQuery:
                        outfile.write("Y")
                        previousQseqId = splitLine[0]

                    outfile.write("\n")

                nHitsForQseqId += 1
            i += 1

    print(core.successStyle + "File " + core.metaXplorAssignments + " was successfully written" + core.normalStyle)

    exportMatchedSubjectFasta = core.promptUser(
        "Also export FASTA for aligned subject intervals (sstart/send) from this hit-table? Enter yes or no",
        None,
        ["yes", "no"],
        1,
        None,
        "",
    )
    if exportMatchedSubjectFasta == "yes":
        api_key = input(core.promptStyle + "NCBI API key (optional)" + core.normalStyle + " (press Enter to skip): ").strip()
        if api_key == "":
            api_key = os.environ.get("NCBI_API_KEY", "")

        output_fasta_path = str(Path(core.current_dir) / "outputs" / (Path(blastTextHitTable).stem + "_subject_intervals.fasta"))
        try:
            requested, retrieved, output_taxonomy_path = _download_subject_aligned_regions_fasta(
                core, lines, database, output_fasta_path, api_key
            )
            if retrieved == 0:
                print(core.warningStyle + "No aligned subject interval sequences could be retrieved from NCBI." + core.normalStyle)
            elif retrieved < requested:
                print(
                    core.warningStyle
                    + f"Retrieved {retrieved}/{requested} deduplicated aligned subject interval sequences into "
                    + output_fasta_path
                    + " and wrote taxonomy into "
                    + output_taxonomy_path
                    + core.normalStyle
                )
            else:
                print(
                    core.successStyle
                    + f"Retrieved {retrieved} deduplicated aligned subject interval sequences into "
                    + output_fasta_path
                    + " and wrote taxonomy into "
                    + output_taxonomy_path
                    + core.normalStyle
                )
        except Exception as e:
            print(core.errorStyle + "Unable to export subject FASTA: " + str(e) + core.normalStyle)

    menu4d(core, False)


def menu4c(core):
    prepare_with_previous_params(core)

    print(
        f"\n\nYou must now provide a tabulated metadata file for your samples. {core.warningStyle}A header field named 'Sample' is expected for the column featuring sample names{core.normalStyle}"
    )
    print(
        f"Any field with '{core.warningStyle}date{core.normalStyle}' in its header name will be considered to be the sample {core.warningStyle}collection date{core.normalStyle} (recommended format: YYYY-MM-DD)"
    )
    print(f"{core.warningStyle}Collection location{core.normalStyle} may be specified:")
    print("\t- either as commma-separated decimal-format values in a single field named 'LatLon' (e.g. -17.7127, -67.9905)")
    print("\t- or in separate columns named 'Latitude' and 'Longitude', in decimal format (e.g. -17.7127) or DMS format (e.g. 16°42'45.6\"S)")
    print("Any additional columns will remain named as provided")

    sampleMetadataFile = core.promptUser("Enter path to tabulated sample metadata file", None, ["back", "home", "exit"], 3, core.main_menu4, "")
    with open(sampleMetadataFile.strip(), "r") as infile:
        lines = infile.read().splitlines()

    if len(lines) == 0:
        print(core.errorStyle + "Provided sample metadata file is empty!" + core.normalStyle)
        menu4c(core)

    headerCols = re.sub(r"\s\s+", "\t", lines[0]).split("\t")
    sampleIndex = None
    collectionDateIndex = None
    latitudeIndex = None
    longitudeIndex = None
    latLonIndex = None

    i = 0
    while i < len(headerCols):
        headerCol = headerCols[i].lower()
        if "sample" in headerCol:
            if sampleIndex is not None:
                print(core.errorStyle + "Ambiguity identifying sample name column between '" + headerCols[sampleIndex] + "'' and '" + headerCols[i] + "'" + core.normalStyle)
                menu4c(core)
            sampleIndex = i
        elif "date" in headerCol:
            if collectionDateIndex is not None:
                print(core.errorStyle + "Ambiguity identifying collection date column between '" + headerCols[collectionDateIndex] + "'' and '" + headerCols[i] + "'" + core.normalStyle)
                menu4c(core)
            collectionDateIndex = i
        elif headerCol.startswith("lat"):
            if "lon" in headerCol:
                if latLonIndex is not None:
                    print(core.errorStyle + "Ambiguity identifying LatLon column between '" + headerCols[latLonIndex] + "'' and '" + headerCols[i] + "'" + core.normalStyle)
                    menu4c(core)
                latLonIndex = i
            else:
                if latitudeIndex is not None:
                    print(core.errorStyle + "Ambiguity identifying latitude column between '" + headerCols[latitudeIndex] + "'' and '" + headerCols[i] + "'" + core.normalStyle)
                    menu4c(core)
                latitudeIndex = i
        elif headerCol.startswith("lon"):
            if "lat" in headerCol:
                if latLonIndex is not None:
                    print(core.errorStyle + "Ambiguity identifying LatLon column between '" + headerCols[latLonIndex] + "'' and '" + headerCols[i] + "'" + core.normalStyle)
                    menu4c(core)
                latLonIndex = i
            else:
                if longitudeIndex is not None:
                    print(core.errorStyle + "Ambiguity identifying longitude column between '" + headerCols[longitudeIndex] + "'' and '" + headerCols[i] + "'" + core.normalStyle)
                    menu4c(core)
                longitudeIndex = i
        i += 1

    if sampleIndex is None:
        print(core.errorStyle + "Unable to identify sample name column! Please read instructions carefully and submit a corrected file" + core.normalStyle)
        menu4c(core)

    if latLonIndex is None and (not None not in [latitudeIndex, longitudeIndex]):
        print(core.warningStyle + "Unable to identify latitude and/or longitude column(s)! Generated file will contain empty values for this field" + core.normalStyle)
    if collectionDateIndex is None:
        print(core.warningStyle + "Unable to identify collection date column! Generated file will contain empty values for this field" + core.normalStyle)

    specialIndexes = [sampleIndex, latitudeIndex, longitudeIndex, latLonIndex, collectionDateIndex]
    with open(core.metaXplorSamples, "w") as outfile:
        j = 0
        outfile.write("sample_name\tlat_lon\tcollection_date")
        while j < len(headerCols):
            if j not in specialIndexes:
                outfile.write("\t" + headerCols[j])
            j += 1
        outfile.write("\n")

        i = 1
        samplesToProcess = core.samples.copy()
        while i < len(lines):
            splitLine = re.sub("  +", "\t", lines[i]).split("\t")

            if splitLine[sampleIndex] not in core.samples:
                print(core.warningStyle + "Skipping unknown sample: " + splitLine[sampleIndex] + core.normalStyle)
            else:
                while len(splitLine) < len(headerCols):
                    splitLine.append("")

                outfile.write(
                    splitLine[sampleIndex]
                    + "\t"
                    + determine_lat_lon(core, splitLine, latLonIndex, latitudeIndex, longitudeIndex, sampleIndex)
                    + "\t"
                )
                collDate = splitLine[collectionDateIndex].strip() if collectionDateIndex is not None else None
                if collDate is not None and collDate.strip() != "":
                    try:
                        outfile.write(str(parse_date(collDate)).replace(" 00:00:00", ""))
                    except Exception as e:
                        print(core.errorStyle + "Unable to parse date at line " + str(i) + " \"" + collDate + "\": " + str(e) + core.normalStyle)
                        menu4c(core)

                j = 0
                while j < len(splitLine):
                    if j not in specialIndexes:
                        outfile.write("\t" + splitLine[j])
                    j += 1
                outfile.write("\n")
                samplesToProcess.remove(splitLine[sampleIndex])
            i += 1

    if len(samplesToProcess) > 0:
        print(core.errorStyle + "Provided file lacks lines for the following sample(s): " + ", ".join(samplesToProcess) + core.normalStyle)
        os.remove(core.metaXplorSamples)
    else:
        print(core.successStyle + "File " + core.metaXplorSamples + " was successfully written" + core.normalStyle)
    menu4d(core, False)


def menu4d(core, invokedByUser):
    if not os.path.isfile(core.metaXplorFasta) or os.path.getsize(core.metaXplorFasta) == 0:
        if invokedByUser:
            print(core.errorStyle + "\nFile " + core.metaXplorFasta + " is missing or empty. Please run step 4a" + core.normalStyle)
        core.rerun(core.main_menu4)
    if not os.path.isfile(core.metaXplorSequenceComposition) or os.path.getsize(core.metaXplorSequenceComposition) == 0:
        if invokedByUser:
            print(core.errorStyle + "\nFile " + core.metaXplorSequenceComposition + " is missing or empty. Please run step 4a" + core.normalStyle)
        core.rerun(core.main_menu4)
    if not os.path.isfile(core.metaXplorAssignments) or os.path.getsize(core.metaXplorAssignments) == 0:
        if invokedByUser:
            print(core.errorStyle + "\nFile " + core.metaXplorAssignments + " is missing or empty. Please run step 4b" + core.normalStyle)
        core.rerun(core.main_menu4)
    if not os.path.isfile(core.metaXplorSamples) or os.path.getsize(core.metaXplorSamples) == 0:
        if invokedByUser:
            print(core.errorStyle + "\nFile " + core.metaXplorSamples + " is missing or empty. Please run step 4c" + core.normalStyle)
        core.rerun(core.main_menu4)

    print("\n\nAll metaXplor files seem to be ready for zipping.")
    core.zipNow = core.promptUser(
        "Zip them now to create the final import file? " + core.normalStyle + "Enter yes or no" + core.promptStyle,
        None,
        ["yes", "no"],
        1,
        None,
        "",
    )
    if core.zipNow == "yes":
        b = io.BytesIO()
        zf = zipfile.ZipFile(b, mode="w")
        zf.write(core.metaXplorSamples, os.path.basename(core.metaXplorSamples))
        zf.write(core.metaXplorAssignments, os.path.basename(core.metaXplorAssignments))
        zf.write(core.metaXplorFasta, os.path.basename(core.metaXplorFasta))
        zf.write(core.metaXplorSequenceComposition, os.path.basename(core.metaXplorSequenceComposition))
        zf.close()
        zipFileName = "mbctools_metaXplor_export_" + core.date.strftime("%Y%m%d") + ".zip"
        open(zipFileName, "wb").write(b.getbuffer())
        print(core.successStyle + "\n\nmetaXplor import archive was successfully created as " + core.current_dir + core.fileSep + zipFileName + core.normalStyle)
        core.rerun(core.main_menu4)
    else:
        main_menu4(core)


def parse_date(date_str):
    match = re.match(r"(\d{4}[/-]?[0-9]{1,2}[/-]?[0-9]{1,2})", date_str)
    if match:
        separator = match.group(0)[4] if len(match.group(0)) > 4 else "-"
        date_parts = date_str.split(separator)
        year = int(date_parts[0])
        month = int(date_parts[1])
        day = int(date_parts[2])
        return datetime.datetime(year, month, day)
    return datetime.datetime.strptime(date_str, "%B %d, %Y")


def dms_to_decimal(dmsString):
    deg, minutes, seconds, direction = re.split('[°\'\"]', dmsString.replace("''", '"').replace(" ", ""))
    return round((float(deg) + float(minutes) / 60 + float(seconds) / (60 * 60)) * (-1 if direction.upper() in ["W", "S"] else 1), 6)


def determine_lat_lon(core, cellArray, latLonIndex, latitudeIndex, longitudeIndex, sampleIndex):
    gotSeparateLatAndLong = None not in [latitudeIndex, longitudeIndex]
    if latLonIndex is not None:
        splitCoords = cellArray[latLonIndex].replace(";", ",").split(",")
        try:
            return str(round(float(splitCoords[0]), 6)) + ", " + str(round(float(splitCoords[1]), 6))
        except Exception:
            pass

    latLon = ""
    if gotSeparateLatAndLong:
        try:
            latLon += str(round(float(cellArray[latitudeIndex]), 6))
        except Exception:
            try:
                latLon += str(dms_to_decimal(cellArray[latitudeIndex]))
            except Exception:
                latLon = ""
        if latLon != "":
            try:
                latLon += ", " + str(round(float(cellArray[longitudeIndex]), 6))
            except Exception:
                try:
                    latLon += ", " + str(dms_to_decimal(cellArray[longitudeIndex]))
                except Exception:
                    latLon = ""

    if latLon == "":
        msg = "Sample " + cellArray[sampleIndex] + ":"
        if latLonIndex is not None:
            msg += " Unable to parse LatLon field" + ("" if gotSeparateLatAndLong else (" '" + cellArray[latLonIndex] + "'")) + "."
        if gotSeparateLatAndLong:
            msg += " Unable to parse latitude / longitude fields."
        print(core.warningStyle + msg + core.normalStyle)
    return latLon
