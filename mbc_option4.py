from mbc_option_common import prepare_with_previous_params
import datetime
import csv
import io
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path


BLAST_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
INITIAL_WAIT = 15
POLL_INTERVAL = 10
MAX_WAIT = 3600
MIAN_RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


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


def _run_live_blastn_and_save_hit_table(core):
    fasta_path_obj = Path(core.metaXplorFasta)
    if not fasta_path_obj.is_file():
        print(
            core.errorStyle
            + f"Input FASTA file not found: {core.metaXplorFasta}. Please run step 4a first."
            + core.normalStyle
        )
        return None
    if fasta_path_obj.stat().st_size == 0:
        print(
            core.errorStyle
            + f"Input FASTA file is empty: {core.metaXplorFasta}. Please run step 4a first."
            + core.normalStyle
        )
        return None

    fasta_text = fasta_path_obj.read_text(encoding="utf-8", errors="replace")
    if not fasta_text.strip().startswith(">"):
        print(core.errorStyle + "Input file does not look like FASTA (missing leading '>')." + core.normalStyle)
        return None

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
        return None
    if not _blast_wait_for_results(core, rid, rtoe):
        return None

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
        return None
    out_tab.write_bytes(tab_bytes)

    # print(core.successStyle + f"Saved BLAST XML to {out_xml}" + core.normalStyle)
    print(core.successStyle + f"Saved BLAST hit table to {out_tab}" + core.normalStyle)
    return str(out_tab)


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


def _fetch_taxid_from_accession(accession, api_key):
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "gb",
        "retmode": "xml",
        "tool": "mbctools",
    }
    if api_key:
        params["api_key"] = api_key

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url, timeout=20) as resp:
            root = ET.fromstring(resp.read())
        for feature in root.iter("GBFeature"):
            key = feature.find("GBFeature_key")
            if key is not None and key.text == "source":
                for qualifier in feature.iter("GBQualifier"):
                    qname = qualifier.find("GBQualifier_name")
                    qvalue = qualifier.find("GBQualifier_value")
                    if qname is not None and qname.text == "db_xref" and qvalue is not None and qvalue.text and "taxon:" in qvalue.text:
                        return qvalue.text.split("taxon:", 1)[1].strip()
    except Exception:
        return None
    return None


def _fetch_lineage_from_taxid(taxid, api_key):
    result = {rank: "" for rank in MIAN_RANKS}
    params = {
        "db": "taxonomy",
        "id": taxid,
        "retmode": "xml",
        "tool": "mbctools",
    }
    if api_key:
        params["api_key"] = api_key

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url, timeout=20) as resp:
            root = ET.fromstring(resp.read())

        rank_elem = root.find(".//Taxon/Rank")
        current_rank = (rank_elem.text or "").lower() if rank_elem is not None else ""
        org_elem = root.find(".//Taxon/ScientificName")
        organism_name = org_elem.text or "" if org_elem is not None else ""

        for taxon in root.iter("Taxon"):
            rank_tag = taxon.find("Rank")
            name_tag = taxon.find("ScientificName")
            if rank_tag is None or name_tag is None:
                continue
            rank = (rank_tag.text or "").lower()
            if rank == "superkingdom":
                rank = "kingdom"
            if rank in MIAN_RANKS:
                result[rank] = name_tag.text or ""

        if result["species"] == "" and current_rank == "species" and organism_name != "":
            parts = organism_name.split()
            result["species"] = parts[1] if len(parts) > 1 else organism_name
    except Exception:
        pass

    return result


def _build_mian_taxonomy(core, assignments_path, output_path, api_key):
    all_hits = _load_all_hits(assignments_path)
    cache = {}
    sleep_between_requests = 0.11 if api_key else 0.34
    total = len(all_hits)

    if total == 0:
        raise ValueError("No query hits found in assignments file")

    print(core.warningStyle + f"Building MIAN taxonomy for {total} query sequences (LCA across top-scoring hits)...")

    with open(output_path, "w", encoding="utf-8", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=["OTU"] + MIAN_RANKS, delimiter="\t")
        writer.writeheader()

        for i, (qseqid, accessions) in enumerate(all_hits.items(), start=1):
            lineages = []
            for accession in accessions:
                if accession in cache:
                    lineages.append(cache[accession])
                else:
                    time.sleep(sleep_between_requests)
                    taxid = _fetch_taxid_from_accession(accession, api_key)
                    if taxid is None:
                        lineage = {rank: "" for rank in MIAN_RANKS}
                    else:
                        time.sleep(sleep_between_requests)
                        lineage = _fetch_lineage_from_taxid(taxid, api_key)
                    cache[accession] = lineage
                    lineages.append(lineage)

            lca = _compute_lca(lineages)
            row = {"OTU": qseqid}
            row.update({rank: lca.get(rank, "") for rank in MIAN_RANKS})
            writer.writerow(row)

            # Print visible progress updates in terminals that do not render carriage-return rewrites.
            if i == 1 or i == total or i % max(1, total // 20) == 0:
                percent = int((i * 100) / total)
                print(f"Taxonomy file building progress: {i}/{total} ({percent}%)", flush=True)

    print("Taxonomy progress: done (100%)" + core.normalStyle, flush=True)


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
    runLiveBlast = core.promptUser(
        "Do you want to run a guided remote NCBI BLASTn now? Enter yes or no",
        None,
        ["yes", "no", "back", "home", "exit"],
        1,
        core.main_menu4,
        "",
    )
    if runLiveBlast == "yes":
        generatedHitTable = _run_live_blastn_and_save_hit_table(core)
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

    maxHits = None
    print()
    while maxHits is None or not maxHits.isnumeric() or int(maxHits) < 1:
        if maxHits is not None:
            print(core.errorStyle + "\n--> WRONG INPUT: " + maxHits + core.normalStyle)
        maxHits = input(core.promptStyle + "Enter maximum number of retained hits per query." + core.normalStyle + " Default is 5: ")
        if maxHits == "":
            maxHits = "5"
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
        core.rerun(None)
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
