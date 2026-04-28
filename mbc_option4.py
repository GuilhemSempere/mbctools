from mbc_option_common import prepare_with_previous_params
import datetime
import io
import os
import re
import zipfile


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
        "4d -> Compresses all metaXplor files into a final, ready to import, zip archive\n"
        + core.normalStyle
    )

    core.rmenu = core.promptUser(
        "Please select an option among those listed above",
        None,
        ["4a", "4b", "4c", "4d", "back", "home", "exit"],
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
            + ", download all results as 'Hit table (text)' (format #7), then come back and launch step 4b"
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

    blastTextHitTable = core.promptUser("Enter path to blastn hit-table (text format #7)", None, ["back", "home", "exit"], 3, core.main_menu4, "")
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
