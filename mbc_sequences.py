import glob
import os
import re
import shutil
import subprocess


def concat_sequences_by_locus(
    loci,
    lociPEs,
    lociSEs,
    samples,
    current_dir,
    fileSep,
    promptUser,
    main,
    derep_several_fasta_files,
    warningStyle,
    errorStyle,
    successStyle,
    normalStyle,
    promptStyle,
):
    """Compiles selected sample sequences by locus into one file and optionally dereplicates."""
    nb_samples = len(samples)
    all_loci = lociPEs + list(set(lociSEs) - set(lociPEs))
    if os.path.exists("outputs/Stats_option_3.txt"):
        os.remove("outputs/Stats_option_3.txt")

    skippedLoci = []
    locusIndex = 0
    anySeqsFound = False
    while loci is None or locusIndex < len(loci):
        loc2cat = (
            promptUser(
                "For which LOCUS do you want to generate a unique sequence file?",
                None,
                all_loci + ["back", "home", "exit"],
                1,
                main,
                f"Results of concatenation session --> {current_dir}{fileSep}outputs{fileSep}Stats_option_3.txt",
            )
            if loci is None
            else loci[locusIndex]
        )

        stat_3 = open(f"{current_dir}{fileSep}outputs{fileSep}Stats_option_3.txt", "a")
        os.chdir(f"{current_dir}{fileSep}results_by_locus{fileSep}{loc2cat}")
        files2cat = glob.glob("*_select.fas")
        if len(files2cat) == 0:
            print(
                errorStyle
                + f"\nSequences for locus {loc2cat} have not been filtered using an abundance threshold."
                + warningStyle
                + " Please run step 2 on all loci for which you want to run step 3"
                + normalStyle
            )
            locusIndex = locusIndex + 1
            skippedLoci.append(loc2cat)
        else:
            with open(f"./{loc2cat}_allseq_select.fasta", "w") as out:
                for file in files2cat:
                    if os.path.exists(file):
                        with open(file, "r") as out2:
                            out.write(out2.read())
            tot = open("./" + loc2cat + "_allseq_select.fasta")
            nb_tot = tot.read().count(">")
            if nb_tot > 0:
                anySeqsFound = True
            stat_3.writelines(f"Locus {loc2cat} has {nb_tot} sequences in total\n")
            print(
                successStyle
                + f"\nLocus {loc2cat}: {nb_tot} sequences from {nb_samples} samples were added to a unique fasta file\n"
                + normalStyle
            )
            if nb_tot > 0:
                print(
                    f"Results in --> {current_dir}{fileSep}results_by_locus{fileSep}{loc2cat}{fileSep}{loc2cat}_allseq_select.fasta\n"
                )

            if loci is not None:
                locusIndex = locusIndex + 1
            elif (
                nb_tot > 0
                and "yes"
                == promptUser(
                    "Do you want to generate dereplicated versions of the above mentioned fasta file? "
                    + normalStyle
                    + "Enter yes or no"
                    + promptStyle,
                    None,
                    ["yes", "no"],
                    1,
                    None,
                    "",
                )
            ):
                print()
                logFile = open(f"{current_dir}{fileSep}outputs{fileSep}res3.log", "w")
                outFasta = (
                    f"{current_dir}{fileSep}results_by_locus{fileSep}{loc2cat}{fileSep}{loc2cat}"
                    + "_allseq_select_derep.fasta"
                )
                outTsv = (
                    f"{current_dir}{fileSep}results_by_locus{fileSep}{loc2cat}{fileSep}{loc2cat}"
                    + "_allseq_select_derep.tsv"
                )
                derepResults = derep_several_fasta_files(
                    {loc2cat: loc2cat + "_allseq_select.fasta"},
                    loc2cat + "_allseq_select_derep.fasta",
                    loc2cat + "_allseq_select_derep.tsv",
                    logFile,
                    samples,
                    warningStyle,
                    normalStyle,
                    fileSep,
                )
                print(
                    successStyle
                    + str(derepResults[0])
                    + " distinct sequences were dereplicated into fasta and tsv files: "
                    + normalStyle
                    + outFasta
                    + ", "
                    + outTsv
                    + "\n"
                )
                logFile.close()
        stat_3.close()

    if loci is not None and anySeqsFound is True and "yes" == promptUser(
        "Do you want to generate dereplicated versions of the above mentioned fasta files? "
        + normalStyle
        + "Enter yes or no"
        + promptStyle,
        None,
        ["yes", "no"],
        1,
        None,
        "",
    ):
        print()
        logFile = open(f"{current_dir}{fileSep}outputs{fileSep}res3.log", "w")
        for locus in loci:
            if locus not in skippedLoci:
                outFasta = (
                    f"{current_dir}{fileSep}results_by_locus{fileSep}{locus}{fileSep}{locus}"
                    + "_allseq_select_derep.fasta"
                )
                outTsv = (
                    f"{current_dir}{fileSep}results_by_locus{fileSep}{locus}{fileSep}{locus}"
                    + "_allseq_select_derep.tsv"
                )
                derepResults = derep_several_fasta_files(
                    {
                        locus: f"{current_dir}{fileSep}results_by_locus{fileSep}{locus}{fileSep}{locus}"
                        + "_allseq_select.fasta"
                    },
                    outFasta,
                    outTsv,
                    logFile,
                    samples,
                    warningStyle,
                    normalStyle,
                    fileSep,
                )
                if derepResults[0] > 0:
                    print(
                        successStyle
                        + str(derepResults[0])
                        + " distinct sequences were dereplicated into fasta and tsv files: "
                        + normalStyle
                        + outFasta
                        + ", "
                        + outTsv
                        + normalStyle
                        + "\n"
                    )
                else:
                    print(warningStyle + "No sequences to process" + normalStyle + "\n")
        logFile.close()

    input("\nPress ENTER to continue ")


def derep_based_on_ids(
    locusToFastaDict,
    outFastaName,
    outTsvName,
    samples,
    warningStyle,
    normalStyle,
    fileSep,
):
    readTypePattern = r"(merged|R1[\-+])"
    fastaContents = ""
    with open(outFastaName, "w") as fastaFile, open(outTsvName, "w") as seqCompositionFile:
        i = 0
        seqCompositionFile.write("qseqid")
        while i < len(samples):
            seqCompositionFile.write("\t" + samples[i])
            i += 1

        seqSampleAbundances = {}
        seqHashToSamplesDict = {}
        seqOrientations = {}
        skippedLoci = []
        totalDistinctSeqCount = 0
        for locus in locusToFastaDict:
            print("Dereplicating sample sequences by locus for " + locus + "...")
            try:
                file1 = open(locusToFastaDict[locus], "r")
                lines = file1.readlines()

                j = 0
                skipActualSequence = False
                while j < len(lines):
                    line = lines[j].replace("\r\n", "").replace("\n", "")
                    if line.startswith(">"):
                        splitIdLine = line[1:].split(" ")
                        sampleAndCount = re.sub(
                            r"_(merged|R1).*;", ";", splitIdLine[1].replace("sample=", "")
                        ).split(";size=")

                        if splitIdLine[0] not in seqSampleAbundances:
                            fastaContents += ">" + splitIdLine[0] + "\n"
                            seqSampleAbundances[splitIdLine[0]] = ["0"] * len(samples)
                            totalDistinctSeqCount += 1
                            seqOrientations[splitIdLine[0]] = re.search(
                                readTypePattern, splitIdLine[1]
                            ).group(0)
                            skipActualSequence = False
                        else:
                            skipActualSequence = True
                        seqSampleAbundances[splitIdLine[0]][samples.index(sampleAndCount[0])] = sampleAndCount[1]
                        if splitIdLine[0] not in seqHashToSamplesDict:
                            seqHashToSamplesDict[splitIdLine[0]] = []
                        seqHashToSamplesDict[splitIdLine[0]].append(sampleAndCount[0])
                    elif skipActualSequence is False:
                        fastaContents += lines[j]
                    j += 1
            except FileNotFoundError:
                print(
                    warningStyle
                    + "File not found: results_by_locus/"
                    + locus
                    + fileSep
                    + locus
                    + "_allseq_select.fasta: skipping locus "
                    + locus
                    + normalStyle
                )
                skippedLoci.append(locus)
            i += 1

        i = 1
        padLevel = len(str(len(seqSampleAbundances)))
        for seqId in seqSampleAbundances:
            seqName = (
                "seq"
                + str(i).zfill(padLevel)
                + "-"
                + seqOrientations[seqId]
                + "."
                + "_".join(seqHashToSamplesDict[seqId][0:5])
                + ("" if len(seqHashToSamplesDict[seqId]) <= 5 else "...")
            )
            fastaContents = fastaContents.replace(">" + seqId + "\n", ">" + seqName + "\n")
            seqCompositionFile.write("\n" + seqName + "\t" + "\t".join(seqSampleAbundances[seqId]))
            i += 1

        fastaFile.write(fastaContents)

    return [totalDistinctSeqCount, skippedLoci]


def derep_several_fasta_files(
    locusToFastaDict,
    outFastaName,
    outTsvName,
    logFile,
    samples,
    warningStyle,
    normalStyle,
    fileSep,
):
    with open(outFastaName + ".tmp1", "wb") as wfd:
        for sp in locusToFastaDict:
            with open(locusToFastaDict[sp], "rb") as fd:
                shutil.copyfileobj(fd, wfd)

    subprocess.run(
        [
            "vsearch",
            "--fastx_filter",
            outFastaName + ".tmp1",
            "--relabel_sha1",
            "--relabel_keep",
            "--fastaout",
            outFastaName + ".tmp2",
        ],
        stderr=logFile,
    )
    os.remove(outFastaName + ".tmp1")
    retVal = derep_based_on_ids(
        {", ".join(locusToFastaDict.keys()): outFastaName + ".tmp2"},
        outFastaName,
        outTsvName,
        samples,
        warningStyle,
        normalStyle,
        fileSep,
    )
    os.remove(outFastaName + ".tmp2")

    return retVal
