from mbc_option_common import prepare_with_previous_params, show_missing_loci_and_return
import math
import os
import re
import subprocess


def _run_trim_option(core, has_loci, no_loci_message, stats_to_reset=None):
    prepare_with_previous_params(core)
    if not has_loci:
        show_missing_loci_and_return(core, no_loci_message, main_menu2)
        return

    if stats_to_reset is not None and core.os.path.exists(stats_to_reset):
        core.os.remove(stats_to_reset)

    core.trim_2x()
    core.rerun(None)


def main_menu2(core):
    """Displays submenu 2."""
    core.os.system("cls" if core.winOS else "clear")
    print(
        core.titleStyle
        + "\n--- MENU 2: PRIMER REMOVAL, SELECTION OF MINIMUM SEQUENCE ABUNDANCE LEVELS ACCORDING TO USER-DEFINED THRESHOLDS ---"
        + core.normalStyle
        + "\n\n"
        "2a -> Apply the SAME size threshold for ALL SAMPLES for the loci based on PAIRED-END reads "
        "(R1/R2 merged)\n"
        "\ti.e. you want to keep only sequences whose abundance is greater than x% of the total number of"
        "\n\tsequences for a given sample. This threshold of x% can be chosen for each locus.\n\n"
        "2b -> Apply the SAME size threshold for ALL SAMPLES for the loci based on SINGLE-END reads "
        "(R1 only)\n"
        "\tsame as option 2a but only using the R1 reads instead of merged ones.\n\n"
        "2c -> Apply a SPECIFIC size threshold for a given sample, for the loci based on PAIRED-END reads "
        "(R1/R2 merged)\n"
        "\ti.e. you want to modulate the threshold of x% by locus but also by sample within a particular locus.\n\n"
        "2d -> Apply a SPECIFIC size threshold for a given sample, for the loci based on SINGLE-END reads "
        "(R1 only)\n"
        "\tsame as option 2c but only using the R1 sequences instead of merged ones.\n"
        + core.normalStyle
    )

    core.rmenu = core.promptUser(
        "Please select an option among those listed above",
        None,
        ["2a", "2b", "2c", "2d", "back", "home", "exit"],
        1,
        core.main,
        "",
    )

    if core.rmenu == "2a":
        menu2a(core)
    elif core.rmenu == "2b":
        menu2b(core)
    elif core.rmenu == "2c":
        menu2c(core)
    elif core.rmenu == "2d":
        menu2d(core)


def menu2a(core):
    _run_trim_option(
        core,
        has_loci=(len(core.lociPEs) > 0),
        no_loci_message="No paired-end reads in current selection!",
    )


def menu2b(core):
    _run_trim_option(
        core,
        has_loci=(len(core.lociSEs) > 0),
        no_loci_message="No R1/single-end reads in current selection!",
    )


def menu2c(core):
    _run_trim_option(
        core,
        has_loci=(len(core.lociPEs) > 0),
        no_loci_message="No paired-end reads in current selection!",
        stats_to_reset="outputs/Stats_option_2c.txt",
    )


def menu2d(core):
    _run_trim_option(
        core,
        has_loci=(len(core.lociSEs) > 0),
        no_loci_message="No R1/single-end reads in current selection!",
        stats_to_reset="outputs/Stats_option_2d.txt",
    )


def get_single_seq_orient_file_suffix(core, loci, samples, rmenu):
    for locus in loci:
        nTotalPlusOrientedReadCount = 0
        nTotalMinusOrientedReadCount = 0
        nTotalPlusOrientedClusterCount = 0
        nTotalMinusOrientedClusterCount = 0
        for sample in samples:
            plusOrientedClusters = []
            minusOrientedClusters = []
            with open(f"./results_by_locus/{locus}/{sample}_singleEnd_orient.tsv", "r") as orientInfo:
                for line in orientInfo.read().split("\n"):
                    splitLine = line.split("\t")
                    if len(splitLine) == 4:
                        if splitLine[1] == "+":
                            nTotalPlusOrientedReadCount = nTotalPlusOrientedReadCount + int(re.sub(r".*;size=", "", splitLine[0]))
                            plusOrientedClusters.append(splitLine[0])
                        elif splitLine[1] == "-":
                            nTotalMinusOrientedReadCount = nTotalMinusOrientedReadCount + int(re.sub(r".*;size=", "", splitLine[0]))
                            minusOrientedClusters.append(splitLine[0])

            nTotalPlusOrientedClusterCount = nTotalPlusOrientedClusterCount + len(plusOrientedClusters)
            nTotalMinusOrientedClusterCount = nTotalMinusOrientedClusterCount + len(minusOrientedClusters)
            nTotalSeqCount = len(plusOrientedClusters) + len(minusOrientedClusters)
            if nTotalSeqCount > 0:
                if len(plusOrientedClusters) > 0:
                    with open(f"./results_by_locus/{locus}/{sample}_singleEnd_orient_plus.tsv", "w") as orientInfo:
                        for seq in plusOrientedClusters:
                            orientInfo.write(seq + "\n")
                    core.dos2unix(f"./results_by_locus/{locus}/{sample}_singleEnd_orient_plus.tsv")

                if len(minusOrientedClusters) > 0:
                    with open(f"./results_by_locus/{locus}/{sample}_singleEnd_orient_minus.tsv", "w") as orientInfo:
                        for seq in minusOrientedClusters:
                            orientInfo.write(seq + "\n")
                    core.dos2unix(f"./results_by_locus/{locus}/{sample}_singleEnd_orient_minus.tsv")

        if nTotalPlusOrientedReadCount + nTotalMinusOrientedReadCount == 0:
            return None
        print(
            f"\nFor locus {locus}, found "
            + str(nTotalPlusOrientedReadCount)
            + " sense (+) reads distributed in "
            + str(nTotalPlusOrientedClusterCount)
            + " clusters and "
            + str(nTotalMinusOrientedReadCount)
            + " antisense (-) reads distributed in "
            + str(nTotalMinusOrientedClusterCount)
            + " clusters.\n"
            + core.warningStyle
            + "We discourage accounting for both, which would probably generate 2 separate alignments."
            + core.normalStyle
        )
        promptOptions = ["back", "home", "exit"]
        if nTotalPlusOrientedReadCount > 0 and nTotalMinusOrientedReadCount > 0:
            promptOptions = ["*"] + promptOptions
        if nTotalMinusOrientedReadCount > 0:
            promptOptions = ["-"] + promptOptions
        if nTotalPlusOrientedReadCount > 0:
            promptOptions = ["+"] + promptOptions
        strand = core.promptUser(
            "Enter \"+\" to keep only sense clusters, \"-\" to keep only antisense clusters, \"*\" to keep both",
            None,
            promptOptions,
            1,
            core.trim_2x,
            "",
        )
        if strand != "*":
            strand = "plus" if strand == "+" else "minus"
            for sample in samples:
                if os.path.exists(f"./results_by_locus/{locus}/{sample}_singleEnd_orient_" + strand + ".tsv"):
                    logFile = open(f"{core.current_dir}{core.fileSep}outputs{core.fileSep}{rmenu}_{locus}_{sample}_{strand}.log", "a")
                    subprocess.run(
                        [
                            "vsearch",
                            "--fastx_getseqs",
                            f"./results_by_locus/{locus}/{sample}_singleEnd_orient.fas",
                            "--labels",
                            f"./results_by_locus/{locus}/{sample}_singleEnd_orient_" + strand + ".tsv",
                            "--fastaout",
                            f"./results_by_locus/{locus}/{sample}_singleEnd_orient_" + strand + ".fas",
                        ],
                        stderr=logFile,
                    )
                    logFile.close()
            return "_" + strand
        return ""


def trim_2x(core):
    if core.rmenu == "2a":
        core.os.chdir(core.current_dir)
        stat_2a = open(f"{core.current_dir}{core.fileSep}outputs{core.fileSep}Stats_option_2a.txt", "w")
        while True:
            core.loc2trim2a = core.in_loc2trim_2x()
            core.os.chdir(f"./results_by_locus/{core.loc2trim2a}")
            core.trim_left = core.in_trim_left("")
            core.trim_right = core.in_trim_right("")
            ts = core.in_ts()
            stat_2a.write(
                f"Locus {core.loc2trim2a} trimmed {core.trim_left} bp (forward) and {core.trim_right} bp (reverse) with threshold set at {core.ts1}\n"
            )
            for sample in core.samples:
                with open(f"{sample}_pairedEnd_orient.fas", "r") as filin, open("trim-select." + core.scriptExt, "w") as out:
                    targets = [line for line in filin if "size" in line]
                    a = 0
                    for target in targets:
                        size = re.search("size=(.+?)$", target).group(1)
                        a = a + int(size)
                    b = math.ceil(a * float(core.ts1))
                    stat_2a.writelines(f"\tSum of cluster sizes for {sample} = {a}\n\tWith threshold {ts[0]}, sizes >= {b} for {sample} were retained\n")
                    out.write(
                        core.start_log_redirect("./" + core.loc2trim2a + "_" + sample + ".log")
                        + f"vsearch --fastx_filter {sample}_pairedEnd_orient.fas --fastq_stripleft {core.trim_left} --fastq_stripright {core.trim_right} --fastaout tmp\n"
                        + core.localErrorOnStopCmd
                        + "\n"
                        + f"vsearch --derep_fulllength tmp --output tmp2 --sizein --sizeout\n"
                        + core.localErrorOnStopCmd
                        + "\n"
                        + f"vsearch --fastx_filter tmp2 --minsize {b} --fastaout {sample}_pairedEnd_select.fas\n"
                        + core.localErrorOnStopCmd
                        + "\n"
                        + core.end_log_redirect("./" + sample + ".log")
                    )
                subprocess.run([core.shellCmd, "./trim-select." + core.scriptExt])
                if os.path.exists(f"{sample}_pairedEnd_select.fas") and os.path.getsize(f"{sample}_pairedEnd_select.fas") > 0:
                    selected = open(f"{sample}_pairedEnd_select.fas", "r")
                    nb_selected = selected.read().count(">")
                    stat_2a.writelines(f"\tNumber of selected clusters for sample {sample}: {nb_selected}\n\n")
                    core.sys.stdout.write(
                        f"Sum of cluster sizes for {sample} at locus {core.loc2trim2a} = {a}: with threshold set at {ts[0]}%,{core.warningStyle} clusters with size >= {b} were retained{core.successStyle}\n"
                        f"Number of selected clusters for sample {sample}: {nb_selected}\n"
                        + core.normalStyle
                    )
                else:
                    core.sys.stdout.write(core.warningStyle + f"No sequences could be selected for sample {sample} on selected locus\n" + core.normalStyle)
                if os.path.exists("tmp"):
                    os.remove("tmp")
                if os.path.exists("tmp2"):
                    os.remove("tmp2")
            if os.path.exists("trim-select." + core.scriptExt):
                os.remove("trim-select." + core.scriptExt)
            else:
                core.sys.stdout.write(core.warningStyle + "\nNo data found to process for selected locus\n" + core.normalStyle)
            core.os.chdir(core.current_dir)
            stat_2a.flush()

    if core.rmenu == "2b":
        core.os.chdir(core.current_dir)
        stat_2b = open(f"{core.current_dir}{core.fileSep}outputs{core.fileSep}Stats_option_2b.txt", "w")
        while True:
            core.loc2trim2b = core.in_loc2trim_2x()
            orientFileSuffix = get_single_seq_orient_file_suffix(core, [core.loc2trim2b], core.samples, core.rmenu)
            if orientFileSuffix is None:
                core.sys.stdout.write(core.errorStyle + f"\nNo reads found at locus {core.loc2trim2b} for any sample\n" + core.normalStyle)
            else:
                core.os.chdir(f"./results_by_locus/{core.loc2trim2b}")
                core.trim_left = core.in_trim_left(orientFileSuffix)
                core.trim_right = core.in_trim_right(orientFileSuffix)
                ts = core.in_ts()
                stat_2b.write(
                    f"Locus {core.loc2trim2b} trimmed {core.trim_left} bp (forward) and {core.trim_right} bp (reverse) with threshold set at {core.ts1}\n"
                )
                for sample in core.samples:
                    if os.path.exists(f"{sample}_singleEnd_orient{orientFileSuffix}.fas") and os.path.getsize(f"{sample}_singleEnd_orient{orientFileSuffix}.fas") > 0:
                        with open(f"{sample}_singleEnd_orient{orientFileSuffix}.fas", "r") as filin, open("trim-select." + core.scriptExt, "w") as out:
                            targets = [line for line in filin if "size" in line]
                            a = 0
                            for target in targets:
                                size = re.search("size=(.+?)$", target).group(1)
                                a = a + int(size)
                            b = math.ceil(a * float(core.ts1))
                            stat_2b.writelines(f"\tSum of cluster sizes for {sample} = {a}\n\tWith threshold {ts[0]}, sizes >= {b} for {sample} were retained\n")
                            out.write(
                                core.start_log_redirect("./" + core.loc2trim2b + "_" + sample + ".log")
                                + f"vsearch --fastx_filter {sample}_singleEnd_orient{orientFileSuffix}.fas --fastq_stripleft {core.trim_left} --fastq_stripright {core.trim_right} --fastaout tmp\n"
                                + core.localErrorOnStopCmd
                                + "\n"
                                + f"vsearch --derep_fulllength tmp --output tmp2 --sizein --sizeout\n"
                                + core.localErrorOnStopCmd
                                + "\n"
                                + "vsearch --fastx_filter tmp2 "
                                + (f"--minsize {b} " if b >= 1 else "")
                                + f"--fastaout {sample}_singleEnd_select.fas\n"
                                + core.localErrorOnStopCmd
                                + "\n"
                                + core.end_log_redirect("./" + sample + ".log")
                            )
                        subprocess.run([core.shellCmd, "./trim-select." + core.scriptExt])

                        senseReplacements = []
                        if os.path.isfile(sample + "_singleEnd_orient_plus.tsv") and os.path.getsize(sample + "_singleEnd_orient_plus.tsv") > 0:
                            for line in open(sample + "_singleEnd_orient_plus.tsv", "r"):
                                idWithoutSize = line.split(";size=")[0]
                                senseReplacements.append((idWithoutSize, idWithoutSize.replace("_R1.", "_R1+.")))
                        antiSenseReplacements = []
                        if os.path.isfile(sample + "_singleEnd_orient_minus.tsv") and os.path.getsize(sample + "_singleEnd_orient_minus.tsv") > 0:
                            for line in open(sample + "_singleEnd_orient_minus.tsv", "r"):
                                idWithoutSize = line.split(";size=")[0]
                                antiSenseReplacements.append((idWithoutSize, idWithoutSize.replace("_R1.", "_R1-.")))
                        core.replaceInFile(sample + "_singleEnd_orient_plus.tsv", senseReplacements)
                        core.replaceInFile(sample + "_singleEnd_orient_minus.tsv", antiSenseReplacements)
                        core.replaceInFile(sample + "_singleEnd_select.fas", senseReplacements + antiSenseReplacements)

                        nb_selected = (
                            open(sample + "_singleEnd_select.fas", "r").read().count(">")
                            if os.path.exists(sample + "_singleEnd_select.fas")
                            else 0
                        )
                        stat_2b.writelines(f"\tNumber of selected clusters for sample {sample}: {nb_selected}\n\n")
                        core.sys.stdout.write(
                            f"Sum of cluster sizes for {sample} at locus {core.loc2trim2b} = {a}: with threshold set at {ts[0]}%,{core.warningStyle} clusters with size >= {b} were retained{core.successStyle}\n"
                            f"Number of selected clusters for sample {sample}: {nb_selected}\n"
                            + core.normalStyle
                        )
                        if os.path.exists("tmp"):
                            os.remove("tmp")
                        if os.path.exists("tmp2"):
                            os.remove("tmp2")
                    else:
                        core.sys.stdout.write(core.warningStyle + f"No sequences could be selected for sample {sample} on selected locus\n" + core.normalStyle)
                if os.path.exists("trim-select." + core.scriptExt):
                    os.remove("trim-select." + core.scriptExt)
                else:
                    core.sys.stdout.write(core.warningStyle + "\nNo data found to process for selected locus\n" + core.normalStyle)
                core.os.chdir(core.current_dir)
                stat_2b.flush()

    if core.rmenu == "2c":
        core.os.chdir(core.current_dir)
        stat_2c = open(f"{core.current_dir}{core.fileSep}outputs{core.fileSep}Stats_option_2c.txt", "a")
        while True:
            core.loc2trim2c = core.in_loc2trim_2x()
            core.os.chdir(f"{core.current_dir}{core.fileSep}results_by_locus{core.fileSep}{core.loc2trim2c}")
            core.trim_left = core.in_trim_left("")
            core.trim_right = core.in_trim_right("")
            stat_2c.write(f"Locus {core.loc2trim2c} trimmed at {core.trim_left} bp (forward) and {core.trim_right} bp (reverse)\n")
            while True:
                core.sam2trim2c = core.in_trim_sample2c(core.loc2trim2c)
                if not os.path.exists(f"{core.sam2trim2c}_pairedEnd_orient.fas") or os.path.getsize(f"{core.sam2trim2c}_pairedEnd_orient.fas") == 0:
                    core.sys.stdout.write(core.errorStyle + f"\nNo reads found at locus {core.loc2trim2c} for sample {core.sam2trim2c}\n" + core.normalStyle)
                else:
                    ts = core.in_ts()
                    with open(f"{core.sam2trim2c}_pairedEnd_orient.fas", "r") as filin, open("trim-select." + core.scriptExt, "w") as filout:
                        targets = [line for line in filin if "size" in line]
                        a = 0
                        for target in targets:
                            size = re.search("size=(.+?)$", target).group(1)
                            a = a + int(size)
                        b = math.ceil(a * float(core.ts1))
                        stat_2c.write(
                            f"\tSum of cluster sizes for {core.sam2trim2c} at locus {core.loc2trim2c} = {a}\n\tWith threshold {ts[0]}, sizes >= {b} for {core.sam2trim2c} were retained\n"
                        )
                        filout.writelines(
                            core.start_log_redirect("./" + core.loc2trim2c + "_" + core.sam2trim2c + ".log")
                            + f"vsearch --fastx_filter {core.sam2trim2c}_pairedEnd_orient.fas --fastq_stripleft {core.trim_left} --fastq_stripright {core.trim_right} --fastaout tmp\n"
                            + core.localErrorOnStopCmd
                            + "\n"
                            + f"vsearch --derep_fulllength tmp --output tmp2 --sizein --sizeout\n"
                            + core.localErrorOnStopCmd
                            + "\n"
                            + f"vsearch --fastx_filter tmp2 --minsize {b} --fastaout {core.sam2trim2c}_pairedEnd_select.fas\n"
                            + core.localErrorOnStopCmd
                            + "\n"
                            + core.end_log_redirect("./" + core.sam2trim2c + ".log")
                        )
                    subprocess.run([core.shellCmd, "./trim-select." + core.scriptExt])
                    selected = open("./" + core.sam2trim2c + "_pairedEnd_select.fas", "r")
                    nb_selected = selected.read().count(">")
                    stat_2c.write(f"\tNumber of selected clusters for sample {core.sam2trim2c}: {nb_selected}\n\n")
                    core.sys.stdout.write(
                        f"Sum of cluster sizes for {core.sam2trim2c} at locus {core.loc2trim2c} = {a}: with threshold set at {ts[0]}%,{core.warningStyle} clusters with size >= {b} were retained{core.successStyle}\n"
                        f"Number of selected clusters for sample {core.sam2trim2c}: {nb_selected}\n"
                        + core.normalStyle
                    )
                    if os.path.exists("trim-select." + core.scriptExt):
                        os.remove("trim-select." + core.scriptExt)
                    else:
                        core.sys.stdout.write(core.warningStyle + "\nNo data found to process for selected locus\n" + core.normalStyle)
                    if os.path.exists("tmp"):
                        os.remove("tmp")
                    if os.path.exists("tmp2"):
                        os.remove("tmp2")
                    core.os.chdir(core.current_dir)
                    stat_2c.flush()

    if core.rmenu == "2d":
        core.os.chdir(core.current_dir)
        stat_2d = open(f"{core.current_dir}{core.fileSep}outputs{core.fileSep}Stats_option_2d.txt", "a")
        while True:
            core.loc2trim2d = core.in_loc2trim_2x()
            while True:
                core.sam2trim2d = core.in_trim_sample2d(core.loc2trim2d)
                orientFileSuffix = get_single_seq_orient_file_suffix(core, [core.loc2trim2d], [core.sam2trim2d], core.rmenu)
                if orientFileSuffix is None:
                    core.sys.stdout.write(core.errorStyle + f"\nNo reads found at locus {core.loc2trim2d} for sample {core.sam2trim2d}\n" + core.normalStyle)
                else:
                    core.os.chdir(f"{core.current_dir}{core.fileSep}results_by_locus{core.fileSep}{core.loc2trim2d}")
                    core.trim_left = core.in_trim_left(orientFileSuffix)
                    core.trim_right = core.in_trim_right(orientFileSuffix)
                    stat_2d.write(f"Locus {core.loc2trim2d} trimmed at {core.trim_left} bp (forward) and {core.trim_right} bp (reverse)\n")
                    ts = core.in_ts()
                    with open(core.sam2trim2d + f"_singleEnd_orient{orientFileSuffix}.fas", "r") as filin, open("trim-select." + core.scriptExt, "w") as filout:
                        targets = [line for line in filin if "size" in line]
                        a = 0
                        for target in targets:
                            size = re.search("size=(.+?)$", target).group(1)
                            a = a + int(size)
                        b = math.ceil(a * float(core.ts1))
                        stat_2d.write(
                            f"\tSum of cluster sizes for {core.sam2trim2d} at locus {core.loc2trim2d} = {a}\n\tWith threshold {ts[0]}, sizes >= {b} for {core.sam2trim2d} were retained\n"
                        )
                        filout.writelines(
                            core.start_log_redirect("./" + core.loc2trim2d + ".log")
                            + f"vsearch --fastx_filter {core.sam2trim2d}_singleEnd_orient{orientFileSuffix}.fas --fastq_stripleft {core.trim_left} --fastq_stripright {core.trim_right} --fastaout tmp\n"
                            + core.localErrorOnStopCmd
                            + "\n"
                            + f"vsearch --derep_fulllength tmp --output tmp2 --sizein --sizeout\n"
                            + core.localErrorOnStopCmd
                            + "\n"
                            + "vsearch --fastx_filter tmp2 "
                            + (f"--minsize {b} " if b >= 1 else "")
                            + f"--fastaout {core.sam2trim2d}_singleEnd_select.fas\n"
                            + core.localErrorOnStopCmd
                            + "\n"
                            + core.end_log_redirect("./" + core.loc2trim2d + "_" + core.sam2trim2d + ".log")
                        )
                    subprocess.run([core.shellCmd, "./trim-select." + core.scriptExt])

                    senseReplacements = []
                    if os.path.isfile(core.sam2trim2d + "_singleEnd_orient_plus.tsv") and os.path.getsize(core.sam2trim2d + "_singleEnd_orient_plus.tsv") > 0:
                        for line in open(core.sam2trim2d + "_singleEnd_orient_plus.tsv", "r"):
                            idWithoutSize = line.split(";size=")[0]
                            senseReplacements.append((idWithoutSize, idWithoutSize.replace("_R1.", "_R1+.")))
                    antiSenseReplacements = []
                    if os.path.isfile(core.sam2trim2d + "_singleEnd_orient_minus.tsv") and os.path.getsize(core.sam2trim2d + "_singleEnd_orient_minus.tsv") > 0:
                        for line in open(core.sam2trim2d + "_singleEnd_orient_minus.tsv", "r"):
                            idWithoutSize = line.split(";size=")[0]
                            antiSenseReplacements.append((idWithoutSize, idWithoutSize.replace("_R1.", "_R1-.")))
                    core.replaceInFile(core.sam2trim2d + "_singleEnd_orient_plus.tsv", senseReplacements)
                    core.replaceInFile(core.sam2trim2d + "_singleEnd_orient_minus.tsv", antiSenseReplacements)
                    core.replaceInFile(core.sam2trim2d + "_singleEnd_select.fas", senseReplacements + antiSenseReplacements)

                    selected = open("./" + core.sam2trim2d + "_singleEnd_select.fas", "r")
                    nb_selected = selected.read().count(">")
                    stat_2d.write(f"\tNumber of selected clusters for {core.sam2trim2d} is: {nb_selected}\n\n")
                    core.sys.stdout.write(
                        f"Sum of cluster sizes for {core.sam2trim2d} at locus {core.loc2trim2d} = {a}: with threshold set at {ts[0]}%,{core.warningStyle} clusters with size >= {b} were retained{core.successStyle}\n"
                        f"Number of selected clusters for sample {core.sam2trim2d}: {nb_selected}\n"
                        + core.normalStyle
                    )
                    if os.path.exists("trim-select." + core.scriptExt):
                        os.remove("trim-select." + core.scriptExt)
                    else:
                        core.sys.stdout.write(core.warningStyle + "\nNo data found to process for selected locus\n" + core.normalStyle)
                    if os.path.exists("tmp"):
                        os.remove("tmp")
                    if os.path.exists("tmp2"):
                        os.remove("tmp2")
                    core.os.chdir(core.current_dir)
                    stat_2d.flush()
