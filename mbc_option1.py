from mbc_option_common import prepare_with_previous_params, show_missing_loci_and_return


def main_menu1(core):
    """Displays submenu 1."""
    core.os.system("cls" if core.winOS else "clear")
    print(
        core.titleStyle
        + "\n--- MENU 1: BASIC ANALYSIS - only option 1 is strictly mandatory ---"
        + core.normalStyle
        + "\n\n"
        "1  -> INITIAL ANALYSIS ("
        + core.warningStyle
        + "mandatory"
        + core.normalStyle
        + "): read merging, sample-level dereplication, sequence clustering,\n\tchimera detection, affiliation of sequences to loci, and sequence re-orientation\n\n"
        "1a -> Re-analyze all loci, from the clustering step, modifying parameters\n\n"
        "1b -> Re-analyze only one locus of paired-end amplicon (merged reads), modifying parameters\n\n"
        "1c -> Re-analyze only one locus of single-end amplicon (R1 only), modifying parameters\n\n"
        "1d -> Re-analyse a given sample, modifying parameters\n\n"
        "1e -> "
        + core.warningStyle
        + "Optional"
        + core.normalStyle
        + " quality checking of fastq files (slow)\n"
        + core.normalStyle
    )
    core.rmenu = core.promptUser(
        "Please select an option among those listed above",
        None,
        ["1", "1a", "1b", "1c", "1d", "1e", "back", "home", "exit"],
        1,
        core.main,
        "",
    )

    if core.rmenu == "1e":
        menu1e(core)
    elif core.rmenu == "1":
        menu1(core)
    elif core.rmenu == "1a":
        menu1a(core)
    elif core.rmenu == "1b":
        menu1b(core)
    elif core.rmenu == "1c":
        menu1c(core)
    elif core.rmenu == "1d":
        menu1d(core)
    return core.rmenu


def menu1(core):
    """Runs option 1."""
    try:
        core.dir_fastq
    except AttributeError:
        core.in_dir_fastq()
    try:
        core.fastq_R1
    except AttributeError:
        core.in_fastq_R1()
    try:
        core.fastq_R2
    except AttributeError:
        core.in_fastq_R2()
    try:
        core.lociPE
    except AttributeError:
        core.in_lociPE()
    try:
        core.lociSE
    except AttributeError:
        core.in_lociSE()
    try:
        core.Samples
    except AttributeError:
        core.in_Samples()
    try:
        core.minsize
    except AttributeError:
        core.in_minsize()
    try:
        core.minseqlength
    except AttributeError:
        core.in_minseqlength()
    try:
        core.alpha
    except AttributeError:
        core.in_alpha()
    try:
        core.identity
    except AttributeError:
        core.in_identity()
    core.folders()
    core.param_1x()
    if len(core.lociPEs) > 0:
        core.merging()
    core.fastq2fas()
    core.derep_1()
    core.cluster_1x()
    core.chimera_remove()
    if len(core.lociPEs) > 0:
        core.runloc_merged()
    if len(core.lociSEs) > 0:
        core.runloc_r1()
    core.orient_1x()
    core.sys.stdout.write("\n\n")
    core.runs_1x()
    core.stats_1x()
    core.rerun(None)
    if len(core.sys.argv) < 2:
        core.rerun(None)


def menu1a(core):
    prepare_with_previous_params(core)
    core.in_minsize()
    core.in_minseqlength()
    core.in_alpha()
    core.in_identity()
    core.param_1x()
    print()
    core.cluster_1x()
    core.chimera_remove()
    core.runloc_merged()
    core.runloc_r1()
    core.orient_1x()
    core.runs_1x()
    core.stats_1x()
    core.rerun(None)


def menu1b(core):
    prepare_with_previous_params(core)
    if len(core.lociPEs) == 0:
        show_missing_loci_and_return(core, "No paired-end reads in current selection!", main_menu1)
    else:
        core.in_loc_sel_merged()
        core.in_minsize()
        core.in_minseqlength()
        core.in_alpha()
        core.in_identity()
        core.param_1x()
        core.cluster_1x()
        core.chimera_remove()
        core.runlocsel_merged()
        core.orient_1x()
        core.runs_1x()
        core.stats_1x()
        core.rerun(None)


def menu1c(core):
    prepare_with_previous_params(core)
    if len(core.lociSEs) == 0:
        show_missing_loci_and_return(core, "No R1/single-end reads in current selection!", main_menu1)
    else:
        core.in_loc_sel_r1()
        core.in_minsize()
        core.in_minseqlength()
        core.in_alpha()
        core.in_identity()
        core.param_1x()
        core.cluster_1x()
        core.runlocsel_r1()
        core.orient_1x()
        core.runs_1x()
        core.stats_1x()
        core.rerun(None)


def menu1d(core):
    prepare_with_previous_params(core)
    core.in_sam_sel()
    core.in_minsize()
    core.in_minseqlength()
    core.in_alpha()
    core.in_identity()
    core.param_1x()
    core.cluster_1x()
    core.chimera_remove()
    core.runloc_one_sample_1d()
    core.orient_1x()
    core.runs_1x()
    core.stats_1x()
    core.rerun(None)


def menu1e(core):
    prepare_with_previous_params(core)
    core.quality()
    core.runs_1x()
    core.rerun(None)


def _sync_core_globals(core):
    names = [
    "os", "sys", "subprocess",
    "start_log_redirect", "end_log_redirect", "main_stream_message", "logFileMessage",
    "scriptExt", "globalErrorOnStopCmd", "localErrorOnStopCmd", "fileSep",
    "samples", "fastq_R1", "fastq_R2", "fastqr1s", "fastqr2s", "lociPEs", "lociSEs", "dir_fastq",
    "rmenu", "sam_sel", "loc_sel1", "loc_sel2", "minsize", "minseqlength", "alpha", "identity",
    "shellCmd", "current_dir",
    "errorStyle", "warningStyle", "successStyle", "normalStyle",
    "winOS", "customExit", "promptUser", "main",
    ]
    module_globals = globals()
    for name in names:
        if hasattr(core, name):
            module_globals[name] = getattr(core, name)


def quality(core):
    """Tests the quality of each 'fastq' file (option 1e) by the VSEARCH command:

    vsearch --fastq_eestats2 dir_fastq/fastqF-R1 --output ../outputs/SN_R1_quality.txt
    vsearch --fastq_eestats2 dir_fastq/fastqF-R2 --output ../outputs/SN_R2_quality.txt

    dir_fastq = directory containing the fastq files
    fastqF-R1 and  fastqF-R2 = fastq file names for the sample
    SN = sample name
    """
    _sync_core_globals(core)
    global fastqr2s, fastqr1s
    with open("scripts/infor1." + scriptExt, "w") as out1:
        i = 0
        out1.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Quality statistical tests on R1 reads for samples:\n'))
        while i < len(samples):
            sample = samples[i]
            fastqr1 = fastqr1s[i]
            i = i + 1
            out1.write(main_stream_message(f' {sample}...') + logFileMessage(f'Quality statistical tests on R1 reads for sample {sample}'))
            out1.write(f"vsearch --fastq_qmax 93 --fastq_eestats2 \"{dir_fastq}{fileSep}{fastqr1}\" --output "
                       f"../outputs/{sample}_R1_quality.txt" + localErrorOnStopCmd + "\n")
        out1.write(main_stream_message(f'\n\n'))

    if len(lociPEs) > 0:
        with open("scripts/infor2." + scriptExt, "w") as out2:
            i = 0
            out2.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Quality statistical tests on R2 reads for samples:\n'))
            while i < len(samples):
                sample = samples[i]
                fastqr2 = fastqr2s[i]
                i = i + 1
                out2.write(main_stream_message(f' {sample}...') + logFileMessage(f'Quality statistical tests on R2 reads for sample {sample}'))
                out2.write(f"vsearch --fastq_qmax 93 --fastq_eestats2 \"{dir_fastq}{fileSep}{fastqr2}\" --output "
                           f"../outputs/{sample}_R2_quality.txt" + localErrorOnStopCmd + "\n")
            out2.write(main_stream_message(f'\n\n'))


def merging(core):
    """Merges paired-end reads into one sequence, when the length of the expected amplicon allows it (option 1)
    according the VSEARCH command:

    vsearch --fastq_mergepairs dir_fastq/fastqR1 --reverse dir_fastq/fastqR2 --fastaout
    ../tmp_files/SN_pairedEnd.fa --fastq_allowmergestagger --relabel sample=SN_merged.

    dir_fastq =  directory containing the fastq files
    fastqR1 = complete file name of fastq R1 read
    fastqR2 = fastqR2 complete name
    option --fastq_allowmergestagger allows the merging of short fragments
    """
    _sync_core_globals(core)
    with open("scripts/merging." + scriptExt, "w") as out:
        i = 0
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Merging paired-end reads for samples:\n'))
        while i < len(samples):
            sample = samples[i]
            fastqr1 = fastqr1s[i]
            fastqr2 = fastqr2s[i]
            out.write(main_stream_message(f" {sample}...") + logFileMessage(f'Merging paired-end reads for sample {sample}') +
                      f"vsearch --fastq_mergepairs \"{dir_fastq}{fileSep}{fastqr1}\" --reverse \"{dir_fastq}{fileSep}{fastqr2}\" "
                      f"--fastaout ../tmp_files/{sample}_pairedEnd.fa --fastq_allowmergestagger --relabel "
                      f"sample={sample}_merged." + localErrorOnStopCmd + "\n")
            i = i + 1
        out.write(main_stream_message(f'\n\n'))



def fastq2fas(core):
    """When the merging R1/R2 is impossible because of an unadapted size of amplicon, the reads R1 of 301 bp
    (better than R2) are used to search the relevant sequences.
    First, all R1 'fastq' files have to be transformed into 'fasta' files by the VSEARCH command:

    vsearch --fastq_filter dir_fastq/fastaqR1 –fastaout ../tmp_files/SN_singleEnd.fa

    Where : dir_fastq = directory containing the fastq files ; SN = sample name ; fastqR1 = name of the 'fastq' file
    containing the R1 reads
    """
    _sync_core_globals(core)
    with open("scripts/fqtofas." + scriptExt, "w") as out:
        i = 0
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Converting FASTQ files into FASTA format for samples:\n'))
        while i < len(samples):
            sample = samples[i]
            fastqr1 = fastqr1s[i]
            out.write(main_stream_message(f' {sample}...') + logFileMessage(f'Converting FASTQ files into FASTA format for sample {sample}') +
                      f"vsearch --fastq_qmax 93 --fastq_filter \"{dir_fastq}{fileSep}{fastqr1}\" --fastaout ../tmp_files/{sample}_singleEnd.fa --relabel sample={sample}_R1." + localErrorOnStopCmd + "\n")
            i = i + 1
        out.write(main_stream_message(f'\n\n'))


def derep_1(core):
    """Dereplicates merged sequences in a given 'fasta' file with the VSEARCH command:

    vsearch --fastx_uniques ../tmp_files/SN_pairedEnd.fa --fastaout ../tmp_files/SN_pairedEnd_derep.fas --sizeout
   --strand both

    And dereplicates the R1 sequences in a given 'fasta' file with the VSEARCH command:

    vsearch --fastx_uniques ../tmp_files/SN_singleEnd.fa --fastaout ../tmp_files/SN_singleEnd_derep.fas --sizeout
   --strand both

    Both commands dereplicate in both strands (option --strand both) and write abundance annotation (frequency)
    to output (option --sizeout).

    SN = sample name
    """
    _sync_core_globals(core)
    if len(lociPEs) > 0:
        with open("scripts/derep." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Dereplicating merged reads for samples:\n'))
            for sample in samples:
                out.write(main_stream_message(
                    f' {sample}...') + logFileMessage(f'Dereplicating merged reads for sample {sample}') +
                          f"vsearch --fastx_uniques ../tmp_files/{sample}_pairedEnd.fa --fastaout "
                          f"../tmp_files/{sample}_pairedEnd_derep.fas --sizeout --strand both" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

    if len(lociSEs) > 0:
        with open("scripts/derep_r1." + scriptExt, "w") as out1:
            out1.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Dereplicating R1/single-end reads for samples:\n'))
            for sample in samples:
                out1.write(main_stream_message(f' {sample}...') + logFileMessage(f'Dereplicating R1/single-end reads for sample {sample}') +
                           f"vsearch --fastx_uniques ../tmp_files/{sample}_singleEnd.fa --fastaout "
                           f"../tmp_files/{sample}_singleEnd_derep.fas --sizeout --strand both" + localErrorOnStopCmd
                           + "\n")
            out1.write(main_stream_message(f'\n\n'))


def cluster_1x(core):
    """Denoises and clusters Illumina dereplicated merged sequences and gives in output the centroids sequences
    to 'fasta' files (options 1, 1a, 1b, 1c and 1d) with the following VSEARCH commands:

    For loci based in paired-end reads:

    vsearch --cluster_unoise ../tmp_files/SN_pairedEnd_derep.fas --sizein --centroids
    ../tmp_files/SN(or SS)_pairedEnd_cluster.fas --strand both --minsize int --sizeout --unoise_alph int
   --minseqlength int

    For loci based in R1/single-end reads:

    vsearch --cluster_unoise ../tmp_files/SN_singleEnd_derep.fas --sizein --centroids
    ../tmp_files/SN(or SS)_singleEnd_cluster.fas --strand both --minsize int --sizeout --unoise_alph int
   --minseqlength int

    SN = sample name; int = integer; SS = selected sample
    """
    _sync_core_globals(core)
    global rmenu
    if len(lociPEs) > 0 and  rmenu in ["1", "1a", "1b"] and len(lociPEs) > 0:
        with open("scripts/cluster." + scriptExt, "w") as out:
            i = 0
            if rmenu != "1":
                print()
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Clustering merged reads for all samples:\n'))
            while i < len(samples):
                sample = samples[i]
                out.write(main_stream_message(f' {sample}...') + logFileMessage(f'Clustering merged reads for sample {sample}') +
                          f"vsearch --cluster_unoise ../tmp_files/{sample}_pairedEnd_derep.fas --sizein --centroids "
                          f"../tmp_files/{sample}_pairedEnd_cluster.fas --strand both --minsize {minsize} --sizeout "
                          f"--unoise_alph {alpha} --minseqlength {minseqlength}" + localErrorOnStopCmd + "\n")
                i = i + 1
            out.write(main_stream_message(f'\n\n'))

    if len(lociSEs) > 0 and rmenu in ["1", "1a", "1c"]:
        with open("scripts/cluster_r1." + scriptExt, "w") as out:
            i = 0
            if rmenu == "1c":
                print()
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Clustering R1/single-end reads for all samples:\n'))
            while i < len(samples):
                sample = samples[i]
                out.write(main_stream_message(f' {sample}...') + logFileMessage(f'Clustering R1/single-end reads for sample {sample}') +
                          f"vsearch --cluster_unoise ../tmp_files/{sample}_singleEnd_derep.fas --sizein --centroids "
                          f"../tmp_files/{sample}_singleEnd_cluster.fas --strand both --minsize {minsize} --sizeout "
                          f"--unoise_alpha {alpha} --minseqlength {minseqlength}" + localErrorOnStopCmd
                          + "\n")
                i = i + 1
            out.write(main_stream_message(f'\n\n'))

    if rmenu == "1d":
        with open("scripts/cluster_one_sample_1d." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'\nClustering reads for selected sample {sam_sel}...'))
            if len(lociPEs) > 0:
                out.write(logFileMessage(f'Clustering merged reads for selected sample {sam_sel}') + f"vsearch --cluster_unoise ../tmp_files/{sam_sel}_pairedEnd_derep.fas --sizein --centroids "
                  f"../tmp_files/{sam_sel}_pairedEnd_cluster.fas --strand both --minsize {minsize} --sizeout "
                  f"--unoise_alph {alpha} --minseqlength {minseqlength}" + localErrorOnStopCmd + "\n")

            if len(lociSEs) > 0:
                out.write(logFileMessage(f'Clustering R1/single-end reads for selected sample {sam_sel}') + f"vsearch --cluster_unoise ../tmp_files/{sam_sel}_singleEnd_derep.fas --sizein --centroids "
                  f"../tmp_files/{sam_sel}_singleEnd_cluster.fas --strand both --minsize {minsize} "
                  f"--sizeout --unoise_alph {alpha} --minseqlength {minseqlength}" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))


def chimera_remove(core):
    """Detects and removes potential chimeras in denoised merged sequences or single-end or selected sample
    (options 1, 1a, 1b, 1c and 1d), by the VSEARCH commands:

    vsearch --uchime3_denovo ../tmp_files/SN(or SS)_pairedEnd_cluster.fas --nonchimeras
    ../tmp_files/SN(or SS)_pairedEnd_cluster_OK.fas
    vsearch --uchime3_denovo ../tmp_files/SN(or SS)_singleEnd_cluster.fas --nonchimeras
    ../tmp_files/SN(or SS)_singleEnd_cluster_OK.fas

    SN = sample name; SS = selected sample
    """
    _sync_core_globals(core)

    if len(lociPEs) > 0 and rmenu in ["1", "1a", "1b"]:
        with open("scripts/chimera." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd +
                      "\n" + main_stream_message(f'Detecting and removing chimeras within merged reads of samples:\n'))
            for sample in samples:
                out.write(main_stream_message(f' {sample}...') + logFileMessage(f'Detecting and removing chimeras within merged reads of sample {sample}') +
                          f"vsearch --uchime3_denovo ../tmp_files/{sample}_pairedEnd_cluster.fas --nonchimeras "
                          f"../tmp_files/{sample}_pairedEnd_cluster_OK.fas" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

    if len(lociSEs) > 0 and rmenu in ["1", "1a", "1c"]:
        with open("scripts/chimera_r1." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Detecting and removing chimeras within R1/single-end reads of samples:\n'))
            for sample in samples:
                out.write(main_stream_message(f' {sample}...') + logFileMessage(f'Detecting and removing chimeras within R1/single-end reads of sample {sample}') +
                          f"vsearch --uchime3_denovo ../tmp_files/{sample}_singleEnd_cluster.fas --nonchimeras"
                          f" ../tmp_files/{sample}_singleEnd_cluster_OK.fas" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

    if rmenu == "1d":
        with open("scripts/chimera_one_sample_1d." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Detecting and removing chimeras for selected sample {sam_sel}...'))
            if len(lociPEs) > 0:
                out.write(logFileMessage(f'Detecting and removing chimeras within merged reads of selected sample {sam_sel}') +
                    f"vsearch --uchime3_denovo ../tmp_files/{sam_sel}_pairedEnd_cluster.fas --nonchimeras "
                            f"../tmp_files/{sam_sel}_pairedEnd_cluster_OK.fas" + localErrorOnStopCmd + "\n")
            if len(lociSEs) > 0:
                out.write(logFileMessage(f'Detecting and removing chimeras within R1/single-end reads of selected sample {sam_sel}') +
                    f"vsearch --uchime3_denovo ../tmp_files/{sam_sel}_singleEnd_cluster.fas --nonchimeras "
                      f"../tmp_files/{sam_sel}_singleEnd_cluster_OK.fas" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))


def runloc_merged(core):
    """Searches similarities between merged sequences, denoised and non-chimera sequences and the local reference
    database (-db), options 1 and 1a, by the VSEARCH command:

    vsearch --usearch_global ../tmp_files/SN_pairedEnd_cluster_OK.fas --db ../refs/L1.fas --matched
    ../results_by_locus/L1/SN_pairedEnd.fas --id int --strand both

    L1 = locus name for amplicons based on paired-end reads
    SN = sample name
    id = minimum identity accepted (0-1.0)
    """
    _sync_core_globals(core)
    with open("scripts/results_by_locusmerged." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters to loci for merged reads of samples:\n'))
        for lociPEb in lociPEs:
            for sample in samples:
                out.write(main_stream_message(
                    f' {sample} vs {lociPEb}...') + logFileMessage(f'Affiliating clusters to locus {lociPEb} for merged reads of sample {sample}') +
                          f"vsearch --usearch_global ../tmp_files/{sample}_pairedEnd_cluster_OK.fas --db ../refs/{lociPEb}.fas"
                          f" --matched ../results_by_locus/{lociPEb}/{sample}_pairedEnd.fas --id {identity} --strand both"
                          + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))


def runloc_r1(core):
    """Searches similarities between single-end sequences (in case of amplicons where the merging R1/R2 is impossible)
    denoised and non-chimera sequences and the local reference database (-db), options 1 and 1a,
    by the VSEARCH command:

     vsearch --usearch_global ../tmp_files/SN_singleEnd_cluster_OK.fas --db ../refs/L2.fas --matched
     ../results_by_locus/L2/SN_singledEnd.fa --id int --strand both

    L2 = locus for amplicons with no mergeable R1/R2
    SN = sample name
    id = minimum identity accepted (0-1.0)
    """
    _sync_core_globals(core)
    with open("scripts/results_by_locusr1." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters to loci for R1/single-end reads of samples:\n'))
        for locusSEb in lociSEs:
            for sample in samples:
                out.write(main_stream_message(f' {sample} vs {locusSEb}...') + logFileMessage(f'Affiliating clusters to locus {locusSEb} for R1/single-end reads of sample {sample}') +
                          f"vsearch --usearch_global ../tmp_files/{sample}_singleEnd_cluster_OK.fas --db "
                          f"../refs/{locusSEb}.fas --matched ../results_by_locus/{locusSEb}/{sample}_singleEnd.fas --id {identity} "
                          f"--strand both" + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))


def runlocsel_merged(core):
    """Searches similarities between merged sequences, denoised and non-chimera sequences and the local reference
    database (-db), for a selected paired-end based locus, option 1b, by the VSEARCH command:

    vsearch --usearch_global ../tmp_files/SN_pairedEnd_cluster_OK.fas --db ../refs/SL.fas --matched
    ../results_by_locus/SL/SN_merged.fas --id real --strand both

    SN = sample name
    SL = Selected locus based on paired-end reads
    id = minimum identity accepted (0-1.0)
    """
    _sync_core_globals(core)
    with open("scripts/results_by_locus_sel." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters to selected locus {loc_sel1} for merged reads of samples:\n'))
        for sample in samples:
            out.write(main_stream_message(f' {sample}...') + logFileMessage(f'Affiliating clusters to selected locus {loc_sel1} for merged reads of sample {sample}') +
                      f"vsearch --usearch_global ../tmp_files/{sample}_pairedEnd_cluster_OK.fas --db "
                      f"../refs/{loc_sel1}.fas --matched ../results_by_locus/{loc_sel1}/{sample}_pairedEnd.fas --id {identity} "
                      f"--strand both"
                      + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))


def runlocsel_r1(core):
    """Searches similarities between single-end sequences (R1), denoised and non-chimera sequences and the local
    reference database (-db), for a selected single-end based locus, option 1c, by the VSEARCH command:

    vsearch --usearch_global ../tmp_files/SN_singleEnd_cluster_OK.fas --db ../refs/SL.fas --matched
    ../results_by_locus/SL/SN_singleEnd.fas --id real --strand both

    SN = sample name
    SL = Selected locus based on single-end read (R1)
    id = minimum identity accepted (0-1.0)
    """
    _sync_core_globals(core)
    with open("scripts/results_by_locusr1_sel." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters to selected locus {loc_sel2} for R1/single-end reads of samples:\n'))
        for sample in samples:
            out.write(main_stream_message(f' {sample}...') + logFileMessage(f'Affiliating clusters to selected locus {loc_sel2} for R1/single-end reads of sample {sample}') +
                      f"vsearch --usearch_global ../tmp_files/{sample}_singleEnd_cluster_OK.fas --db "
                      f"../refs/{loc_sel2}.fas --matched ../results_by_locus/{loc_sel2}/{sample}_singleEnd.fas --id {identity} "
                      f"--strand both"
                      + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))


def runloc_one_sample_1d(core):
    """ Searches similarities between denoised and non-chimera sequences and local reference database (db)
    only for a selected sample (option 1d) by the VSEARCH command:

    vsearch --usearch_global ../tmp_files/SS_pairedEnd_cluster_OK.fas --db ../refs/L1.fas --matched
    ../results_by_locus/L1/SS_pairedEnd.fas --id real --strand both

    vsearch --usearch_global ../tmp_files/SS_singleEnd_cluster_OK.fas --db ../refs/L2.fas --matched
    ../results_by_locus/L2/SS_singleEnd.fas --id real --strand both

    SS = selected sample
    L1 = locus name for amplicons based on paired-end reads
    L2 = locus name for amplicons with no mergeable R1/R2
    id = minimum identity accepted (0-1.0)
    """
    _sync_core_globals(core)
    with open("scripts/results_by_locus_merged_1d." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters to loci for merged reads of selected sample {sam_sel}:\n'))
        for lociPEb in lociPEs:
            out.write(main_stream_message(f' {lociPEb}...') + logFileMessage(f'Affiliating clusters to locus {lociPEb} for merged reads of selected sample {sam_sel}') +
                      f"vsearch --usearch_global ../tmp_files/{sam_sel}_pairedEnd_cluster_OK.fas --db "
                      f"../refs/{lociPEb}.fas --matched ../results_by_locus/{lociPEb}/{sam_sel}_pairedEnd.fas --id {identity} "
                      f"--strand both" + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))

    with open("scripts/results_by_locus_R1_1d." + scriptExt, "w") as out1:
        out1.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters to loci for R1/single-end reads of selected sample {sam_sel}:\n'))
        for locusSEb in lociSEs:
            out1.write(main_stream_message(f' {locusSEb}...') + logFileMessage(f'Affiliating clusters to locus {locusSEb} for R1/single-end reads of selected sample {sam_sel}') +
                       f"vsearch --usearch_global ../tmp_files/{sam_sel}_singleEnd_cluster_OK.fas --db "
                       f"../refs/{locusSEb}.fas --matched ../results_by_locus/{locusSEb}/{sam_sel}_singleEnd.fas "
                       f"--id {identity} --strand both"
                       + localErrorOnStopCmd + "\n")
            out1.write(main_stream_message(f'\n\n'))


def orient_1x(core):
    """Orients all the sequences in the same direction (forward) than references, options 1, 1a, 1c and 1d,
    with the following script:

    For paired-end merged sequences:
    vsearch --orient ../results_by_locus/L1/SN_pairedEnd.fas --db ../refs/L1.fas --fastaout ../results_by_locus/L1/SN_pairedEnd_orient.fas

    For single-end R1 sequences:
    vsearch --orient ../results_by_locus/L2/SN_singleEnd.fas --db ../refs/L2.fas --fastaout ../results_by_locus/L2/SN_singleEnd_orient.fas

    For selected loci:
    vsearch --orient ../results_by_locus/SL/SN_pairedEnd.fas --db ../refs/SL.fas --fastaout ../results_by_locus/SL/SN_pairedEnd_orient.fas
    or
    vsearch --orient ../results_by_locus/SL/SN_singleEnd.fas --db ../refs/SL.fas --fastaout ../results_by_locus/SL/SN_singleEnd_orient.fas

    For selected sample:
    vsearch --orient ../results_by_locus/L1/SS_pairedEnd.fas --db ../refs/L1.fas --fastaout ../results_by_locus/L1/SS_pairedEnd_orient.fas
    and
    vsearch --orient ../results_by_locus/L2/SS_singleEnd.fas --db ../refs/L2.fas --fastaout ../results_by_locus/L1/SS_singleEnd_orient.fas

    SN = locus name; SS = selected sample; SL = selected locus
    L1 = locus name for amplicons based on paired-end reads
    L2 = locus name for amplicons with no mergeable R1/R2
    """
    _sync_core_globals(core)

    if rmenu in ["1", "1a"]:
        if len(lociPEs) > 0:
            with open("scripts/orient_merged_1a." + scriptExt, "w") as out:
                out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f"Orienting all merged reads according to the reference:\n"))
                for locusPEb in lociPEs:
                    for sample in samples:
                        out.write(main_stream_message(f' {sample} vs {locusPEb}...') + logFileMessage(f'Orienting all merged reads according to the reference for locus {locusPEb} and sample {sample}') +
                                  f"vsearch --orient ../results_by_locus/{locusPEb}/{sample}_pairedEnd.fas --db ../refs/{locusPEb}.fas "
                                  f"--fastaout ../results_by_locus/{locusPEb}/{sample}_pairedEnd_orient.fas" + localErrorOnStopCmd + "\n")
                out.write(main_stream_message(f'\n\n'))

        if len(lociSEs) > 0:
            with open("scripts/orient_R1_1a." + scriptExt, "w") as out1:
                out1.write(globalErrorOnStopCmd + "\n" + main_stream_message(f"Orienting all R1/single-end reads according to the reference:\n"))
                for locusSEb in lociSEs:
                    for sample in samples:
                        out1.write(main_stream_message(f' {sample} vs {locusSEb}...') + logFileMessage(f'Orienting all R1/single-end reads according to the reference for locus {locusSEb} and sample {sample}') +
                                   f"vsearch --orient ../results_by_locus/{locusSEb}/{sample}_singleEnd.fas --db ../refs/{locusSEb}.fas "
                                   f"--fastaout ../results_by_locus/{locusSEb}/{sample}_singleEnd_orient.fas --tabbedout ../results_by_locus/{locusSEb}/{sample}_singleEnd_orient.tsv" + localErrorOnStopCmd + "\n")
                out1.write(main_stream_message(f'\n\n'))

    elif rmenu == "1b":
        with open("scripts/orient_merged_1b." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Orienting all merged reads according to the reference for selected locus {loc_sel1} and all samples:\n'))
            for sample in samples:
                out.write(main_stream_message(f' {sample}...') + logFileMessage(f'Orienting all merged reads according to the reference for selected locus {loc_sel1} and {sample}') +
                          f"vsearch --orient ../results_by_locus/{loc_sel1}/{sample}_pairedEnd.fas --db "
                          f"../refs/{loc_sel1}.fas --fastaout ../results_by_locus/{loc_sel1}/{sample}_pairedEnd_orient.fas" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

    elif rmenu == "1c":
        with open("scripts/orient_R1_1c." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Orienting all R1/single-end reads according to the reference for selected locus {loc_sel2} and all samples:\n'))
            for sample in samples:
                out.write(main_stream_message(f' {sample}...') + logFileMessage(f'Orienting all R1/single-end reads according to the reference for selected locus {loc_sel2} and {sample}') +
                          f"vsearch --orient ../results_by_locus/{loc_sel2}/{sample}_singleEnd.fas --db ../refs/{loc_sel2}.fas "
                          f"--fastaout ../results_by_locus/{loc_sel2}/{sample}_singleEnd_orient.fas --tabbedout ../results_by_locus/{loc_sel2}/{sample}_singleEnd_orient.tsv" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

    if rmenu == "1d":
        if len(lociPEs) > 0:
            with open("scripts/orient_merged_1d." + scriptExt, "w") as out:
                out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Orienting all merged reads of selected sample {sam_sel} for all R1/R2 loci:\n'))
                for locusPEb in lociPEs:
                    out.write(main_stream_message(f' {locusPEb}...') + logFileMessage(f'Orienting all merged reads of selected sample {sam_sel} for locus {locusPEb}') +
                              f"vsearch --orient ../results_by_locus/{locusPEb}/{sam_sel}_pairedEnd.fas --db ../refs/{locusPEb}.fas --fastaout "
                              f"../results_by_locus/{locusPEb}/{sam_sel}_pairedEnd_orient.fas" + localErrorOnStopCmd + "\n")
                out.write(main_stream_message(f'\n\n'))

        if len(lociSEs) > 0:
            with open("scripts/orient_R1_1d." + scriptExt, "w") as out18:
                out18.write(main_stream_message(f'Orienting all R1/single-end reads of selected sample {sam_sel} for all R1/single-end loci:\n'))
                for locusSEb in lociSEs:
                    out18.write(main_stream_message(f' {locusSEb}...') + logFileMessage(f'Orienting all R1/single-end reads of selected sample {sam_sel} for locus {locusSEb}') +
                                f"vsearch --orient ../results_by_locus/{locusSEb}/{sam_sel}_singleEnd.fas --db ../refs/{locusSEb}.fas "
                                f"--fastaout ../results_by_locus/{locusSEb}/{sam_sel}_singleEnd_orient.fas --tabbedout ../results_by_locus/{locusSEb}/{sam_sel}_singleEnd_orient.tsv" + localErrorOnStopCmd + "\n")
                out18.write(main_stream_message(f'\n\n'))


def runs_1x(core):
    """Creation of global scripts for the options 1, 1a, 1b, 1c, Ad and 1e
    """
    _sync_core_globals(core)
    global rmenu, p
    if rmenu == "1":
        with open(f"scripts/runall1.{scriptExt}", "w") as out:
            if not winOS:
                out.write(start_log_redirect('../outputs/res1.log') +
                          "scriptArray=(" + (("'./merging." + scriptExt + "' ") if len(lociPEs) > 0 else "")
                            + "'./fqtofas." + scriptExt + "' "
                            + (("'./derep." + scriptExt + "' ") if len(lociPEs) > 0 else "")
                            + (("'./derep_r1." + scriptExt + "' ") if len(lociSEs) > 0 else "")
                            + (("'./cluster." + scriptExt + "' ") if len(lociPEs) > 0 else "")
                            + (("'./cluster_r1." + scriptExt + "' ") if len(lociSEs) > 0 else "")
                            + (("'./chimera." + scriptExt + "' ") if len(lociPEs) > 0 else "")
                            + (("'./chimera_r1." + scriptExt + "' ") if len(lociSEs) > 0 else "")
                            + (("'./results_by_locusmerged." + scriptExt + "' ") if len(lociPEs) > 0 else "")
                            + (("'./results_by_locusr1." + scriptExt + "' ") if len(lociSEs) > 0 else "")
                            + (("'./orient_merged_1a." + scriptExt + "' ") if len(lociPEs) > 0 else "")
                            + (("'./orient_R1_1a." + scriptExt + "' ") if len(lociSEs) > 0 else "")
                            + ')\nfor script in "${scriptArray[@]}"\n' +
                              '     do\n' +
                              '     if ! ${script}; then\n' +
                              '             printf "Error executing ${script}\\n\\n" >&3\n' +
                              '             exit 1\n' +
                              '     fi\n' +
                              'done\n' +
                          end_log_redirect('../outputs/res1.log'))
            else:
                out.write(start_log_redirect('../outputs/res1.log') +
                            '   $scriptArray = @(' + (("\"./merging." + scriptExt + "\", ") if len(lociPEs) > 0 else "")
                            + '"./fqtofas.' + scriptExt + '", '
                            + (("\"./derep." + scriptExt + "\", ") if len(lociPEs) > 0 else "")
                            + (("\"./derep_r1." + scriptExt + "\", ") if len(lociSEs) > 0 else "")
                            + (("\"./cluster." + scriptExt + "\", ") if len(lociPEs) > 0 else "")
                            + (("\"./cluster_r1." + scriptExt + "\", ") if len(lociSEs) > 0 else "")
                            + (("\"./chimera." + scriptExt + "\", ") if len(lociPEs) > 0 else "")
                            + (("\"./chimera_r1." + scriptExt + "\", ") if len(lociSEs) > 0 else "")
                            + (("\"./results_by_locusmerged." + scriptExt + "\", ") if len(lociPEs) > 0 else "")
                            + (("\"./results_by_locusr1." + scriptExt + "\", ") if len(lociSEs) > 0 else "")
                            + (("\"./orient_merged_1a." + scriptExt + '"' + (', ' if len(lociSEs) > 0 else "")) if len(lociPEs) > 0 else "")
                            + (("\"./orient_R1_1a." + scriptExt + '"') if len(lociSEs) > 0 else "")
                            + ')\n   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
                            '               $script = $scriptArray[$i]\n' +
                            '               & "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; '
                            'exit $LASTEXITCODE }\n' +
                            '   }\n' +
                          end_log_redirect('../outputs/res1.log'))
        os.chdir('scripts')
        if not winOS:
            for file in os.listdir("."):
                os.chmod(file, 0o755)
        p = subprocess.run([shellCmd, "./runall1." + scriptExt])
        if p.returncode > 0:
            print(errorStyle + f"\nMain analysis execution failed with error code {p.returncode}" + normalStyle + f", please check {current_dir}" + fileSep + "outputs" + fileSep + "res1.log")
            customExit(1)
        os.chdir('..')
        print(successStyle + f"Step {rmenu} ended successfully" + normalStyle)

    if rmenu == "1a":
        with open(f"scripts/runall1a.{scriptExt}", "w") as out:
            if not winOS:
                out.write(start_log_redirect('../outputs/res1a.log') + "scriptArray=("
                                                    + (("'./cluster." + scriptExt + "' ") if len(lociPEs) > 0 else "")
                                                    + (("'./cluster_r1." + scriptExt + "' ") if len(lociSEs) > 0 else "")
                                                    + (("'./chimera." + scriptExt + "' ") if len(lociPEs) > 0 else "")
                                                    + (("'./chimera_r1." + scriptExt + "' ") if len(lociSEs) > 0 else "")
                                                    + (("'./results_by_locusmerged." + scriptExt + "' ") if len(lociPEs) > 0 else "")
                                                    + (("'./results_by_locusr1." + scriptExt + "' ") if len(lociSEs) > 0 else "")
                                                    + (("'./orient_merged_1a." + scriptExt + "' ") if len(lociPEs) > 0 else "")
                                                    + (("'./orient_R1_1a." + scriptExt + "' ") if len(lociSEs) > 0 else "")
                                                    + ')\nfor script in "${scriptArray[@]}"\n' +
                                                      '     do\n' +
                                                      '     if ! ${script}; then\n' +
                                                      '             printf "Error executing ${script}\\n\\n" >&3\n' +
                                                      '             exit 1\n' +
                                                      '     fi\n' +
                                                      'done\n' +
                          end_log_redirect('../outputs/res1a.log'))
            else:
                out.write(start_log_redirect('../outputs/res1a.log') +
                          '   $scriptArray = @('
                            + (("\"./cluster." + scriptExt + "\", ") if len(lociPEs) > 0 else "")
                            + (("\"./cluster_r1." + scriptExt + "\", ") if len(lociSEs) > 0 else "")
                            + (("\"./chimera." + scriptExt + "\", ") if len(lociPEs) > 0 else "")
                            + (("\"./chimera_r1." + scriptExt + "\", ") if len(lociSEs) > 0 else "")
                            + (("\"./results_by_locusmerged." + scriptExt + "\", ") if len(lociPEs) > 0 else "")
                            + (("\"./results_by_locusr1." + scriptExt + "\", ") if len(lociSEs) > 0 else "")
                            + (("\"./orient_merged_1a." + scriptExt + '"' + (', ' if len(lociSEs) > 0 else "")) if len(lociPEs) > 0 else "")
                            + (("\"./orient_R1_1a." + scriptExt + '"') if len(lociSEs) > 0 else "")
                            + ')\n   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
                              '             $script = $scriptArray[$i]\n' +
                              '             & "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; '
                              'exit $LASTEXITCODE }\n' +
                              '   }\n' +
                          end_log_redirect('../outputs/res1a.log'))
        os.chdir('scripts')
        if not winOS:
            for file in os.listdir("."):
                os.chmod(file, 0o755)
        p = subprocess.run([shellCmd, "./runall1a." + scriptExt])
        if p.returncode > 0:
            print(errorStyle + f"\nMain analysis execution failed with error code {p.returncode}" + normalStyle + f", please check {current_dir}" + fileSep + "outputs" + fileSep + "res1a.log")
            customExit(1)
        os.chdir('..')
        print(successStyle + f"Step {rmenu} ended successfully" + normalStyle)

    if rmenu == "1b":
        with open(f"scripts/runall1b.{scriptExt}", "w") as out:
            if not winOS:
                out.write(start_log_redirect('../outputs/res1b.log') + "scriptArray=('"
                                                   "./cluster." + scriptExt + "' '"
                                                   "./chimera." + scriptExt + "' '"
                                                   "./results_by_locus_sel." + scriptExt + "' '"
                                                   "./orient_merged_1b." + scriptExt + "')\n" +
                                                      'for script in "${scriptArray[@]}"\n' +
                                                      '     do\n' +
                                                      '     if ! ${script}; then\n' +
                                                      '             printf "Error executing ${script}\\n\\n" >&3\n' +
                                                      '             exit 1\n' +
                                                      '     fi\n' +
                                                      'done\n' +
                          end_log_redirect('../outputs/res1b.log'))
            else:
                out.write(start_log_redirect('../outputs/res1b.log') +
                          '   $scriptArray = @("'
                          './cluster.' + scriptExt + '", "'
                                                     './chimera.' + scriptExt + '", "'
                                                     './results_by_locus_sel.' + scriptExt + '", "'
                                                     './orient_merged_1b.' + scriptExt + '")\n' +
                                                      '   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
                                                      '             $script = $scriptArray[$i]\n' +
                                                      '             & "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; '
                                                      'exit $LASTEXITCODE }\n' +
                                                      '   }\n' +
                          end_log_redirect('../outputs/res1b.log'))
        os.chdir('scripts')
        if not winOS:
            for file in os.listdir("."):
                os.chmod(file, 0o755)
        p = subprocess.run([shellCmd, "./runall1b." + scriptExt])
        if p.returncode > 0:
            print(errorStyle + f"\nMain analysis execution failed with error code {p.returncode}" + normalStyle + f", please check {current_dir}" + fileSep + "outputs" + fileSep + "res1b.log")
            customExit(1)
        print(successStyle + f"Step {rmenu} ended successfully" + normalStyle)

    if rmenu == "1c":
        with open(f"scripts/runall1c.{scriptExt}", "w") as out:
            if not winOS:
                out.write(start_log_redirect('../outputs/res1c.log') + "scriptArray=('"
                                                                       "./cluster_r1." + scriptExt + "' '"
                                                                       "./chimera_r1." + scriptExt + "' '"
                                                                       "./results_by_locusr1_sel." + scriptExt + "' '"
                                                                       "./orient_R1_1c." + scriptExt + "')\n" +
                                                                          'for script in "${scriptArray[@]}"\n' +
                                                                          '     do\n' +
                                                                          '     if ! ${script}; then\n' +
                                                                          '             printf "Error executing ${script}\\n\\n" >&3\n' +
                                                                          '             exit 1\n' +
                                                                          '     fi\n' +
                                                                          'done\n' +
                          end_log_redirect('../outputs/res1c.log'))
            else:
                out.write(start_log_redirect('../outputs/res1c.log') +
                          '   $scriptArray = @("'
                          './cluster_r1.' + scriptExt + '", "'
                                                        './chimera_r1.' + scriptExt + '", "'
                                                        './results_by_locusr1_sel.' + scriptExt + '", "'
                                                        './orient_R1_1c.' + scriptExt + '")\n' +
                                                          '   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
                                                          '             $script = $scriptArray[$i]\n' +
                                                          '             & "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; '
                                                          'exit $LASTEXITCODE }\n' +
                                                          '   }\n' +
                          end_log_redirect('../outputs/res1c.log'))
        os.chdir('scripts')
        if not winOS:
            for file in os.listdir("."):
                os.chmod(file, 0o755)
        p = subprocess.run([shellCmd, "./runall1c." + scriptExt])
        if p.returncode > 0:
            print(errorStyle + f"\nMain analysis execution failed with error code {p.returncode}" + normalStyle + f", please check {current_dir}" + fileSep + "outputs" + fileSep + "res1c.log")
            customExit(1)
        os.chdir('..')
        print(successStyle + f"Step {rmenu} ended successfully" + normalStyle)

    if rmenu == "1d":
        with open(f"scripts/runall1d.{scriptExt}", "w") as out:
            if not winOS:
                out.write(start_log_redirect('../outputs/res1d.log') +
                          "scriptArray=('./cluster_one_sample_1d." + scriptExt + "' './chimera_one_sample_1d." + scriptExt + "' "
                                                                + (("'./results_by_locus_merged_1d." + scriptExt + "' ") if len(lociPEs) > 0 else "")
                                                                + (("'./results_by_locus_R1_1d." + scriptExt + "' ") if len(lociSEs) > 0 else "")
                                                                + (("'./orient_merged_1d." + scriptExt + "' ") if len(lociPEs) > 0 else "")
                                                                + (("'./orient_R1_1d." + scriptExt + "' ") if len(lociSEs) > 0 else "")
                                                                +        ')\nfor script in "${scriptArray[@]}"\n' +
                                                                      '     do\n' +
                                                                      '     if ! ${script}; then\n' +
                                                                      '             printf "Error executing ${script}\\n\\n" >&3\n' +
                                                                      '             exit 1\n' +
                                                                      '     fi\n' +
                                                                      'done\n' +
                          end_log_redirect('../outputs/res1d.log'))
            else:
                out.write(start_log_redirect('../outputs/res1d.log') +
                          '   $scriptArray = @('
                            + '"./cluster_one_sample_1d.' + scriptExt + '", "./chimera_one_sample_1d.' + scriptExt + '", '
                            + (("\"./results_by_locus_merged_1d." + scriptExt + "\", ") if len(lociPEs) > 0 else "")
                            + (("\"./results_by_locus_R1_1d." + scriptExt + "\", ") if len(lociSEs) > 0 else "")
                            + (("\"./orient_merged_1d." + scriptExt + '"' + (', ' if len(lociSEs) > 0 else "")) if len(lociPEs) > 0 else "")
                            + (("\"./orient_R1_1d." + scriptExt + '"') if len(lociSEs) > 0 else "")
                            + ')\n   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
                              '             $script = $scriptArray[$i]\n' +
                              '             & "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; '
                              'exit $LASTEXITCODE }\n' +
                          '   }\n' +
                          end_log_redirect('../outputs/res1d.log'))
        os.chdir('scripts')
        if not winOS:
            for file in os.listdir("."):
                os.chmod(file, 0o755)
        p = subprocess.run([shellCmd, "./runall1d." + scriptExt])
        if p.returncode > 0:
            print(errorStyle + f"\nMain analysis execution failed with error code {p.returncode}" + normalStyle + f", please check {current_dir}" + fileSep + "outputs" + fileSep + "res1d.log")
            customExit(1)
        os.chdir('..')
        print(successStyle + f"Step {rmenu} ended successfully" + normalStyle)

    if rmenu == "1e":
        with open("scripts/runall1e." + scriptExt, "w") as out:
            if not winOS:
                out.write(start_log_redirect('../outputs/res1e.log') +
                          "scriptArray=('"
                          "./infor1." + scriptExt + "' '"
                                                    "./infor2." + scriptExt + "')\n" +
                          'for script in "${scriptArray[@]}"\n' +
                          '     do\n' +
                          '     if ! ${script}; then\n' +
                          '             printf "Error executing ${script}\\n\\n" >&3\n' +
                          '             exit 1\n' +
                          '     fi\n' +
                          'done\n' + end_log_redirect('../outputs/res1e.log'))
            else:
                out.write(start_log_redirect('../outputs/res1e.log') +
                          '   $scriptArray = @("'
                          './infor1.' + scriptExt + '", "'
                                                    './infor2.' + scriptExt + '")\n' +
                          '   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
                          '             $script = $scriptArray[$i]\n' +
                          '             & "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; exit $LASTEXITCODE }\n' +
                          '   }\n' + end_log_redirect('../outputs/res1e.log'))
        os.chdir('scripts')
        if not winOS:
            for file in os.listdir("."):
                os.chmod(file, 0o755)
        sys.stdout.write("Quality checking is being processed " + warningStyle + "(slow procedure, be patient!)\n\n" + normalStyle)
        p = subprocess.run([shellCmd, "./runall1e." + scriptExt])
        if p.returncode > 0:
            print(str(p))
            print(errorStyle + errorStyle + f"\nMain analysis execution failed with error code {p.returncode}" + normalStyle + f", please check {current_dir}" + fileSep + "outputs" + fileSep + "res1e.log" + normalStyle)
            customExit(1)
        print(successStyle + f"Step {rmenu} ended successfully" + normalStyle)
        sys.stdout.write(successStyle + "\nExecution of option 1e is complete\n" + normalStyle + f"Statistical test files (*_quality.txt) are located at ---> {current_dir}{fileSep}outputs\n\n")
        os.chdir(current_dir)


def stats_1x(core):
    """Calculates the number of resulting sequences (reads, merged, dereplicates and clusters) according
    to the selected options 'minsize', 'minseqlength', 'alpha parameter' and 'identity', options 1, 1a, 1b, 1c and 1d
    """
    _sync_core_globals(core)
    if rmenu == "1" or rmenu == "1a":
        os.chdir(current_dir)
        print(f"\nComputing statistics after using option {rmenu}, please wait...\n" + successStyle +
              f"Results in --> {current_dir}{fileSep}outputs{fileSep}Stats_option_{rmenu}.txt" + normalStyle)
        with open(f"outputs/Stats_option_{rmenu}.txt", "w") as out:
            out.write(f"With option {rmenu}, parameters set to:\n\n"
                      f"Directory = {current_dir}\n"
                      f"Fastq R1 file name = {fastq_R1}\n"
                      f"Fastq R2 file name = {fastq_R2}\n"
                      f"Paired-end based loci = {lociPEs}\n"
                      f"Single-end based (R1) loci = {lociSEs}\n"
                      f"Samples = {samples}\n"
                      f"Minimum abundance for clusters = {minsize}\n"
                      f"Minimum length for sequences = {minseqlength}\n"
                      f"Alpha clustering parameter = {alpha}\n"
                      f"Identity for allocating clusters = {identity}\n")
            for sample in samples:
                sample = sample.rstrip()

                r1fa = open(f"tmp_files/{sample}_singleEnd.fa", "rt")
                reads = r1fa.read()
                nb_reads = reads.count(">")

                try:
                    derepfasr1 = open(f"tmp_files/{sample}_singleEnd_derep.fas", "rt")
                    dfr1 = derepfasr1.read()
                    nb_derpr1 = dfr1.count("sample")

                    clusfasr1 = open(f"tmp_files/{sample}_singleEnd_cluster.fas", "rt")
                    clusr1 = clusfasr1.read()
                    nb_clusr1 = clusr1.count("sample")

                    clusfasr1ok = open(f"tmp_files/{sample}_singleEnd_cluster_OK.fas", "rt")
                    clusr1ok = clusfasr1ok.read()
                    nb_clusr1ok = clusr1ok.count("sample")
                except FileNotFoundError:
                    nb_derpr1 = 0
                    nb_clusr1 = 0
                    nb_clusr1ok = 0

                try:
                    merged = open(f"tmp_files/{sample}_pairedEnd.fa", "rt")
                    mgd = merged.read()
                    nb_merged = mgd.count(">")
                    percent_merging = (nb_merged / nb_reads) * 100
                    percent = '{:.2f}%'.format(percent_merging)

                    derepfas = open(f"tmp_files/{sample}_pairedEnd_derep.fas", "rt")
                    df = derepfas.read()
                    nb_derpm = df.count("sample")

                    clusmerged = open(f"tmp_files/{sample}_pairedEnd_cluster.fas", "rt")
                    clusmer = clusmerged.read()
                    nb_clusm = clusmer.count("sample")

                    clusmergedok = open(f"tmp_files/{sample}_pairedEnd_cluster_OK.fas", "rt")
                    clusmerok = clusmergedok.read()
                    nb_clusmok = clusmerok.count("sample")
                except FileNotFoundError:
                    nb_merged = 0
                    percent = 0
                    nb_derpm = 0
                    nb_clusm = 0
                    nb_clusmok = 0

                out.writelines(f"\nSample {sample} has:\n"
                               f"\t{nb_reads} reads\n"
                               f"\t{nb_merged} merged sequences\n"
                               f"\tThe percentage of merging is {percent}\n"
                               f"\t{nb_derpm} dereplicated merged sequences\n"
                               f"\t{nb_clusm} merged clusters\n"
                               f"\t{nb_clusmok} merged clusters without chimera (OK)\n\n"
                               f"\t{nb_derpr1} dereplicated R1 sequences\n"
                               f"\t{nb_clusr1} R1 clusters\n"
                               f"\t{nb_clusr1ok} R1 clusters without chimera (R1_OK)\n\n")
                for locusPE in lociPEs:
                    os.chdir(f"results_by_locus/{locusPE}")
                    refs = open(f"{sample}_pairedEnd.fas")
                    refsloc = refs.read()
                    nb_ref = refsloc.count("sample")
                    out.writelines(f"\t{nb_ref} clusters of merged sequences of {sample} affiliated "
                                   f"to locus {locusPE}\n")
                    os.chdir(current_dir)
                for locusSE in lociSEs:
                    os.chdir(f"results_by_locus/{locusSE}")
                    refs2 = open(sample + "_singleEnd.fas")
                    refsloc2 = refs2.read()
                    nb_ref2 = refsloc2.count("sample")
                    out.writelines(f"\t{nb_ref2} clusters of single-end (R1) sequences of {sample} affiliated "
                                   f"to locus {locusSE}\n")
                    os.chdir(current_dir)
            print(successStyle + "\n" + ("The MAIN MANDATORY ANALYSIS option 1" if rmenu == "1" else "Execution of option 1a") + " is complete\n" + normalStyle)
        out.close()

    elif rmenu == "1b":
        os.chdir(current_dir)
        print("\nComputing statistics after using option 1b, please wait...\n"
              f"Results in --> {current_dir}{fileSep}outputs{fileSep}Stats_option_1b.txt:")
        with open("outputs/Stats_option_1b.txt", "w") as out:
            out.write("With option 1b, parameters set to:\n\n"
                      f"Directory = {current_dir}\n"
                      f"Fastq R1 file name = {fastq_R1}\n"
                      f"Fastq R2 file name = {fastq_R2}\n"
                      f"Paired-end based loci = {lociPEs}\n"
                      f"Single-end based (R1) loci = {lociSEs}\n"
                      f"Samples = {samples}\n"
                      f"Minimum abundance for clusters = {minsize}\n"
                      f"Minimum length for sequences = {minseqlength}\n"
                      f"Alpha clustering parameter = {alpha}\n"
                      f"Identity for allocating clusters = {identity}\n\n"
                      f"The selected locus is {loc_sel1}\n"
                      f"For {loc_sel1}:\n")
            os.chdir(f"./results_by_locus/{loc_sel1}")
            for sample in samples:
                a = open(sample + "_pairedEnd.fas")
                b = a.read()
                c = b.count("sample")
                out.write(f"\t{sample} has {c} clusters\n")
            os.chdir(current_dir)
        print(successStyle + "\nExecution of option 1b is complete\n" + normalStyle)
        out.close()

    elif rmenu == "1c":
        os.chdir(current_dir)
        print("\nComputing statistics after using option 1c, please wait...\n"
              f"Results in --> {current_dir}{fileSep}outputs{fileSep}Stats_option_1c.txt")
        with open("outputs/Stats_option_1c.txt", "w") as out:
            out.write("With option 1c, parameters set to:\n\n"
                      f"Directory = {current_dir}\n"
                      f"Fastq R1 file name = {fastq_R1}\n"
                      f"Fastq R2 file name = {fastq_R2}\n"
                      f"Paired-end based loci = {lociPEs}\n"
                      f"Single-end based (R1) loci = {lociSEs}\n"
                      f"Samples = {samples}\n"
                      f"Minimum abundance for clusters = {minsize}\n"
                      f"Minimum length for sequences = {minseqlength}\n"
                      f"Alpha clustering parameter = {alpha}\n"
                      f"Identity for allocating clusters = {identity}\n\n"
                      f"The selected locus is {loc_sel2}\n"
                      f"For {loc_sel2}:\n")
            os.chdir(f"./results_by_locus/{loc_sel2}")
            for sample in samples:
                a = open(sample + "_singleEnd.fas")
                b = a.read()
                c = b.count("sample")
                out.write(f"\t{sample} has {c} clusters\n")
            os.chdir(current_dir)
        print(successStyle + "\nExecution of option 1c is complete\n" + normalStyle)
        out.close()

    elif rmenu == "1d":
        os.chdir(current_dir)
        print("\nComputing statistics after using option 1d, wait\n"
              f"Results in --> {current_dir}{fileSep}outputs{fileSep}Stats_option_1d.txt:")
        with open("outputs/Stats_option_1d.txt", "w") as out:
            out.write("With option 1d, parameters set to:\n\n"
                      f"Directory = {current_dir}\n"
                      f"Fastq R1 file name = {fastq_R1}\n"
                      f"Fastq R2 file name = {fastq_R2}\n"
                      f"Paired-end based loci = {lociPEs}\n"
                      f"Single-end based (R1) loci = {lociSEs}\n"
                      f"Samples = {samples}\n"
                      f"Minimum abundance for clusters = {minsize}\n"
                      f"Minimum length for sequences = {minseqlength}\n"
                      f"Alpha clustering parameter = {alpha}\n"
                      f"Identity for allocating clusters = {identity}\n\n"
                      f"The selected sample is {sam_sel}\n")

            for locusPE in lociPEs:
                os.chdir(f"./results_by_locus/{locusPE}")
                a = open(sam_sel + "_pairedEnd.fas")
                b = a.read()
                c = b.count("sample")
                out.write(f"\t{sam_sel} has {c} clusters for locus {locusPE} based on paired-end reads\n")
                os.chdir(current_dir)

            for locusSE in lociSEs:
                os.chdir(f"./results_by_locus/{locusSE}")
                a = open(sam_sel + "_singleEnd.fas")
                b = a.read()
                c = b.count("sample")
                out.write(f"\t{sam_sel} has {c} clusters for locus {locusSE} based on R1/single-end reads\n")
                os.chdir(current_dir)
        print(successStyle + "\nExecution of option 1d is complete\n" + normalStyle)
        out.close()
