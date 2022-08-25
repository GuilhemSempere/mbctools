#!/usr/bin/python3

"""This Python program is a toolkit to make the use of VSEARCH easier and interactive, to analyze your
metabarcoding NGS data in the best conditions. It proposes the following MAIN MENU:

0  -> Initial and optional quality checking of your fastq files (slow - be patient!)
1  -> NEW COMPLETE ANALYSIS
1a -> Re-analyze all loci, from the clustering step, modifying parameters
1b -> Re-analyze only one locus of paired-end amplicon (merged reads), modifying parameters
1c -> Re-analyze only one locus of single-end amplicon (R1 only), modifying parameters
1d -> Re-analyse only one sample, modifying parameters
2  -> SELECTION OF MINIMUM SIZES ACCORDING TO THRESHOLDS SET BY THE USER
3  -> CLUSTERING OF ALL SAMPLES BY LOCUS FOR PHYLOGENETIC ANALYZES
end -> Exit the program

VSEARCH reference:
Rognes T, Flouri T, Nichols B, Quince C, Mahe F (2016)
VSEARCH: a versatile open source tool for metagenomics
PeerJ 4:e2584 doi: 10.7717/peerj.2584 https://doi.org/10.7717/peerj.2584

Usage:
=====
    python mbctools.py
"""

__authors__ = "Christian Barnabe"
__contact__ = "christian.barnabe@ird.fr"
__date__ = "2022"
__version__ = "0.1"
__copyright__ = " copyleft "

import sys
import datetime
from pathlib import Path
import os
import platform
import subprocess
import re
import glob


# BEFORE ANYTHING ELSE: MAKE SURE VSEARCH IS INSTALLED ################################################################################################
try:
    p = subprocess.run(["vsearch"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
except FileNotFoundError:
    print("Please install vsearch, then try again...")
    exit(1)


global dir_fastq, fastqr1_user, fastqr2_user, loci1_user, loci2_user, sample_user, minsize_user, minseqlength, alpha
global identity, loc_sel, fastqr1s, fastqr2s, sam_sel, loci1s, loci2s, samples

global shellCmd, scriptExt, fileSep
winOS = "Windows" in platform.uname()
shellCmd = "powershell" if winOS else "bash"
scriptExt = "ps1" if winOS else "sh"
fileSep = "\\" if winOS else "/"
globalErrorOnStopCmd = "" if winOS else "set -e"
localErrorOnStopCmd = "; If ($LASTEXITCODE -gt 0) { exit $LASTEXITCODE }" if winOS else ""


# FUNCTIONS FOR MULTIPLATFORM SUPPORT #################################################################################
def start_log_redirect(filepath):
    if not winOS:
        return "exec 3>&1 4>&2 >" + filepath + " 2>&1\n"
    else:
        return "&{\n"


def end_log_redirect():
    if not winOS:
        return "exec 1>&3 2>&4\n"
    else:
        return "} 2> '../infiles/results.log'\n"


def main_stream_message(message):
    if not winOS:
        return "printf \"" + message + "\" >&3\n"
    else:
        return "Write-Host -NoNewline \"" + message + "\"\n"


# FOLDER CREATION #####################################################################################################
def folder_infiles():
    """Creates folder 'infiles'
    """
    sys.stdout.write("")
    folder = "infiles"
    path = os.path.join(current_dir, folder)
    if Path(path).is_dir():
        sys.stdout.write(f"\nThe folder {folder} already exists and will be used in the current analysis")
    else:
        os.mkdir(path)


def folder_scripts():
    """Creates folder 'scripts'
    """
    sys.stdout.write("")
    folder = "scripts"
    path = os.path.join(current_dir, folder)
    if Path(path).is_dir():
        sys.stdout.write(f"\nThe folder {folder} already exists and will be used in the current analysis")
    else:
        os.mkdir(path)


def folder_loci():
    sys.stdout.write("")
    for loci1 in loci1s:
        path = os.path.join(current_dir, loci1)
        if Path(path).is_dir():
            sys.stdout.write(f"\nThe folder {loci1} already exists and will be used in the current analysis")
        else:
            os.mkdir(path)
    for loci2 in loci2s:
        path = os.path.join(current_dir, loci2)
        if Path(path).is_dir():
            sys.stdout.write(f"\nThe R1 folder {loci2} already exists and will be used in the current analysis")
        else:
            os.mkdir(path)


def folder_refs():
    sys.stdout.write("")
    folder = "refs"
    path = os.path.join(current_dir, folder)
    if Path(path).is_dir():
        sys.stdout.write(f"\nThe folder {folder} already exists and will be used in the current analysis")
    else:
        os.mkdir(path)


# INPUTS AND VARIABLES ################################################################################################
def input_dir_fastq():
    global dir_fastq
    dir_fastq = input("\nEnter the full path to the folder where fastq files are located"
                      "\ne.g. C:/Users/barnabe/Documents/NGSdata/run190823"
                      "\ndefault = <current_directory>/fastq: ")
    while Path(dir_fastq).is_dir() is False and dir_fastq != '':
        dir_fastq = input("\nYour path is not valid. Please enter a valid name for the path\n"
                          "default = <current_directory>/fastq: ")
    else:
        dir_fastq = current_dir + "/fastq"


def input_fastqr1_user():
    global fastqr1_user, fastqr1s
    fastqr1_user = input("\nEnter the name of your fastqR1 file, default = fastqR1.txt: ")
    while os.path.isfile(fastqr1_user) is False and fastqr1_user != '':
        fastqr1_user = input("\nYour file name is not valid, please enter a valid name\n"
                             "default = fastqR1.txt: ")
    else:
        fastqr1_user = 'fastqR1.txt'
    with open(fastqr1_user, "r") as out1:
        fastqr1s = out1.read().splitlines()


def input_fastqr2_user():
    global fastqr2_user, fastqr2s
    fastqr2_user = input("\nEnter the name of your fastqR2 file, default = fastqR2.txt: ")
    while os.path.isfile(fastqr2_user) is False and fastqr2_user != '':
        fastqr2_user = input("\nYour file name is not valid, please enter a valid name\n"
                             "default = fastqR2.txt: ")
    else:
        fastqr2_user = 'fastqR2.txt'
    with open(fastqr2_user, "r") as out2:
        fastqr2s = out2.read().splitlines()


def input_loci1_user():
    global loci1_user, loci1s
    loci1_user = input("\nEnter the name of your file containing the loci based on paired-end reads\n"
                       "default = locus1.txt: ")
    while os.path.isfile(loci1_user) is False and loci1_user != '':
        loci1_user = input("\nYour file name is not valid, please enter a valid name: ")
    else:
        loci1_user = 'locus1.txt'
    with open(loci1_user, "r") as out3:
        loci1s = out3.read().splitlines()


def input_loci2_user():
    global loci2_user, loci2s
    loci2_user = input("\nEnter the name of your file containing the loci based on single-end reads only (R1),\n"
                       "default = locus2.txt: ")
    while os.path.isfile(loci2_user) is False and loci2_user != '':
        loci2_user = input("\nYour name is not valid, please enter a valid name: ")
    else:
        loci2_user = 'locus2.txt'

    with open(loci2_user, "r") as out4:
        loci2s = out4.read().splitlines()


def input_sample_user():
    global sample_user
    sample_user = input("\nEnter the name of your file containing the sample names, default = samples.txt: ")
    while os.path.isfile(sample_user) is False and sample_user != '':
        sample_user = input("\nYour file name is not valid, please enter a valid name\n"
                            "default = samples.txt: ")
    else:
        sample_user = 'samples.txt'
    global samples
    with open(sample_user, "r") as out5:
        samples = out5.read().splitlines()


def input_minsize_user():
    global minsize_user
    minsize_user = input("\nEnter the minsize option value for your clusters,\n"
                         "i.e. the minimum abundance of the retained clusters, default = 8: ")
    while minsize_user.isnumeric() is False and minsize_user != '':
        minsize_user = input("\nYour minsize option may be an integer, e.g. 2, 10, 50... default = 8: : ")
    if minsize_user == '':
        minsize_user = '8'


def input_minseqlength():
    global minseqlength
    minseqlength = input("\nEnter the minimum length of sequences to keep for any locus, default = 100: ")
    while minseqlength.isnumeric() is False and minseqlength != '':
        minseqlength = input("\nYour minimum length may be an integer, e.g. 100, 150, 180... default = 100: ")
    if minseqlength == '':
        minseqlength = '100'


def input_alpha():
    global alpha
    alpha = input("\nEnter your alpha parameter for the clustering, default = 2: ")
    while alpha.isnumeric() is False and alpha != '':
        alpha = input("\nYour alpha parameter may be an integer, e.g. 0, 1, 2, 3... default = 2: ")
    if alpha == '':
        alpha = '2'


def input_identity():
    global identity
    identity = input("\nEnter your identity parameter to BLAST the clusters against the references, default = 0.7: ")
    while identity > '1' and identity != '':
        identity = input("\nYour identity parameter may be a decimal < 1, e.g. 0.85, 0.75, 0.95... default = 0.7: ")
    if identity == '':
        identity = '0.7'


def input_loc_sel_merged():
    global loc_sel
    loc_sel = input("\nEnter the name of the locus analysed by paired-end reads you want to rerun, no default: ")
    while loc_sel not in loci1s:
        loc_sel = input("\nYour locus name is not valid. Please enter the valid name of the locus: ")


def input_loc_sel_r1():
    global loc_sel
    loc_sel = input("\nEnter the name of the locus analysed by only single-end (R1) reads "
                    "you want to rerun, no default: ")
    while loc_sel not in loci2s:
        loc_sel = input("\nYour locus name is not valid. Please enter the valid name of the locus: ")


def input_sam_sel():
    global sam_sel, samples
    sam_sel = input("\nEnter the sample name you want to rerun, no default: ")
    while sam_sel not in samples:
        sam_sel = input("\nThe sample name is not valid, please enter a valid sample name, no default: ")


# PARAMETER FILES #####################################################################################################
def param_1():
    """ Creates a file with one parameter by line for options 1
    """
    with open("my_parameters_option_1.txt", "w") as outfile:
        outfile.writelines(f"Run option 1: {date}\n{dir_fastq}\n{fastqr1_user}\n{fastqr2_user}\n{loci1_user}\n"
                           f"{loci2_user}\n{sample_user}\n{minsize_user}\n{minseqlength}\n{alpha}\n{identity}\n"
                           f"Samples used = {samples}\nLoci paired-end used = {loci1s}\nLoci single-end (R1) "
                           f"used = {loci2s}\n"
                           f"-------------------------------------------------------------\n")


def prev_param():
    """ Recalls global variables for different options
    """
    with open("my_parameters_option_1.txt", "r") as infile:
        lines = infile.read().splitlines()
        global dir_fastq
        dir_fastq = lines[1]
        global fastqr1_user
        fastqr1_user = lines[2]
        global fastqr2_user
        fastqr2_user = lines[3]
        global loci1_user
        loci1_user = lines[4]
        global loci2_user
        loci2_user = lines[5]
        global sample_user
        # sample_user = lines[6]
        sample_user = re.search("File of sample names:(.*)$", lines[6]).group(1)

    global fastqr1s
    with open(fastqr1_user, "r") as out1:
        fastqr1s = out1.read().splitlines()

    global fastqr2s
    with open(fastqr2_user, "r") as out1:
        fastqr2s = out1.read().splitlines()

    global loci1s
    with open(loci1_user, "r") as out1:
        loci1s = out1.read().splitlines()

    global loci2s
    with open(loci2_user, "r") as out1:
        loci2s = out1.read().splitlines()

    global samples
    with open(sample_user, "r") as out1:
        samples = out1.read().splitlines()


def param_1a():
    """ Creates a file with one parameter by line for option 1a
    """
    with open("my_parameters_option_1.txt", "a") as outfile:
        outfile.writelines(f"Run option 1a: {date}\n{dir_fastq}\n{fastqr1_user}\n{fastqr2_user}\n"
                           f"{loci1_user}\n{loci2_user}\n{sample_user}\n{minsize_user}\n{minseqlength}"
                           f"\n{alpha}\n{identity}\n"
                           f"-------------------------------------------------------------\n")


def param_1b():
    """Creates a file with one parameter by line for option 1b
    """
    with open("my_parameters_option_1.txt", "a") as outfile:
        outfile.writelines(f"Run option 1b: {date}\n{dir_fastq}\n{fastqr1_user}\n{fastqr2_user}\n{loc_sel}\n"
                           f"{sample_user}\n{minsize_user}\n{minseqlength}\n{alpha}\n{identity}\n"
                           f"-------------------------------------------------------------\n")


def param_1c():
    """Creates a file with one parameter by line for option 1c
    """
    with open("my_parameters_option_1.txt", "a") as outfile:
        outfile.writelines(f"Run option 1c: {date}\n{dir_fastq}\n{fastqr1_user}\n{fastqr2_user}\n{loc_sel}\n"
                           f"{sample_user}\n{minsize_user}\n{minseqlength}\n{alpha}\n{identity}\n"
                           f"-------------------------------------------------------------\n")


def param_1d():
    """Creates a file with one parameter by line for option 1d
    """
    with open("my_parameters_option_1.txt", "a") as outfile:
        outfile.writelines(f"Run option 1d: {date}\n{dir_fastq}\n{fastqr1_user}\n{fastqr2_user}\n{sam_sel}\n"
                           f"{minsize_user}\n{minseqlength}\n{alpha}\n{identity}\n"
                           f"-------------------------------------------------------------\n")


# STATISTICS ON FASTQ FILES ###########################################################################################
def stat():
    """Tests the quality of each 'fastq' file by the VSEARCH command:
    vsearch --fastq_eestats2 ../infiles/fastqF-R1 --output ../infiles/SN_R1info.txt
    vsearch --fastq_eestats2 ../infiles/fastqF-R2 --output ../infiles/SN_R2info.txt

    fastqF = fastq file name
    SN = sample name
    """
    with open("scripts/infor1." + scriptExt, "w") as out3:
        i = 0
        out3.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Checking R1 data quality:'))
        while i < len(samples):
            sample = samples[i]
            fastqr1 = fastqr1s[i]
            i = i + 1
            out3.write(main_stream_message(f' {sample}...'))
            out3.write(f"vsearch --fastq_eestats2 {dir_fastq}{fileSep}{fastqr1} --output "
                       f"../infiles/{sample}_R1info.txt" + localErrorOnStopCmd + "\n")
        out3.write(main_stream_message(f'\n\n'))

    with open("scripts/infor2." + scriptExt, "w") as out4:
        i = 0
        out4.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Checking R2 data quality:'))
        while i < len(samples):
            sample = samples[i]
            fastqr2 = fastqr2s[i]
            i = i + 1
            out4.write(main_stream_message(f' {sample}...'))
            out4.write(f"vsearch --fastq_eestats2 {dir_fastq}{fileSep}{fastqr2} --output "
                       f"../infiles/{sample}_R2info.txt" + localErrorOnStopCmd + "\n")
        out4.write(main_stream_message(f'\n\n'))


# STEP 1: MERGING FOR AMPLICONS BASED ON PAIRED-END READS #############################################################
def merging():
    """Merges paired-end reads into one sequence, when the length of the expected amplicon allows it
    according the VSEARCH command:
    vsearch --fastq_mergepairs ../infiles/fastqR1 --reverse ../infiles/fastqR2 --fastaout
    ./infiles/SN.fa --fastq_allowmergestagger --relabel sample='sample=SN'_merged

    fastqR1 = fastqR1 complete name
    fastqR2 = fastqR2 complete name
    option --fastq_allowmergestagger allows the merging of short fragments
    """
    with open("scripts/merging." + scriptExt, "w") as out5:
        i = 0
        out5.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Merging paired-end reads:'))
        while i < len(samples):
            sample = samples[i]
            fastqr1 = fastqr1s[i]
            fastqr2 = fastqr2s[i]
            out5.write(main_stream_message(f' {sample}...') +
                       f"vsearch --fastq_mergepairs {dir_fastq}{fileSep}{fastqr1} "
                       f"--reverse {dir_fastq}{fileSep}{fastqr2} --fastaout "
                       f"../infiles/{sample}.fa --fastq_allowmergestagger --fastq_maxee 1 --relabel "
                       f"sample={sample}_merged" + localErrorOnStopCmd + "\n")
            i = i + 1
        out5.write(main_stream_message(f'\n\n'))


# FASTQ TO FASTA FOR AMPLICONS BASED ON SINGLE-END READS (R1 only) #####################################################
def fastq2fas():
    """When the merging R1/R2 is impossible because of an unadapted size of amplicon, the reads R1 of 301 bp
    (better than R2) are used to search the relevant sequences.
    First, all R1 'fastq' files have to be transformed into 'fasta' files by the VSEARCH command:
    vsearch --fastq_filter ../infiles/fastaqR1 --fastaout ../infiles/SN_R1.fa

    SN = sample name
    fastqR1 = fastqR1 file name
    """
    with open("scripts/fqtofas." + scriptExt, "w") as out6:
        i = 0
        out6.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Converting FASTQ files into FASTA format:'))
        while i < len(samples):
            sample = samples[i]
            fastqr1 = fastqr1s[i]
            out6.write(main_stream_message(f' {sample}...') +
                       f"vsearch --fastq_filter {dir_fastq}{fileSep}{fastqr1} --fastaout ../infiles/{sample}_R1.fa" + localErrorOnStopCmd + "\n")
            i = i + 1
        out6.write(main_stream_message(f'\n\n'))


# DEREPLICATION #######################################################################################################
def derep_merged():
    """Dereplicates merged sequences in a given 'fasta' file with the VSEARCH command:
    vsearch --derep_fulllength ../infiles/SN.fa --output ../infiles/SN_derep.fas --sizeout
    --strand both

    SN = sample name
    Dereplicates in both strands and writes abundance annotation (frequency) to output.
    """
    with open("scripts/derep." + scriptExt, "w") as out7:
        out7.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Dereplicating paired-end reads:'))
        for sample in samples:
            out7.write(main_stream_message(f' {sample}...') +
                       f"vsearch --derep_fulllength ../infiles/{sample}.fa --output ../infiles/{sample}_derep.fas "
                       f"--sizeout --strand both" + localErrorOnStopCmd + "\n")
        out7.write(main_stream_message(f'\n\n'))


def derep_r1():
    """Dereplicates R1 sequences in a given 'fasta' file with the VSEARCH command:
    vsearch --derep_fulllength ../infiles/_R1_SN'.fa --output ../infiles/SN_derep.fas --sizeout
    --strand both --relabel sample='SN'_merged.

    SN = sample name
    """
    with open("scripts/derep_r1." + scriptExt, "w") as out8:
        out8.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Dereplicating single-end reads:'))
        for sample in samples:
            out8.write(main_stream_message(f' {sample}...') +
                       f"vsearch --derep_fulllength ../infiles/{sample}_R1.fa --output ../infiles/{sample}"
                       f"_derep_R1.fas --sizeout --strand both --relabel sample={sample}_R1." + localErrorOnStopCmd + "\n")
        out8.write(main_stream_message(f'\n\n'))


# UNOISING AND CLUSTERING #############################################################################################
def cluster_merged():
    """Denoises and clusters Illumina dereplicated merged sequences and gives in output the centroids sequences
    to 'fasta' files with the VSEARCH command:
    vsearch --cluster_unoise ../infiles/SN_derep.fas --sizein --centroids ../infiles/SN_cluster.fas
    --strand both --minsize int --sizeout --sizeorder --unoise_alph int --minseqlength int

    SN = sample name
    int = integer

    NB: It's the crucial step of metabarcoding analyze. Three options are used 'minsize' - all sequences of frequency
    < minsize will be discarded; 'alpha parameter' see doc VSEARCH for interpretation and 'minseqlength' - all
    sequences of length < minseqlenght will be discarded. This step gives all the relevant sequences for both
     frequency and length criteria.
    """
    with open("scripts/cluster." + scriptExt, "w") as out9:
        i = 0
        out9.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Clustering paired-end reads:'))
        while i < len(samples):
            sample = samples[i]
            out9.write(main_stream_message(f' {sample}...') +
                       f"vsearch --cluster_unoise ../infiles/{sample}_derep.fas --sizein --centroids "
                       f"../infiles/{sample}_cluster.fas --strand both --minsize {minsize_user} --sizeout --sizeorder "
                       f"--unoise_alph {alpha} --minseqlength {minseqlength}" + localErrorOnStopCmd + "\n")
            i = i + 1
        out9.write(main_stream_message(f'\n\n'))


def cluster_r1():
    """Denoises and clusters Illumina dereplicated R1 sequences and gives in output the centroids sequences
    to 'fasta' files with the VSEARCH command:
    vsearch --cluster_unoise ../infiles/SN_derep.fas --sizein --centroids ../infiles/SN_cluster.fas
    --strand both --minsize int --sizeout --sizeorder --unoise_alph int --minseqlength int

    SN = sample name
    int = integer
    Similar to cluster_merged but for R1 sequences
    """
    with open("scripts/cluster_r1." + scriptExt, "w") as out10:
        i = 0
        out10.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Clustering single-end reads:'))
        while i < len(samples):
            sample = samples[i]
            out10.write(main_stream_message(f' {sample}...') +
                        f"vsearch --cluster_unoise ../infiles/{sample}_derep_R1.fas --sizein --centroids "
                        f"../infiles/{sample}_cluster_R1.fas --strand both --minsize {minsize_user} --sizeout "
                        f"--sizeorder --unoise_alpha {alpha} --minseqlength {minseqlength}" + localErrorOnStopCmd + "\n")
            i = i + 1
        out10.write(main_stream_message(f'\n\n'))


def cluster_one_sample_1d():
    """Denoises and clusters Illumina dereplicated for a selected sample and gives in output the centroids sequences
    to 'fasta' files, for a selected sample, with the VSEARCH command:
    vsearch --cluster_unoise ../infiles/SS_derep.fas --sizein --centroids ../infiles/SS_cluster.fas
    --strand both --minsize int --sizeout --sizeorder --unoise_alph int --minseqlength int

    SS = selected sample
    """
    with open("scripts/cluster_one_sample_1d." + scriptExt, "w") as out1:
        out1.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Clustering reads for sample {sam_sel}:') +
                   f"vsearch --cluster_unoise ../infiles/{sam_sel}_derep.fas --sizein --centroids "
                   f"../infiles/{sam_sel}_cluster.fas --strand both --minsize {minsize_user} --sizeout "
                   f"--sizeorder --unoise_alph {alpha} --minseqlength {minseqlength}" + localErrorOnStopCmd + "\n"
                   f"vsearch --cluster_unoise ../infiles/{sam_sel}_derep_R1.fas --sizein --centroids "
                   f"../infiles/{sam_sel}_cluster_R1.fas --strand both --minsize {minsize_user} "
                   f"--sizeout --sizeorder --unoise_alph {alpha} --minseqlength {minseqlength}" + localErrorOnStopCmd + "\n")
        out1.write(main_stream_message(f'\n\n'))


# REMOVING CHIMERA ####################################################################################################
def chimera_merged():
    """Detects and removes potential chimeras in denoised merged sequences by the VSEARCH command:
    vsearch --uchime3_denovo ../infiles/SN_cluster.fas --nonchimeras ../infiles/SN_cluster_OK.fas

    SN = sample name
    After denoising, the chimeras are very scarce.
    """
    with open("scripts/chimera." + scriptExt, "w") as out11:
        out11.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Detecting chimeras within paired-end reads:'))
        for sample in samples:
            out11.write(main_stream_message(f' {sample}...') +
                        f"vsearch --uchime3_denovo ../infiles/{sample}_cluster.fas --nonchimeras "
                        f"../infiles/{sample}_cluster_OK.fas" + localErrorOnStopCmd + "\n")
        out11.write(main_stream_message(f'\n\n'))


def chimera_r1():
    """Detects and removes potential chimeras in R1 sequences by the VSEARCH command:
    vsearch --uchime3_denovo ../infiles/SN_cluster-R1.fas --nonchimeras ../infiles/SN_cluster_R1_OK.fas

    SN = sample name
    After denoising, the chimeras are very scarce.
    """
    with open("scripts/chimera_r1." + scriptExt, "w") as out12:
        out12.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Detecting chimeras within single-end reads:'))
        for sample in samples:
            out12.write(main_stream_message(f' {sample}...') +
                        f"vsearch --uchime3_denovo ../infiles/{sample}_cluster_R1.fas --nonchimeras"
                        f" ../infiles/{sample}_cluster_R1_OK.fas" + localErrorOnStopCmd + "\n")
        out12.write(main_stream_message(f'\n\n'))


def chimera_one_sample_1d():
    """Detects potential chimeras in sequences for a selected sample by the VSEARCH command:
    vsearch --uchime3_denovo ../infiles/SS_cluster.fas --nonchimeras ../infiles/SS_cluster_OK.fas
    vsearch --uchime3_denovo ../infiles/SS_cluster_R1.fas --nonchimeras ../infiles/SS_cluster_R1_OK.fas

    ss = selected sample
    After denoising, the chimeras are very scarce
    """
    with open("scripts/chimera_one_sample_1d." + scriptExt, "w") as out1:
        out1.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Detecting chimeras for sample {sam_sel}') +
                   f"vsearch --uchime3_denovo ../infiles/{sam_sel}_cluster.fas --nonchimeras ../infiles/{sam_sel}"
                   f"_cluster_OK.fas" + localErrorOnStopCmd + "\n"
                   f"vsearch --uchime3_denovo ../infiles/{sam_sel}_cluster_R1.fas --nonchimeras "
                   f"../infiles/{sam_sel}_cluster_R1_OK.fas" + localErrorOnStopCmd + "\n")
        out1.write(main_stream_message(f'\n\n'))


# DISTRIBUTION OF CLUSTERS INTO THE DIFFERENT LOCI ####################################################################
def runloc_merged():
    """Searches similarities between merged, denoised and non-chimera sequences and the local reference
    database (-db) by the VSEARCH command:
    vsearch --usearch_global ../infiles/SN_cluster_OK.fas --db ../refs/L1.fas --matched ../L1/SN_merged.fas
    --id int --strand both

    L1 = locus name for amplicons based on paired-end reads
    This is the second crucial step metabarcoding analyze, depending on a required minimal identity setup by
    the option 'identity'.
    """
    with open("scripts/locimerged." + scriptExt, "w") as out13:
        out13.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters to loci for paired-end reads:'))
        for loci1b in loci1s:
            for sample in samples:
                out13.write(main_stream_message(f' {sample} VS {loci1b}...') +
                            f"vsearch --usearch_global ../infiles/{sample}_cluster_OK.fas --db ../refs/{loci1b}.fas"
                            f" --matched ../{loci1b}/{sample}_merged.fas --id {identity} --strand both" + localErrorOnStopCmd + "\n")
        out13.write(main_stream_message(f'\n\n'))


def runloc_r1():
    """Similar to runloc_merged but with the R1 denoised and non-chimera sequences
    in case of amplicons where the merging R1/R2 is impossible, with the VSEARCH command:
    vsearch --usearch_global ../infiles/SN_cluster_R1_OK.fas --db ../refs/L2.fas --matched ../L2/SN_R1.fas
     --id real --strand both

     L1 = locus for amplicons with no mergeable R1/R2
     SN = sample name
     real = a real from 0 to 1, generally around 0.7
    """
    with open("scripts/locir1." + scriptExt, "w") as out14:
        out14.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters to loci for single-end reads:'))
        for locus2b in loci2s:
            for sample in samples:
                out14.write(main_stream_message(f' {sample} VS {locus2b}...') +
                            f"vsearch --usearch_global ../infiles/{sample}_cluster_R1_OK.fas --db ../refs/{locus2b}.fas"
                            f" --matched ../{locus2b}/{sample}_R1.fas --id {identity} --strand both" + localErrorOnStopCmd + "\n")
        out14.write(main_stream_message(f'\n\n'))


# DISTRIBUTION OF CLUSTERS FOR ONE SELECTED LOCUS #####################################################################
def runlocsel_merged():
    """Similar to runloc_merged but with only one selected locus with the VSEARCH command:
    vsearch --usearch_global ../infiles/SN_cluster_OK.fas --db ../refs/L1.fas --matched
    ../L2/SN_merged.fas --id real --strand both

    SN = sample name
    L1 = locus for amplicons based on paired-end reads
    """
    with open("scripts/loci_sel." + scriptExt, "w") as out15:
        out15.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'locsel:'))
        for sample in samples:
            out15.write(main_stream_message(f' {sample}...') +
                        f"vsearch --usearch_global ../infiles/{sample}_cluster_OK.fas --db ../refs/{loc_sel}.fas"
                        f" --matched ../{loc_sel}/{sample}_merged.fas --id {identity} --strand both" + localErrorOnStopCmd + "\n")
        out15.write(main_stream_message(f'\n\n'))


def runlocsel_r1():
    """identical to runlocsel_merged but with R1 denoised and non-chimera sequences with the VSEARCH command:
    vsearch --usearch_global ../infiles/SN_cluster_R1_OK.fas --db ../refs/L2.fas --matched
    ../L2/SN_R1.fas --id real --strand both

    SN = sample name
    L2 = locus name for amplicons with no mergeable R1/R2
    """
    with open("scripts/locir1_sel." + scriptExt, "w") as out16:
        out16.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'loc_selR1:'))
        for sample in samples:
            out16.write(main_stream_message(f' {sample}...') +
                        f"vsearch --usearch_global ../infiles/{sample}_cluster_R1_OK.fas --db ../refs/{loc_sel}.fas"
                        f" --matched ../{loc_sel}/{sample}_R1.fas --id {identity} --strand both" + localErrorOnStopCmd + "\n")
        out16.write(main_stream_message(f'\n\n'))


# DISTRIBUTION OF CLUSTERS FOR ONE SELECTED SAMPLE ####################################################################
def runloc_one_sample_1d():
    """ Searches similarities between denoised and non-chimera sequences and your local reference database (db)
    by the VSEARCH command, but only for a selected sample:
    vsearch --usearch_global ../infiles/SS_cluster_OK.fas --db ../refs/L1.fas --matched ../L1/SS_merged.fas
     --id real --strand both
    vsearch --usearch_global ../infiles/SS_cluster_OK.fas --db ../refs/L2.fas --matched ../L2/SS_merged.fas
     --id real --strand both

    SN = locus name
    L1 = locus name for amplicons based on paired-end reads
    L2 = locus name for amplicons with no mergeable R1/R2
    real = real from 0 to 1
    This is the second crucial step metabarcoding analyze, depending on a required minimal identity setup by
    the option 'identity'
    """
    with open("scripts/loci_merged_1d." + scriptExt, "w") as out13:
        out13.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'runloci:'))
        for loci1b in loci1s:
            out13.write(main_stream_message(f' {sam_sel}...') +
                        f"vsearch --usearch_global ../infiles/{sam_sel}_cluster_OK.fas --db ../refs/{loci1b}.fas"
                        f" --matched ../{loci1b}/{sam_sel}_merged.fas --id {identity} --strand both" + localErrorOnStopCmd + "\n")
        out13.write(main_stream_message(f'\n\n'))

    with open("scripts/loci_R1_1d." + scriptExt, "w") as out14:
        out14.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'runlociR1:'))
        for locus2b in loci2s:
            out14.write(main_stream_message(f' {sam_sel}...') +
                        f"vsearch --usearch_global ../infiles/{sam_sel}_cluster_R1_OK.fas --db ../refs/{locus2b}.fas"
                        f" --matched ../{locus2b}/{sam_sel}_R1.fas --id {identity} --strand both" + localErrorOnStopCmd + "\n")
            out14.write(main_stream_message(f'\n\n'))


# ORIENTATION OF SEQUENCES INTO FORWARD DIRECTION ####################################################################
def orient_merged():
    """Orients all the merged sequences in the same direction (forward) than references with the following script:
    vsearch --orient ../L1/SN_merged.fas --db ../refs/L1.fas --fastaout ../L1/SN_orient.fas

    SN = locus name
    L1 = locus name for amplicons based on paired-end reads
    """
    with open("scripts/orientloc1." + scriptExt, "w") as out17:
        out17.write(globalErrorOnStopCmd + "\n" + main_stream_message(f"Correcting paired-end reads' orientation:"))
        for locus1b in loci1s:
            for sample in samples:
                out17.write(main_stream_message(f' {sample} VS {locus1b}...') +
                            f"vsearch --orient ../{locus1b}/{sample}_merged.fas --db ../refs/{locus1b}.fas --fastaout "
                            f"../{locus1b}/{sample}_orient.fas" + localErrorOnStopCmd + "\n")
        out17.write(main_stream_message(f'\n\n'))


def orient_r1():
    """Orients all the R1 sequences in the same direction (forward) than references with the following script:
    vsearch --orient ../L2/SN_merged.fas --db ../refs/L2.fas --fastaout ../L2/SN_orient.fas

    SN = locus name
    L2 = locus name for amplicons with no mergeable R1/R2
    """
    with open("scripts/orientloc2." + scriptExt, "w") as out18:
        out18.write(globalErrorOnStopCmd + "\n" + main_stream_message(f"Correcting single-end reads' orientation:"))
        for locus2b in loci2s:
            for sample in samples:
                out18.write(main_stream_message(f' {sample} VS {locus2b}...') +
                            f"vsearch --orient ../{locus2b}/{sample}_R1.fas --db ../refs/{locus2b}.fas --fastaout "
                            f"../{locus2b}/{sample}_R1_orient.fas" + localErrorOnStopCmd + "\n")
        out18.write(main_stream_message(f'\n\n'))


def orient_one_loc_1b():
    """Orients all the merged sequences in the same direction (forward) for a selected locus based on paired-end
    R1/R2 mergeable reads with the following script:
    vsearch --orient ../SL/SN_merged.fas --db ../refs/SL.fas --fastaout ../SL/SN_orient.fas

    SL = selected locus
    SN = sample name
    """
    with open("scripts/orient_merged_1b." + scriptExt, "w") as out17:
        out17.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'orientR1:'))
        for sample in samples:
            out17.write(main_stream_message(f' {sample}...') +
                        f"vsearch --orient ../{loc_sel}/{sample}_merged.fas --db ../refs/{loc_sel}.fas --fastaout "
                        f"../{loc_sel}/{sample}_orient.fas" + localErrorOnStopCmd + "\n")
        out17.write(main_stream_message(f'\n\n'))


def orient_one_loc_1c():
    """Orients all the R1 sequences in the same direction (forward) for loci based on single-end R1 reads
    with the following script:
    vsearch --orient ../SL/SN_R1.fas --db ../refs/SL.fas --fastaout ../SL/SN_R1_orient.fas

    SL = selected locus
    SN = sample name
    """
    with open("scripts/orient_R1_1c." + scriptExt, "w") as out18:
        out18.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'orientR1:'))
        for sample in samples:
            out18.write(main_stream_message(f' {sample}...') +
                        f"vsearch --orient ../{loc_sel}/{sample}_R1.fas --db ../refs/{loc_sel}.fas --fastaout "
                        f"../{loc_sel}/{sample}_R1_orient.fas" + localErrorOnStopCmd + "\n")
        out18.write(main_stream_message(f'\n\n'))


def orient_one_sample_1d():
    """Orients all the merged and R1 sequences in the same direction (forward) for a selected sample with
    the following script:
    vsearch --orient ../L1/SS_merged.fas --db ../refs/L1.fas --fastaout ../L1/SS_orient.fas
    vsearch --orient ../L2/SS_merged.fas --db ../refs/L2.fas --fastaout ../L1/SS_orient.fas

    L1 = locus name for amplicons based on paired-end mergeable reads
    L2 = locus name for amplicons with no mergeable R1/R2 (R1 only)
    SS = selected sample
    """
    with open("scripts/orient_merged_1d." + scriptExt, "w") as out17:
        out17.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'loci1:'))
        for locus1b in loci1s:
            out17.write(main_stream_message(f' {sam_sel}...') +
                        f"vsearch --orient ../{locus1b}/{sam_sel}_merged.fas --db ../refs/{locus1b}.fas --fastaout "
                        f"../{locus1b}/{sam_sel}_orient.fas" + localErrorOnStopCmd + "\n")
        out17.write(main_stream_message(f'\n\n'))

    with open("scripts/orient_R1_1d." + scriptExt, "w") as out18:
        out18.write(main_stream_message(f'lociR1:'))
        for locus2b in loci2s:
            out18.write(main_stream_message(f' {sam_sel}...') +
                        f"vsearch --orient ../{locus2b}/{sam_sel}_merged.fas --db ../refs/{locus2b}.fas --fastaout "
                        f"../{locus2b}/{sam_sel}_R1_orient.fas" + localErrorOnStopCmd + "\n")
        out18.write(main_stream_message(f'\n\n'))


# CREATION OF GLOBAL SCRIPTS ##########################################################################################
def runall_0():
    """ Defines the run order for option 1
    """
    with open("scripts/runall0." + scriptExt, "w") as out19:
        if not winOS:
            out19.write(start_log_redirect('../infiles/results.log') + 
            "scriptArray=('./infor1." + scriptExt + "' './infor2." + scriptExt + "')\n" +
            'for script in "${scriptArray[@]}"\n' +
            '    do\n' +
            '    if ! ${script}; then\n' +
            '        printf "Error executing ${script}\\n\\n" >&3\n' +
            '        exit 1\n' +
            '    fi\n' +
            'done\n' +
            end_log_redirect())
        else:
            out19.write(start_log_redirect('../infiles/results.log') +
            '   $scriptArray = @("./infor1.' + scriptExt + '", "./infor2.' + scriptExt + '")\n' + 
            '   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
            '        $script = $scriptArray[$i]\n' +
            '        & "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; exit $LASTEXITCODE }\n' +
            '   }\n' +
            end_log_redirect())

    os.chdir('scripts')
    if not winOS:
        for file in os.listdir("."):
            os.chmod(file, 0o755)

    sys.stdout.write("Quality checking is being processed by 'vsearch', be patient!\n\n")
    return subprocess.run([shellCmd, "./runall0." + scriptExt]).returncode


def runall_1():
    """ Defines the run order for option 1
    """
    with open("scripts/runall1." + scriptExt, "w") as out19:
        if not winOS:
            out19.write(start_log_redirect('../infiles/results.log') + 
            "scriptArray=('./merging." + scriptExt + "' './fqtofas." + scriptExt + "' './derep." + scriptExt + "' './derep_r1." + scriptExt + "' './cluster." + scriptExt + "' './cluster_r1." + scriptExt + "' './chimera." + scriptExt + "' './chimera_r1." + scriptExt + "' './locimerged." + scriptExt + "' './locir1." + scriptExt + "' './orientloc1." + scriptExt + "' './orientloc2." + scriptExt + "')\n" +
            'for script in "${scriptArray[@]}"\n' +
            '    do\n' +
            '    if ! ${script}; then\n' +
            '        printf "Error executing ${script}\\n\\n" >&3\n' +
            '        exit 1\n' +
            '    fi\n' +
            'done\n' +
            end_log_redirect())
        else:
            out19.write(start_log_redirect('../infiles/results.log') +
            '   $scriptArray = @("./merging.' + scriptExt + '", "./fqtofas.' + scriptExt + '", "./derep.' + scriptExt + '", "./derep_r1.' + scriptExt + '", "./cluster.' + scriptExt + '", "./cluster_r1.' + scriptExt + '", "./chimera.' + scriptExt + '", "./chimera_r1.' + scriptExt + '", "./locimerged.' + scriptExt + '", "./locir1.' + scriptExt + '", "./orientloc1.' + scriptExt + '", "./orientloc2.' + scriptExt + '")\n' + 
            '   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
            '        $script = $scriptArray[$i]\n' +
            '        & "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; exit $LASTEXITCODE }\n' +
            '   }\n' +
            end_log_redirect())


    os.chdir('scripts')
    if not winOS:
        for file in os.listdir("."):
            os.chmod(file, 0o755)
    return subprocess.run([shellCmd, "./runall1." + scriptExt]).returncode


def runall_1a():
    """ Defines the run order for option 1a i.e. analysis of all loci from the clustering step.
    """
    with open("scripts/runall1a." + scriptExt, "w") as out19:
        out19.write(start_log_redirect('../infiles/res1a.log') +
                    "./cluster." + scriptExt + "\n"
                    "./cluster_r1." + scriptExt + "\n"
                    "./chimera." + scriptExt + "\n"
                    "./chimera_r1." + scriptExt + "\n"
                    "./locimerged." + scriptExt + "\n"
                    "./locir1." + scriptExt + "\n"
                    "./orientloc1." + scriptExt + "\n"
                    "./orientloc2." + scriptExt + "\n"
                    + end_log_redirect())
    os.chdir('scripts')
    subprocess.run([shellCmd, "./runall1a." + scriptExt])


def runall_1b():
    """ Defines the run order for option 1b
    """
    with open("scripts/runall1b." + scriptExt, "w") as out19:
        out19.write(start_log_redirect('../infiles/res1b.log') +
                    "./cluster." + scriptExt + "\n"
                    "./chimera." + scriptExt + "\n"
                    "./loci_sel." + scriptExt + "\n"
                    "./orient_merged_1b." + scriptExt + "\n"
                    + end_log_redirect())
    os.chdir('scripts')
    subprocess.run([shellCmd, "./runall1b." + scriptExt])


def runall_1c():
    """ Defines the run order for option 5
    """
    with open("scripts/runall1c." + scriptExt, "w") as out19:
        out19.write(start_log_redirect('../infiles/res1c.log') +
                    "./cluster_r1." + scriptExt + "\n"
                    "./chimera_r1." + scriptExt + "\n"
                    "./locir1_sel." + scriptExt + "\n"
                    "./orient_R1_5." + scriptExt + "\n"
                    + end_log_redirect())
    os.chdir('scripts')
    subprocess.run([shellCmd, "./runall1c." + scriptExt])


def runall_1d():
    """ Defines the run order for option 1d
    """
    with open("scripts/runall1d." + scriptExt, "w") as out19:
        out19.write(start_log_redirect('../infiles/res1d.log') +
                    "./cluster_one_sample_1d." + scriptExt + "\n"
                    "./chimera_one_sample_1d." + scriptExt + "\n"
                    "./loci_merged_1d." + scriptExt + "\n"
                    "./loci_R1_1d." + scriptExt + "\n"
                    "./orient_merged_1d." + scriptExt + "\n"
                    "./orient_R1_1d." + scriptExt + "\n"
                    + end_log_redirect())
    os.chdir('scripts')
    subprocess.run([shellCmd, "./runall1d." + scriptExt])


# ANALYSIS OF NUMBER OF READS CLUSTERS BY LOCI ETC..###################################################################
def nb_clus():
    """Calculates the number of resulting sequences (reads, merged, dereplicates and clusters) according
    to the selected options 'minsize', minseqlength, alpha parameter' and 'identity'.
    """
    os.chdir("../")
    with open(sample_user, "rt") as out32, open("infiles/counts.txt", "w") as out33:
        lines = out32.readlines()
        for line in lines:
            line = line.rstrip()
            # Nb READS Calculated on R1.fa
            r1fa = open("infiles/" + line + "_R1.fa", "rt")
            reads = r1fa.read()
            nb_reads = reads.count(">")
            # Nb merged .fa
            merged = open("infiles/" + line + ".fa", "rt")
            mgd = merged.read()
            nb_merged = mgd.count(">")
            percent_merging = (nb_merged/nb_reads)*100
            percent = '{:.2f}%'.format(percent_merging)
            # Number of dereplicated sequences calculated on _derep.fas
            derepfas = open("infiles/" + line + "_derep.fas", "rt")
            df = derepfas.read()
            nb_derpm = df.count("sample")
            # Nb clustermerged Calculated on _cluster.fas
            clusmerged = open("infiles/" + line + "_cluster.fas", "rt")
            clusmer = clusmerged.read()
            nb_clusm = clusmer.count("sample")
            # Number of merged clusters without chimeras calculated on _cluster_OK.fas
            clusmergedok = open("infiles/" + line + "_cluster_OK.fas", "rt")
            clusmerok = clusmergedok.read()
            nb_clusmok = clusmerok.count("sample")
            # Number of dereplicated sequences (R1) calculated on _derep_R1.fas
            derepfasr1 = open("infiles/" + line + "_derep_R1.fas", "rt")
            dfr1 = derepfasr1.read()
            nb_derpr1 = dfr1.count("sample")
            # Nb clusterR1 Calculated on _cluster_R1.fas
            clusfasr1 = open("infiles/" + line + "_cluster_R1.fas", "rt")
            clusr1 = clusfasr1.read()
            nb_clusr1 = clusr1.count("sample")
            # Nb clusterR1OK Calculated on _cluster_R1_OK.fas
            clusfasr1ok = open("infiles/" + line + "_cluster_R1_OK.fas", "rt")
            clusr1ok = clusfasr1ok.read()
            nb_clusr1ok = clusr1ok.count("sample")
            out33.writelines(f"The sample {line} has:\n"
                             f"\t{nb_reads} reads\n"
                             f"\t{nb_merged} merged sequences\n"
                             f"\tThe percentage of merging is {percent}\n"
                             f"\t{nb_derpm} dereplicated merged sequences\n"
                             f"\t{nb_clusm} merged clusters\n"
                             f"\t{nb_clusmok} merged clusters without chimera (OK)\n"
                             f"\t---------------------\n"
                             f"\t{nb_derpr1} dereplicated R1 sequences\n"
                             f"\t{nb_clusr1} R1 clusters\n"
                             f"\t{nb_clusr1ok} R1 clusters without chimera (R1_OK)\n\n\n")


def nb_seqbyloc_merged():
    """Calculates the number of relevant cluster sequences by locus according
    to the selected options 'minsize', minseqlength, alpha parameter' and 'identity' for merged
    sequences (amplicons based on paired-end reads).
    """
    with open("infiles/nbseq_bylocmerged.txt", "w") as out34:
        for sample in samples:
            out34.write(f"\nThe sample {sample} showed:\n")
            for locus1 in loci1s:
                os.chdir(locus1)
                refs = open(sample + "_merged.fas")
                refsloc = refs.read()
                nb_ref = refsloc.count("sample")
                out34.writelines(f"\t{nb_ref} clusters of merged sequences for the locus {locus1}\n")
                os.chdir("../")


def nb_seqbyloc_r1():
    """Calculates the number of relevant cluster sequences by locus according
    to the selected options 'minsize', minseqlength, alpha parameter' and 'identity' for R1
    sequences (amplicons with no mergeable R1/R2 reads).
    """
    with open("infiles/nbseq_bylocR1.txt", "w") as out34:
        for sample in samples:
            out34.write(f"\nThe sample {sample} showed:\n")
            for locus2 in loci2s:
                os.chdir(locus2)
                refs = open(sample + "_R1.fas")
                refsloc = refs.read()
                nb_ref = refsloc.count("sample")
                out34.writelines(f"\t{nb_ref} clusters of R1 sequences for the locus {locus2}\n")
                os.chdir("../")


# TRIM PRIMERS AND SELECT SEQUENCES ACCORDING THRESHOLDS ###############################################################
def trim_select_alls():
    """Select the final clusters of merged sequences by locus and by sample according to the size threshold retained
    for each of them.
    """
    dirloc = input("--------------------------------------------------------\n"
                   "Which LOCUS based on paired-end mergeable reads do you want to analyze? \n"
                   "Check that the corresponding folder does exists, e.g. 16S: ")
    while dirloc != "end":
        while Path(dirloc).is_dir() is False:
            dirloc = input("--------------------------------------------------------\n"
                           "Your locus name is not valid. Please enter a valid name for the locus: ")
        os.chdir(dirloc)
        trim_left = input(f"--------------------------------------------------------\n"
                          f"What is the number of bp of the left primer for {dirloc}? e.g. 20: ")
        trim_right = input(f"--------------------------------------------------------\n"
                           f"What is the number of bp of the right primer for {dirloc}? e.g. 22: ")
        ts = input(f"-----------------------------------------------\n"
                   f"What THRESHOLD you want for this locus {dirloc}\n"
                   f"e.g. for 5%, enter 0.05: ")
        for sample in samples:
            with open(sample + "_orient.fas", 'r') as filin, open("trim-select." + scriptExt, "w") as out20:
                targets = [line for line in filin if "size" in line]
                a = 0
                for target in targets:
                    size = re.search('size=(.+?)$', target).group(1)
                    a = a + int(size)
                b = int(a * float(ts) + 1)
                out20.writelines(f"" + start_log_redirect('./' + sample + '.log') +
                                 f"#Sum of sizes for {sample} = {a}'\n"
                                 f"#Threshold set to: {ts}'\n"
                                 f"#The sizes > {b} were conserved'\n"
                                 f"#Trimming left and right primers:'\n"
                                 f"#---------------------------------'\n"
                                 f"vsearch --fastx_filter {sample}_orient.fas --fastq_stripleft {trim_left} "
                                 f" --fastq_stripright {trim_right} --fastaout tmp --minsize {b}" + localErrorOnStopCmd + "\n"
                                 f"#Dereplication after trimming:'\n"
                                 f"#---------------------------------'\n"
                                 f"vsearch --derep_fulllength ./tmp --output {sample}_select.fas --sizein --sizeout" + localErrorOnStopCmd + "\n"
                                 f"" + end_log_redirect())
            sys.stdout.write("----------------------")
            sys.stdout.write(f"\nSum of sizes for {sample} = {a}")
            sys.stdout.write(f"\nThe sizes > {b} were conserved")
            sys.stdout.write("----------------------")
            subprocess.run([shellCmd, "./trim-select." + scriptExt])
            os.remove("tmp")
            os.remove("trim-select." + scriptExt)
        os.chdir('../')
        dirloc = input("--------------------------------------------------------\n"
                       "Which new LOCUS based on paired-end mergeable reads do you want to analyze? \n"
                       "Check that the corresponding folder does exists, e.g. nd4\n"
                       "Enter the NEW locus name or if you are done enter: 'end': ")


def trim_select_alls_r1():
    """Select the final clusters of merged sequences by locus and by sample according to the size threshold retained
    for each of them.
    """
    dirloc = input("--------------------------------------------------------\n"
                   "Which LOCUS based on single-end R1 reads do you want to analyze? \n"
                   "Check that the corresponding folder does exists, e.g. gpi: ")
    while dirloc != "end":
        while Path(dirloc).is_dir() is False:
            dirloc = input("--------------------------------------------------------\n"
                           "Your locus name is not valid. Please enter a valid name for the locus: ")
        os.chdir(dirloc)
        trim_left = input(f"--------------------------------------------------------\n"
                          f"What is the number of bp of the left primer for {dirloc}? e.g. 20: ")
        trim_right = input(f"--------------------------------------------------------\n"
                           f"What is the number of nucleotides to trim at the right of the "
                           f"fragment {dirloc}? e.g. 5: ")
        ts = input(f"-----------------------------------------------\n"
                   f"What THRESHOLD you want for this locus {dirloc}\n"
                   f"e.g. for 5%, enter 0.05: ")
        for sample in samples:
            with open(sample + "_R1_orient.fas", 'r') as filin, open("trim-select." + scriptExt, "w") as out20:
                targets = [line for line in filin if "size" in line]
                a = 0
                for target in targets:
                    size = re.search('size=(.+?)$', target).group(1)
                    a = a + int(size)
                b = int(a * float(ts) + 1)
                out20.writelines(f"" + start_log_redirect('./{sample}.log') +
                                 f"#Sum of sizes for {sample} = {a}'\n"
                                 f"#Threshold set to: {ts}'\n"
                                 f"#The sizes > {b} were conserved'\n"
                                 f"#Trimming left and right primers:'\n"
                                 f"#---------------------------------'\n"
                                 f"vsearch --fastx_filter {sample}_R1_orient.fas --fastq_stripleft {trim_left} "
                                 f" --fastq_stripright {trim_right} --fastaout tmp --minsize {b}" + localErrorOnStopCmd + "\n"
                                 f"#Dereplication after trimming:'\n"
                                 f"#---------------------------------'\n"
                                 f"vsearch --derep_fulllength ./tmp --output {sample}_R1_select.fas --sizein "
                                 f"--sizeout" + localErrorOnStopCmd + "\n"
                                 f"" + end_log_redirect())
            sys.stdout.write("----------------------")
            sys.stdout.write(f"\nSum of sizes for {sample}= {a}")
            sys.stdout.write(f"\nThe sizes > {b} were conserved")
            sys.stdout.write("----------------------")
            subprocess.run([shellCmd, "./trim-select." + scriptExt])
            os.remove("tmp")
            os.remove("trim-select." + scriptExt)
        os.chdir('../')
        dirloc = input("--------------------------------------------------------\n"
                       "Which new LOCUS based on single-end R1 reads do you want to analyze? \n"
                       "Check that the corresponding folder does exists, e.g. nd4\n"
                       "Enter the NEW locus name or if you are done enter: 'end': ")


def trim_select():
    """Select the final clusters of merged sequences by locus and by sample according to the size threshold retained
    for each of them.
    """
    dirloc = input("--------------------------------------------------------\n"
                   "Which LOCUS based on paired-end mergeable reads do you want to analyze? \n"
                   "Check that the corresponding folder does exists, e.g. 16S: ")
    while dirloc != "end":
        while Path(dirloc).is_dir() is False:
            dirloc = input("--------------------------------------------------------\n"
                           "Your locus name is not valid. Please enter a valid name for the locus: ")
        os.chdir(dirloc)
        trim_left = input(f"--------------------------------------------------------\n"
                          f"What is the number of bp of the left primer for {dirloc}? e.g. 20: ")
        trim_right = input(f"--------------------------------------------------------\n"
                           f"What is the number of bp of the right primer for {dirloc}? e.g. 22: ")
        orient = input(f"------------------------------------------------------\n"
                       f"Which SAMPLE x_orient.fas do you want to trim for the locus {dirloc}\n"
                       "only enter the sample name x, e.g. sample1: ")
        while orient != "end":
            while os.path.isfile(orient + "_orient.fas") is False:
                orient = input(f"--------------------------------------------------------\n"
                               f"Your sample name '{orient}' is not valid, please enter a valid name: ")
            ts = input(f"--------------------------------------------------------\n"
                       f"What THRESHOLD you want for this sample {orient} of locus {dirloc}\n"
                       f"e.g. for 5%, enter 0.05: ")
            with open(orient + "_orient.fas", 'r') as filin, open("trim-select." + scriptExt, "w") as out20:
                targets = [line for line in filin if "size" in line]
                a = 0
                for target in targets:
                    size = re.search('size=(.+?)$', target).group(1)
                    a = a + int(size)
                b = int(a * float(ts) + 1)
                out20.writelines(f'' + start_log_redirect('./' + orient + '.log') +
                                 f'#Sum of sizes for {orient} = {a}\n'
                                 f'#Threshold set to: {ts}\n'
                                 f'#The sizes > {b} were conserved\n'
                                 f'#Trimming left and right primers:\n'
                                 f'#---------------------------------\n'
                                 f'vsearch --fastx_filter {orient}_orient.fas --fastq_stripleft {trim_left} '
                                 f' --fastq_stripright {trim_right} --fastaout tmp --minsize {b}\n'
                                 f'#Dereplication after trimming:\n'
                                 f'#---------------------------------\n'
                                 f'vsearch --derep_fulllength ./tmp --output {orient}_select.fas --sizein --sizeout\n'
                                 f'' + end_log_redirect())
            sys.stdout.write("----------------------")
            sys.stdout.write(f"\nSum of sizes for {orient} = {a}")
            sys.stdout.write(f"\nThe sizes > {b} were conserved")
            sys.stdout.write("----------------------")
            subprocess.run([shellCmd, "./trim-select." + scriptExt])
            os.remove("tmp")
            os.remove("trim-select." + scriptExt)
            orient = input(f"--------------------------------------------------------\n"
                           f"Which NEW FILE do you want to trim for the locus {dirloc}?\n"
                           f"enter the sample name or, if you terminated with the locus {dirloc} enter 'end': ")
        os.chdir('../')
        dirloc = input("--------------------------------------------------------\n"
                       "Which new LOCUS based on paired-end mergeable reads do you want to analyze? \n"
                       "Check that the corresponding folder does exists, e.g. nd4\n"
                       "Enter the NEW locus name or if you are done enter: 'end': ")


def trim_select_r1():
    """Select the final clusters of merged sequences by locus and by sample according to the size threshold retained
    for each of them.
    """
    dirloc = input("--------------------------------------------------------\n"
                   "Which LOCUS based on single-end R1 reads do you want to analyze? \n"
                   "Check that the corresponding folder does exists, e.g. gpi: ")
    while dirloc != "end":
        while Path(dirloc).is_dir() is False:
            dirloc = input("--------------------------------------------------------\n"
                           "Your locus name is not valid. Please enter a valid name for the locus: ")
        os.chdir(dirloc)
        trim_left = input(f"--------------------------------------------------------\n"
                          f"What is the number of bp of the left primer for {dirloc}? e.g. 20: ")
        trim_right = input(f"--------------------------------------------------------\n"
                           f"What is the number of nucleotides to trim at the right of the fragment"
                           f" {dirloc}? e.g. 5: ")
        orient = input(f"------------------------------------------------------\n"
                       f"Which SAMPLE x_R1_orient.fas do you want to trim for the locus {dirloc}\n"
                       "only enter the sample name x, e.g. sample1: ")
        while orient != "end":
            while os.path.isfile(orient + "_R1_orient.fas") is False:
                orient = input(f"--------------------------------------------------------\n"
                               f"Your sample name '{orient}' is not valid, please enter a valid name: ")
            ts = input(f"--------------------------------------------------------\n"
                       f"What THRESHOLD you want for this sample {orient} of locus {dirloc}\n"
                       f"e.g. for 5%, enter 0.05: ")
            with open(orient + "_R1_orient.fas", 'r') as filin, open("trim-select." + scriptExt, "w") as out20:
                targets = [line for line in filin if "size" in line]
                a = 0
                for target in targets:
                    size = re.search('size=(.+?)$', target).group(1)
                    a = a + int(size)
                b = int(a * float(ts) + 1)
                out20.writelines(f'' + start_log_redirect('./' + orient + '.log') +
                                 f'#Sum of sizes for {orient} = {a}\n'
                                 f'#Threshold set to: {ts}\n'
                                 f'#The sizes > {b} were conserved\n'
                                 f'#Trimming left and right primers:\n'
                                 f'#---------------------------------\n'
                                 f'vsearch --fastx_filter {orient}_R1_orient.fas --fastq_stripleft {trim_left} '
                                 f' --fastq_stripright {trim_right} --fastaout tmp --minsize {b}\n'
                                 f'#Dereplication after trimming:\n'
                                 f'#---------------------------------\n'
                                 f'vsearch --derep_fulllength ./tmp --output {orient}_R1_select.fas --sizein '
                                 f'--sizeout\n'
                                 f'' + end_log_redirect())
            sys.stdout.write("----------------------")
            sys.stdout.write(f"\nSum of sizes for {orient} = {a}")
            sys.stdout.write(f"\nThe sizes > {b} were conserved")
            sys.stdout.write("----------------------")
            subprocess.run([shellCmd, "./trim-select." + scriptExt])
            os.remove("tmp")
            os.remove("trim-select." + scriptExt)
            orient = input(f"--------------------------------------------------------\n"
                           f"Which NEW FILE do you want to trim for the locus {dirloc}?\n"
                           f"enter the sample name or, if you are done with the locus {dirloc} enter 'end': ")
        os.chdir('../')
        dirloc = input("--------------------------------------------------------\n"
                       "Which new LOCUS based on single-end R1 reads do you want to analyze? \n"
                       "Check that the corresponding folder does exists, e.g. nd4\n"
                       "Enter the NEW locus name or if you are done enter: 'end': ")


# CLUSTERING OF ALL SAMPLES BY LOCUS FOR PHYLOGENETIC ANALYZES ########################################################
def concat():
    """CLUSTERING OF ALL SAMPLES BY LOCUS FOR PHYLOGENETIC ANALYZES
    """
    loc2cat = input("----------------------------------------------------------\n"
                    "\nFor which LOCUS do you want to cluster all the sample sequences?: ")
    while loc2cat != "end":
        while Path(loc2cat).is_dir() is False:
            loc2cat = input("\nYour locus name is not valid. Please enter a valid name for the locus: ")
        os.chdir(loc2cat)
        files2cat = glob.glob('*_select.fas')
        with open("./" + loc2cat + "_allseq_select.fasta", "w") as out:
            for file in files2cat:
                with open(file, "r") as out2:
                    out.write(out2.read())
        sys.stdout.write(f"\nThe files {files2cat} have been clustered for locus {loc2cat}\n")
        os.chdir('../')
        loc2cat = input("----------------------------------------------------\n"
                        "\nFor which NEW LOCUS do you want to cluster all the sample sequences?\n"
                        "Enter the NEW locus name or 'end' if you terminated: ")
    else:
        sys.stdout.write("\n\n**** YOUR CLUSTERING SESSION FOR PHYLOGENETIC ANALYZES is COMPLETE ****\n\n")
        exit()


#######################################################################################################################
if __name__ == "__main__":
    date = datetime.datetime.now()
    current_dir = os.getcwd()
    print("\n#################################################################################################\n"
          "STEPS 1, 2 and 3 of the MAIN MENU are mandatory for a complete analysis while others are optional\n"
          "#################################################################################################")
    rmenu = input("\n************************************** MAIN MENU *************************************\n\n"
                  "0  -> Initial and optional quality checking of your fastq files (slow - be patient!)\n"
                  "1  -> NEW COMPLETE ANALYSIS\n"
                  "1a -> Re-analyze all loci, from the clustering step, modifying parameters\n"
                  "1b -> Re-analyze only one locus of paired-end amplicon (merged reads), modifying parameters\n"
                  "1c -> Re-analyze only one locus of single-end amplicon (R1 only), modifying parameters\n"
                  "1d -> Re-analyse only one sample, modifying parameters\n"
                  "2  -> SELECTION OF MINIMUM SIZES ACCORDING TO THRESHOLDS SET BY THE USER\n"
                  "3  -> CLUSTERING OF ALL SAMPLES BY LOCUS FOR PHYLOGENETIC ANALYZES\n"
                  "end -> Exit the program\n\n"
                  " ********* Type 0, 1, 1a, 1b, 1c, 1d, 2, 3 or end: ")

    # Optional statistical tests on fastq files
    if rmenu == "0":
        input_dir_fastq()
        input_fastqr1_user()
        input_fastqr2_user()
        input_sample_user()
        folder_infiles()
        folder_scripts()
        stat()

        sys.stdout.write("\n\n")
        if runall_0() > 0:
            print("\nQuality checking execution failed, please check infiles/results.log")
            exit(1)

        sys.stdout.write("\n**** PROCEDURE 0 (quality checking) EXECUTION IS COMPLETE ****\n\n")

    # NEW ANALYSIS
    if rmenu == "1":
        input_dir_fastq()
        input_fastqr1_user()
        input_fastqr2_user()
        input_loci1_user()
        input_loci2_user()
        input_sample_user()
        input_minsize_user()
        input_minseqlength()
        input_alpha()
        input_identity()
        param_1()
        folder_infiles()
        folder_scripts()
        folder_loci()
        folder_refs()
        merging()
        fastq2fas()
        derep_merged()
        derep_r1()
        cluster_merged()
        cluster_r1()
        chimera_merged()
        chimera_r1()
        runloc_merged()
        runloc_r1()
        orient_merged()
        orient_r1()

        sys.stdout.write("\n\n")
        if runall_1() > 0:
            print("\nMain analysis execution failed, please check infiles/results.log")
            exit(1)

        nb_clus()
        nb_seqbyloc_merged()
        nb_seqbyloc_r1()
        sys.stdout.write("\n**** PROCEDURE 1 (main analysis) EXECUTION IS COMPLETE ****\n\n")

    # Re-analyze all loci, from the clustering step, modifying parameters
    if rmenu == "1a":
        prev_param()
        input_minsize_user()
        input_minseqlength()
        input_alpha()
        input_identity()
        param_1a()
        cluster_merged()
        cluster_r1()
        chimera_merged()
        chimera_r1()
        runloc_merged()
        runloc_r1()
        orient_merged()
        orient_r1()
        runall_1a()
        nb_clus()
        nb_seqbyloc_merged()
        nb_seqbyloc_r1()
        sys.stdout.write("\n\n**** YOUR RUN OPTION 1a IS COMPLETE ****\n\n")

    # Re-analyze only one locus of paired-end amplicon (merged reads), modifying parameters
    if rmenu == "1b":
        prev_param()
        input_loc_sel_merged()
        input_minsize_user()
        input_minseqlength()
        input_alpha()
        input_identity()
        param_1b()
        cluster_merged()
        chimera_merged()
        runlocsel_merged()
        orient_one_loc_1b()
        runall_1b()
        sys.stdout.write("\n\n**** YOUR RUN OPTION 1b IS COMPLETE ****\n\n")

    # Re-analyze only one locus of single-end amplicon (R1 only), modifying parameters
    if rmenu == "1c":
        prev_param()
        input_loc_sel_r1()
        input_minsize_user()
        input_minseqlength()
        input_alpha()
        input_identity()
        param_1c()
        cluster_r1()
        chimera_r1()
        runlocsel_r1()
        orient_one_loc_1c()
        runall_1c()
        sys.stdout.write("\n\n**** YOUR RUN OPTION 1c IS COMPLETE ****\n\n")

    # Re-analyse only one sample, modifying parameters
    if rmenu == "1d":
        prev_param()
        input_sam_sel()
        input_minsize_user()
        input_minseqlength()
        input_alpha()
        input_identity()
        param_1d()
        cluster_one_sample_1d()
        chimera_one_sample_1d()
        runloc_one_sample_1d()
        orient_one_sample_1d()
        runall_1d()
        sys.stdout.write("\n\n**** YOUR RUN OPTION 1d IS COMPLETE ****\n\n")

    # SELECTION OF MINIMUM SIZES ACCORDING TO THRESHOLDS SET BY THE USER
    if rmenu == '2':
        submenu = input("\n***************** MINIMUM SIZE SELECTION MENU *****************\n\n"
                        "1- Apply the SAME size threshold for ALL SAMPLES for the loci based on PAIRED-END reads "
                        "(R1/R2 merged)\n"
                        "i.e. you want to keep only sequences whose abundance (size)\n"
                        "is greater than x% of the total number of sequences for a given sample.\n"
                        "This threshold of x% can be chosen for each locus.\n\n"
                        "2- Apply the SAME size threshold for ALL SAMPLES for the loci based on SINGLE-END reads "
                        "(R1 only)\n"
                        "idem than option 1 but only using the R1 reads instead of merged ones.\n\n"
                        "3- Apply a SPECIFIC size threshold for EACH SAMPLE, for the loci based on PAIRED-END reads "
                        "(R1/R2 merged)\n"
                        "i.e. you want to modulate the threshold of x% by locus but also by sample\n"
                        "within a particular locus.\n\n"
                        "4- Apply a SPECIFIC size threshold for EACH SAMPLE, for the loci based on SINGLE-END reads "
                        "(R1 only)\n"
                        "idem than option 3 but only using the R1 sequences instead of merged ones.\n\n"
                        " ************************************** Type 1, 2, 3, 4 or end: ")
        if submenu == "1":
            prev_param()
            trim_select_alls()
            sys.stdout.write("\n\n**** YOUR SIZE SELECTION for loci based on PAIRED-END reads is COMPLETE ****\n\n")

        if submenu == "2":
            prev_param()
            trim_select_alls_r1()
            sys.stdout.write("\n\n**** YOUR SIZE SELECTION for loci base on SINGLE-END reads (R1 only) "
                             "is COMPLETE ****\n\n")

        if submenu == "3":
            prev_param()
            trim_select()
            sys.stdout.write("\n\n**** YOUR SIZE SELECTION by sample and by loci based on PAIRED-END is "
                             "COMPLETE ****\n\n")

        if submenu == "4":
            prev_param()
            trim_select_r1()
            sys.stdout.write("\n\n**** YOUR SIZE SELECTION by sample and by loci based on SINGLE-END is "
                             "COMPLETE ****\n\n")

        if submenu == "end":
            exit()

    # CLUSTERING OF ALL SAMPLES BY LOCUS FOR PHYLOGENETIC ANALYZES
    if rmenu == "3":
        concat()

    # Exit the program
    if rmenu == "end":
        sys.stdout.write("\n\n**** THANK YOU FOR USING MBCTOOLS ****\n\n")
        exit()
