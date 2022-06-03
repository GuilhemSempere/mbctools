#!/usr/bin/python3

"""This program is a toolkit to make the use of VSEARCH easier and interactive, to analyze your
Metabarcoding NGS data in best conditions. It proposes the following MAIN MENU:

1- NEW complete analyze (slow - be patient!)
2- NEW analyze WITHOUT the unnecessary statistical step (faster)
3- Re-analyze all loci, from the clustering step and with other parameters (faster)
4- Re-analyze only one locus of amplicon < 550bp, with other parameters
5- Re-analyze only one locus of amplicon > 550bp, with other parameters
6- Re-analyse only one sample, with other parameters
7- Clean unnecessary files '*.fa' files to free up space
8- Blast at NCBI a cluster of sequences for a sample
9- Retrieve a particular sequence from GenBank knowing its gb code (accession number)
10- Select a size threshold for each sample and locus, trim the sequences of loci < 550bp (merged)
11- Select a size threshold for each sample and locus, trim the sequences of loci > 550bp (R1)
12- Concatenate all the files '*_select.fas' for a final analysing in MEGA7
13- Exit the program

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

global dir_fastq, fastqr1_user, fastqr2_user, loci1_user, loci2_user, sample_user, minsize_user, minseqlength, alpha
global identity, loc_sel, fastqr1s, fastqr2s, sam_sel, loci1s, loci2s, samples

global winOS, shellCmd, scriptExt, fileSep
winOS = "Windows" in platform.uname()
shellCmd = "powershell" if winOS else "sh"
scriptExt = "ps1" if winOS else "sh"
fileSep = "\\" if winOS else "/"


# ------------------------- FUNCTIONS FOR MULTIPLATFORM SUPPORT -----------------------------------

def startLogRedirect(filePath):
    if not winOS:
        return "exec 3>&1 4>&2 >" + filePath + " 2>&1\n"
    else:
        return "&{\n"


def endLogRedirect():
    if not winOS:
        return "exec 1>&3 2>&4\n"
    else:
        return "} 2> '../infiles/res2.log'\n"


def mainStreamMessage(message):
    if not winOS:
        return "printf \"" + message + "\" >&3\n"
    else:
        return "Write-Host -NoNewline \"" + message + "\"\n"


# ------------------------------ FOLDER CREATION -----------------------------
def folder_1_2():
    """ Creates relevant folders for your Metabarcoding analysis: one folder for each locus, one for the Vsearch
    scripts and another for intermediate files named 'infiles'
    """
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

    folder = "infiles"
    path = os.path.join(current_dir, folder)
    if Path(path).is_dir():
        sys.stdout.write(f"\nThe folder {folder} already exists and will be used in the current analysis")
    else:
        os.mkdir(path)

    folder = "scripts"
    path = os.path.join(current_dir, folder)
    if Path(path).is_dir():
        sys.stdout.write(f"\nThe folder {folder} already exists and will be used in the current analysis")
    else:
        os.mkdir(path)

    folder = "refs"
    path = os.path.join(current_dir, folder)
    if Path(path).is_dir():
        sys.stdout.write(f"\nThe folder {folder} already exists and will be used in the current analysis")
    else:
        os.mkdir(path)

    sys.stdout.write("\n\n")


# ----------------------------- INPUTS AND VARIABLES -------------------------
def inputs_1_2():
    """ Defines global variables for options 1 and 2: new analyses with or without statistics on fastq files
    """
    global dir_fastq
    dir_fastq = input("\nEnter the entire path of the folder where are the fastq files "
                      "\ne.g. C:/Users/barnabe/Documents/NGSdata/run190823"
                      "\n default = current_directory/fastq : ")
    if dir_fastq == '':
        dir_fastq = current_dir + "/fastq"
    else:
        while Path(dir_fastq).is_dir() is False:
            dir_fastq = input("\nYour path is not valid. Please enter a valid name for the path : ")

    global fastqr1_user, fastqr1s
    fastqr1_user = input("\nEnter the name of your fastqR1 file, default = fastqR1.txt : ")
    if fastqr1_user == '':
        fastqr1_user = 'fastqR1.txt'
    else:
        while os.path.isfile(fastqr1_user) is False:
            fastqr1_user = input("\nYour name is not valid, please enter a valid name : ")
    with open(fastqr1_user, "r") as out1:
        fastqr1s = out1.read().splitlines()

    global fastqr2_user, fastqr2s
    fastqr2_user = input("\nEnter the name of your fastqR2 file, default = fastqR2.txt : ")
    if fastqr2_user == '':
        fastqr2_user = 'fastqR2.txt'
    else:
        while os.path.isfile(fastqr2_user) is False:
            fastqr2_user = input("\nYour name is not valid, please enter a valid name : ")
    with open(fastqr2_user, "r") as out2:
        fastqr2s = out2.read().splitlines()

    global loci1_user, loci1s
    loci1_user = input("\nEnter the name of your loci file with amplicons < 550bp, default = locus1.txt "
                       "(one locus by line) : ")
    if loci1_user == '':
        loci1_user = 'locus1.txt'
    else:
        while os.path.isfile(loci1_user) is False:
            loci1_user = input("\nYour name is not valid, please enter a valid name : ")
    with open(loci1_user, "r") as out3:
        loci1s = out3.read().splitlines()

    global loci2_user, loci2s
    loci2_user = input("\nEnter the name of your loci file with amplicons > 550bp, default = locus2.txt"
                       " (one locus by line) : ")
    if loci2_user == '':
        loci2_user = 'locus2.txt'
    else:
        while os.path.isfile(loci2_user) is False:
            loci2_user = input("\nYour name is not valid, please enter a valid name : ")
    with open(loci2_user, "r") as out4:
        loci2s = out4.read().splitlines()

    global sample_user
    sample_user = input("\nEnter the name of your sample file, default = samples.txt : ")
    if sample_user == '':
        sample_user = 'samples.txt'
    else:
        while os.path.isfile(sample_user) is False:
            sample_user = input("\nYour name for sample file is not valid, please enter a valid name : ")

    global samples
    with open(sample_user, "r") as out5:
        samples = out5.read().splitlines()

    global minsize_user
    minsize_user = input("\nEnter the minsize option value for your clusters, default = 8 : ")
    if minsize_user == '':
        minsize_user = '8'

    global minseqlength
    minseqlength = input("\nEnter your minseqlength parameter for clustering, default = 100 : ")
    if minseqlength == '':
        minseqlength = '100'

    global alpha
    alpha = input("\nEnter your alpha parameter for clustering, default = 2 : ")
    if alpha == '':
        alpha = '2'

    global identity
    identity = input("\nEnter your id parameter for usearch global, default = 0.7 : ")
    if identity == '':
        identity = '0.7'


def inputs_3():
    """ Defines global variables necessary for option 3: re-analyse all loci from the clustering step
    """
    global minsize_user
    minsize_user = input("\nEnter the minsize option value for your clusters, default = 8 : ")
    if minsize_user == '':
        minsize_user = '8'

    global minseqlength
    minseqlength = input("\nEnter your minseqlength parameter for clustering, default = 100 : ")
    if minseqlength == '':
        minseqlength = '100'

    global alpha
    alpha = input("\nEnter your alpha parameter for clustering, default = 2 : ")
    if alpha == '':
        alpha = '2'

    global identity
    identity = input("\nEnter your id parameter for usearch global, default = 0.7 : ")
    if identity == '':
        identity = '0.7'


def inputs_4_5():
    """ Defines global variables necessary for option 4 or 5: re-analysis of only one locus < 550 bp or > 550 bp
    """
    global loc_sel
    loc_sel = input("\nEnter the name of the locus you want rerun, no default : ")
    while Path(loc_sel).is_dir() is False:
        loc_sel = input("\nYour locus name is not valid. Please enter the valid name of the locus : ")

    global minsize_user
    minsize_user = input("\nEnter the minsize option value for your clusters, default = 8 : ")
    if minsize_user == '':
        minsize_user = '8'

    global minseqlength
    minseqlength = input("\nEnter your minseqlength parameter for clustering, default = 100 : ")
    if minseqlength == '':
        minseqlength = '100'

    global alpha
    alpha = input("\nEnter your alpha parameter for clustering, default = 2 : ")
    if alpha == '':
        alpha = '2'
    global identity
    identity = input("\nEnter your id parameter for usearch global, default = 0.7 : ")
    if identity == '':
        identity = '0.7'


def inputs_6():
    """ Defines global variables for option 6: re-analysis of a single sample
    """
    global sam_sel
    sam_sel = input("\nEnter the sample you want to rerun, no default :  ")

    global minsize_user
    minsize_user = input("\nEnter the minsize option value for your clusters, default = 8 : ")
    if minsize_user == '':
        minsize_user = '8'

    global minseqlength
    minseqlength = input("\nEnter your minseqlength parameter for clustering, default = 100 : ")
    if minseqlength == '':
        minseqlength = '100'

    global alpha
    alpha = input("\nEnter your alpha parameter for clustering, default = 2  : ")
    if alpha == '':
        alpha = '2'

    global identity
    identity = input("\nEnter your id parameter for usearch global (0.7 - 0.9  may be convenient) : ")
    if identity == '':
        identity = '0.7'


# -------------------------------- PARAMETER FILES ------------------------
def param_1_2():
    """ Creates a file with one parameter by line for options 1 and 2
    """
    with open("my_parameters.txt", "w") as outfile:  # CREATION OF MY PARAMETERS
        outfile.writelines(f"Menu option 1: {date}\n{dir_fastq}\n{fastqr1_user}\n{fastqr2_user}\n{loci1_user}\n"
                           f"{loci2_user}\n{sample_user}\n{minsize_user}\n{minseqlength}\n{alpha}\n{identity}")


def prev_param():
    """ Recalls global variables for different options
    """
    with open("my_parameters.txt", "r") as infile:
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
        sample_user = lines[6]

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


def param_3():
    """ Creates a file with one parameter by line for option 3
    """
    with open("my_parameters_option3.txt", "w") as outfile:
        outfile.writelines(f"Menu option 3: {date}\n{dir_fastq}\n{fastqr1_user}\n{fastqr2_user}\n"
                           f"{loci1_user}\n{loci2_user}\n{sample_user}\n{minsize_user}\n{minseqlength}\n"
                           f"{identity}\n{alpha}")


def param_4():
    """Creates a file with one parameter by line for option 4
    """
    with open("my_parameters_option4.txt", "w") as outfile:
        outfile.writelines(f"Menu option 4: {date}\n{dir_fastq}\n{fastqr1_user}\n{fastqr2_user}\n{loc_sel}\n"
                           f"{sample_user}\n{minsize_user}\n{minseqlength}\n{alpha}\n{identity}")


def param_5():
    """Creates a file with one parameter by line for option 5
    """
    with open("my_parameters_option5.txt", "w") as outfile:
        outfile.writelines(f"Menu option 5: {date}\n{dir_fastq}\n{fastqr1_user}\n{fastqr2_user}\n{loc_sel}\n"
                           f"{sample_user}\n{minsize_user}\n{minseqlength}\n{alpha}\n{identity}")


def param_6():
    """Creates a file with one parameter by line for option 6
    """
    with open("my_parameters_option6.txt", "w") as outfile:
        outfile.writelines(f"Menu option 6: {date}\n{dir_fastq}\n{fastqr1_user}\n{fastqr2_user}\n{sam_sel}\n"
                           f"{minsize_user}\n{minseqlength}\n{alpha}\n{identity}")


# ------------------------------------- STAT ON FASTQ FILES ---------------------
def stat():
    """Tests the quality of each 'fastq' file by the VSEARCH command:
    vsearch --fastq_eestats2 ../infiles/fastqF-R1 --output ../infiles/SN_R1info.txt
    vsearch --fastq_eestats2 ../infiles/fastqF-R2 --output ../infiles/SN_R2info.txt
    fastqF = fastq file name

    SN = sample name
    """
    with open("scripts/infor1." + scriptExt, "w") as out3:
        i = 0
        out3.write(mainStreamMessage(f'Checking R1 stats:'))
        while i < len(samples):
            sample = samples[i]
            fastqr1 = fastqr1s[i]
            i = i + 1
            out3.write(mainStreamMessage(f' {sample}'))
            out3.write(f"vsearch --fastq_eestats2 {dir_fastq}{fileSep}{fastqr1} --output ../infiles/{sample}_R1info.txt\n")
        out3.write(mainStreamMessage(f'\n\n'))

    with open("scripts/infor2." + scriptExt, "w") as out4:
        i = 0
        out4.write(mainStreamMessage(f'Checking R2 stats:'))
        while i < len(samples):
            sample = samples[i]
            fastqr2 = fastqr2s[i]
            i = i + 1
            out4.write(mainStreamMessage(f' {sample}'))
            out4.write(f"vsearch --fastq_eestats2 {dir_fastq}{fileSep}{fastqr2} --output ../infiles/{sample}_R2info.txt\n")
        out4.write(mainStreamMessage(f'\n\n'))


# ------------------------------------- MERGING FOR AMPLICONS <  550 BP --------
def merging():
    """Merges paired-end reads into one sequence, when the length of the amplicon allows it (< 550bp)
    according the VSEARCH command:
    vsearch --fastq_mergepairs ../infiles/fastqR1 --reverse ../infiles/fastqR2 --fastaout
    ./infiles/SN.fa --fastq_allowmergestagger --relabel sample='sample=SN'_merged

    fastqR1 = fastqR1 complete name
    fastqR2 = fastqR2 complete name
    option --fastq_allowmergestagger allows the merging of short fragments
    """
    with open("scripts/merging." + scriptExt, "w") as out5:
        i = 0
        out5.write(mainStreamMessage(f'merging sample reads:'))
        while i < len(samples):
            sample = samples[i]
            fastqr1 = fastqr1s[i]
            fastqr2 = fastqr2s[i]
            out5.write(mainStreamMessage(f' {sample}...') +
                       f"vsearch --fastq_mergepairs {dir_fastq}{fileSep}{fastqr1} --reverse {dir_fastq}{fileSep}{fastqr2} --fastaout "
                       f"../infiles/{sample}.fa --fastq_allowmergestagger --fastq_maxee 1 --relabel "
                       f"sample={sample}_merged\n")
            i = i + 1
        out5.write(mainStreamMessage(f'\n\n'))


# ----------------------------------- FASTQR TO FAS FOR AMPLICONS >  550 BP -----
def fastq2fas():
    """When amplicons >550 bp, merging is impossible, the reads R1 of 301 bp (better than R2)
    are used to search the relevant sequences. First, all R1 'fastq' files have to be transformed into 'fasta' files
    by the VSEARCH command:
    vsearch --fastq_filter ../infiles/fastaqR1 --fastaout ../infiles/SN_R1.fa

    SN = sample name
    fastqR1 = fastqR1 file name
    """
    with open("scripts/fqtofas." + scriptExt, "w") as out6:
        i = 0
        out6.write(mainStreamMessage(f'fastqtofas:'))
        while i < len(samples):
            sample = samples[i]
            fastqr1 = fastqr1s[i]
            out6.write(mainStreamMessage(f' {sample}...') +
                       f"vsearch --fastq_filter {dir_fastq}{fileSep}{fastqr1} --fastaout ../infiles/{sample}_R1.fa\n")
            i = i + 1
        out6.write(mainStreamMessage(f'\n\n'))


# ----------------------------------------- DEREPLICATION -----------------------
def derep_merged():
    """Dereplicates merged sequences in a given 'fasta' file with the VSEARCH command:
    vsearch --derep_fulllength ../infiles/SN.fa --output ../infiles/SN_derep.fas --sizeout
    --strand both

    SN = sample name
    Dereplicates in both strands and writes abundance annotation (frequency) to output.
    """
    with open("scripts/derep." + scriptExt, "w") as out7:
        out7.write(mainStreamMessage(f'dereplication:'))
        for sample in samples:
            out7.write(mainStreamMessage(f' {sample}...') +
                       f"vsearch --derep_fulllength ../infiles/{sample}.fa --output ../infiles/{sample}_derep.fas "
                       f"--sizeout --strand both\n")
        out7.write(mainStreamMessage(f'\n\n'))


def derep_r1():
    """Dereplicates R1 sequences in a given 'fasta' file with the VSEARCH command:
    vsearch --derep_fulllength ../infiles/_R1_SN'.fa --output ../infiles/SN_derep.fas --sizeout
    --strand both --relabel sample='SN'_merged.

    SN = sample name
    """
    with open("scripts/derep_r1." + scriptExt, "w") as out8:
        out8.write(mainStreamMessage(f'dereplicationR1:'))
        for sample in samples:
            out8.write(mainStreamMessage(f' {sample}...') +
                       f"vsearch --derep_fulllength ../infiles/{sample}_R1.fa --output ../infiles/{sample}"
                       f"_derep_R1.fas --sizeout --strand both --relabel sample={sample}_R1.\n")
        out8.write(mainStreamMessage(f'\n\n'))


# --------------------------------------- UNOISE CLUSTERING ---------------------
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
        out9.write(mainStreamMessage(f'cluster:'))
        while i < len(samples):
            sample = samples[i]
            out9.write(mainStreamMessage(f' {sample}...') +
                       f"vsearch --cluster_unoise ../infiles/{sample}_derep.fas --sizein --centroids "
                       f"../infiles/{sample}_cluster.fas --strand both --minsize {minsize_user} --sizeout --sizeorder "
                       f"--unoise_alph {alpha} --minseqlength {minseqlength}\n")
            i = i + 1
        out9.write(mainStreamMessage(f'\n\n'))


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
        out10.write(mainStreamMessage(f'clusteringR1:'))
        while i < len(samples):
            sample = samples[i]
            out10.write(mainStreamMessage(f' {sample}...') +
                        f"vsearch --cluster_unoise ../infiles/{sample}_derep_R1.fas --sizein --centroids "
                        f"../infiles/{sample}_cluster_R1.fas --strand both --minsize {minsize_user} --sizeout "
                        f"--sizeorder --unoise_alpha {alpha} --minseqlength {minseqlength}\n")
            i = i + 1
        out10.write(mainStreamMessage(f'\n\n'))


def cluster_one_sample_6():
    """Denoises and clusters Illumina dereplicated for a selected sample and gives in output the centroids sequences
    to 'fasta' files, for a selected sample, with the VSEARCH command:
    vsearch --cluster_unoise ../infiles/SS_derep.fas --sizein --centroids ../infiles/SS_cluster.fas
    --strand both --minsize int --sizeout --sizeorder --unoise_alph int --minseqlength int

    SS = selected sample
    """
    with open("scripts/cluster_one_sample_6." + scriptExt, "w") as out1:
        out1.write(mainStreamMessage(f'clustering {sam_sel}') +
                   f"vsearch --cluster_unoise ../infiles/{sam_sel}_derep.fas --sizein --centroids "
                   f"../infiles/{sam_sel}_cluster.fas --strand both --minsize {minsize_user} --sizeout "
                   f"--sizeorder --unoise_alph {alpha} --minseqlength {minseqlength}\n"
                   f"vsearch --cluster_unoise ../infiles/{sam_sel}_derep_R1.fas --sizein --centroids "
                   f"../infiles/{sam_sel}_cluster_R1.fas --strand both --minsize {minsize_user} "
                   f"--sizeout --sizeorder --unoise_alph {alpha} --minseqlength {minseqlength}\n")
        out1.write(mainStreamMessage(f'\n\n'))


# ----------------------------------------- REMOVING CHIMERA ---------------------------
def chimera_merged():
    """Detects potential chimeras in denoised merged sequences by the VSEARCH command:
    vsearch --uchime3_denovo ../infiles/SN_cluster.fas --nonchimeras ../infiles/SN_cluster_OK.fas

    SN = sample name
    After demoising, the chimeras are very scarce.
    """
    with open("scripts/chimera." + scriptExt, "w") as out11:
        out11.write(mainStreamMessage(f'chimera:'))
        for sample in samples:
            out11.write(mainStreamMessage(f' {sample}...') +
                        f"vsearch --uchime3_denovo ../infiles/{sample}_cluster.fas --nonchimeras "
                        f"../infiles/{sample}_cluster_OK.fas\n")
        out11.write(mainStreamMessage(f'\n\n'))


def chimera_r1():
    """Detects potential chimeras in R1 sequences by the VSEARCH command:
    vsearch --uchime3_denovo ../infiles/SN_cluster-R1.fas --nonchimeras ../infiles/SN_cluster_R1_OK.fas

    SN = sample name
    After demoising, the chimeras are very scarce.
    """
    with open("scripts/chimera_r1." + scriptExt, "w") as out12:
        out12.write(mainStreamMessage(f'chimeraR1:'))
        for sample in samples:
            out12.write(mainStreamMessage(f' {sample}...') +
                        f"vsearch --uchime3_denovo ../infiles/{sample}_cluster_R1.fas --nonchimeras"
                        f" ../infiles/{sample}_cluster_R1_OK.fas\n")
        out12.write(mainStreamMessage(f'\n\n'))


def chimera_one_sample_6():
    """Detects potential chimeras in sequences for a selected sample by the VSEARCH command:
    vsearch --uchime3_denovo ../infiles/SS_cluster.fas --nonchimeras ../infiles/SS_cluster_OK.fas
    vsearch --uchime3_denovo ../infiles/SS_cluster_R1.fas --nonchimeras ../infiles/SS_cluster_R1_OK.fas

    ss = selected sample
    After demoising, the chimeras are very scarce
    """
    with open("scripts/chimera_one_sample_6." + scriptExt, "w") as out1:
        out1.write(mainStreamMessage(f'chimera {sam_sel}') +
                   f"vsearch --uchime3_denovo ../infiles/{sam_sel}_cluster.fas --nonchimeras ../infiles/{sam_sel}"
                   f"_cluster_OK.fas\n"
                   f"vsearch --uchime3_denovo ../infiles/{sam_sel}_cluster_R1.fas --nonchimeras "
                   f"../infiles/{sam_sel}_cluster_R1_OK.fas\n")
        out1.write(mainStreamMessage(f'\n\n'))


# ------------------- DISTRIBUTION OF CLUSTERS INTO THE DIFFERENT LOCI -----------
def runloc_merged():
    """Searches similarities between merged, denoised and non-chimera sequences and the local reference
    database (-db) by the VSEARCH command:
    vsearch --usearch_global ../infiles/SN_cluster_OK.fas --db ../refs/L1.fas --matched ../L1/SN_merged.fas
    --id int --strand both

    L1 = locus name for amplicons < 550 bp
    This is the second crucial step metabarcoding analyze, depending on a required minimal identity setup by
    the option 'identity'.
    """
    with open("scripts/locimerged." + scriptExt, "w") as out13:
        out13.write(mainStreamMessage(f'runloci:'))
        for loci1b in loci1s:
            for sample in samples:
                out13.write(mainStreamMessage(f' {sample} VS {loci1b}...') +
                            f"vsearch --usearch_global ../infiles/{sample}_cluster_OK.fas --db ../refs/{loci1b}.fas"
                            f" --matched ../{loci1b}/{sample}_merged.fas --id {identity} --strand both\n")
        out13.write(mainStreamMessage(f'\n\n'))


def runloc_r1():
    """Similar to runloc_merged but with the R1 denoised and non-chimera sequences
    in case of amplicons > 550bp, with the VSEARCH command:
    vsearch --usearch_global ../infiles/SN_cluster_R1_OK.fas --db ../refs/L2.fas --matched ../L2/SN_R1.fas
     --id real --strand both

     L1 = locus for amplicons > 550 bp
     SN = sample name
     real = a real from 0 to 1, generally around 0.7
    """
    with open("scripts/locir1." + scriptExt, "w") as out14:
        out14.write(mainStreamMessage(f'lociR1:'))
        for locus2b in loci2s:
            for sample in samples:
                out14.write(mainStreamMessage(f' {sample} VS {locus2b}...') +
                            f"vsearch --usearch_global ../infiles/{sample}_cluster_R1_OK.fas --db ../refs/{locus2b}.fas"
                            f" --matched ../{locus2b}/{sample}_R1.fas --id {identity} --strand both\n")
        out14.write(mainStreamMessage(f'\n\n'))


# ------------------- DISTRIBUTION OF CLUSTERS FOR ONE SELECTED LOCUS -----------
def runlocsel_merged():
    """Similar to runloc_merged but with only one selected locus with the VSEARCH command:
    vsearch --usearch_global ../infiles/SN_cluster_OK.fas --db ../refs/L1.fas --matched
    ../L2/SN_merged.fas --id real --strand both

    SN = sample name
    L1 = locus for amplicons < 550 bp
    """
    with open("scripts/loci_sel." + scriptExt, "w") as out15:
        out15.write(mainStreamMessage(f'locsel:'))
        for sample in samples:
            out15.write(mainStreamMessage(f' {sample}...') +
                        f"vsearch --usearch_global ../infiles/{sample}_cluster_OK.fas --db ../refs/{loc_sel}.fas"
                        f" --matched ../{loc_sel}/{sample}_merged.fas --id {identity} --strand both\n")
        out15.write(mainStreamMessage(f'\n\n'))


def runlocsel_r1():
    """identical to runlocsel_merged but with R1 denoised and non-chimera sequences with the VSEARCH command:
    vsearch --usearch_global ../infiles/SN_cluster_R1_OK.fas --db ../refs/L2.fas --matched
    ../L2/SN_R1.fas --id real --strand both

    SN = sample name
    L2 = locus name for amplicons > 550 bp
    """
    with open("scripts/locir1_sel." + scriptExt, "w") as out16:
        out16.write(mainStreamMessage(f'loc_selR1:'))
        for sample in samples:
            out16.write(mainStreamMessage(f' {sample}...') +
                        f"vsearch --usearch_global ../infiles/{sample}_cluster_R1_OK.fas --db ../refs/{loc_sel}.fas"
                        f" --matched ../{loc_sel}/{sample}_R1.fas --id {identity} --strand both\n")
        out16.write(mainStreamMessage(f'\n\n'))


# ------------------- DISTRIBUTION OF CLUSTERS FOR ONE SELECTED SAMPLE -----------

def runloc_one_sample_6():
    """ Searches similarities between denoised and non-chimera sequences and your local reference database (db)
    by the VSEARCH command, but only for a selected sample:
    vsearch --usearch_global ../infiles/SS_cluster_OK.fas --db ../refs/L1.fas --matched ../L1/SS_merged.fas
     --id real --strand both
    vsearch --usearch_global ../infiles/SS_cluster_OK.fas --db ../refs/L2.fas --matched ../L2/SS_merged.fas
     --id real --strand both

    SN = locus name
    L1 = locus name for amplicons < 550 bp
    L2 = locus name for amplicons > 550 bp
    real = real from 0 to 1
    This is the second crucial step metabarcoding analyze, depending on a required minimal identity setup by
    the option 'identity'
    """
    with open("scripts/loci_merged_6." + scriptExt, "w") as out13:
        out13.write(mainStreamMessage(f'runloci:'))
        for loci1b in loci1s:
            out13.write(mainStreamMessage(f' {sam_sel}...') +
                        f"vsearch --usearch_global ../infiles/{sam_sel}_cluster_OK.fas --db ../refs/{loci1b}.fas"
                        f" --matched ../{loci1b}/{sam_sel}_merged.fas --id {identity} --strand both\n")
        out13.write(mainStreamMessage(f'\n\n'))

    with open("scripts/loci_R1_6." + scriptExt, "w") as out14:
        out13.write(mainStreamMessage(f'runlociR1:'))
        for locus2b in loci2s:
            out14.write(mainStreamMessage(f' {sam_sel}...') +
                        f"vsearch --usearch_global ../infiles/{sam_sel}_cluster_R1_OK.fas --db ../refs/{locus2b}.fas"
                        f" --matched ../{locus2b}/{sam_sel}_R1.fas --id {identity} --strand both\n")
            out14.write(mainStreamMessage(f'\n\n'))


# ----------------------------- ORIENTATION OF SEQ INTO FORWARD -------------------------
def orient_merged():
    """Orients all the merged sequences in the same direction (forward) with the following script:
    vsearch --orient ../L1/SN_merged.fas --db ../refs/L1.fas --fastaout ../L1/SN_orient.fas

    SN = locus name
    L1 = locus name for amplicons < 550 bp
    """
    with open("scripts/orientloc1." + scriptExt, "w") as out17:
        out17.write(mainStreamMessage(f'orient:'))
        for locus1b in loci1s:
            for sample in samples:
                out17.write(mainStreamMessage(f' {sample} VS {locus1b}...') +
                            f"vsearch --orient ../{locus1b}/{sample}_merged.fas --db ../refs/{locus1b}.fas --fastaout "
                            f"../{locus1b}/{sample}_orient.fas\n")
        out17.write(mainStreamMessage(f'\n\n'))


def orient_r1():
    """Orients all the R1 sequences in the same direction (forward) with the following script:
    vsearch --orient ../L2/SN_merged.fas --db ../refs/L2.fas --fastaout ../L2/SN_orient.fas

    SN = locus name
    L2 = locus name for amplicons > 550 bp
    """
    with open("scripts/orientloc2." + scriptExt, "w") as out18:
        out18.write(mainStreamMessage(f'orientR1:'))
        for locus2b in loci2s:
            for sample in samples:
                out18.write(mainStreamMessage(f' {sample} VS {locus2b}...') +
                            f"vsearch --orient ../{locus2b}/{sample}_R1.fas --db ../refs/{locus2b}.fas --fastaout "
                            f"../{locus2b}/{sample}_R1_orient.fas\n")
        out18.write(mainStreamMessage(f'\n\n'))


def orient_one_loc_4():
    """Orients all the merged sequences in the same direction (forward) for a selected locus < 550 bp
    with the following script:
    vsearch --orient ../SL/SN_merged.fas --db ../refs/SL.fas --fastaout ../SL/SN_orient.fas

    SL = selected locus
    SN = sample name
    """
    with open("scripts/orient_merged_4." + scriptExt, "w") as out17:
        out17.write(mainStreamMessage(f'orientR1:'))
        for sample in samples:
            out17.write(mainStreamMessage(f' {sample}...') +
                        f"vsearch --orient ../{loc_sel}/{sample}_merged.fas --db ../refs/{loc_sel}.fas --fastaout "
                        f"../{loc_sel}/{sample}_orient.fas\n")
        out17.write(mainStreamMessage(f'\n\n'))


def orient_one_loc_5():
    """Orients all the R1 sequences in the same direction (forward) for loci > 550 bp with the following script:
    vsearch --orient ../SL/SN_R1.fas --db ../refs/SL.fas --fastaout ../SL/SN_R1_orient.fas

    SL = selected locus
    SN = sample name
    """
    with open("scripts/orient_R1_5." + scriptExt, "w") as out18:
        out18.write(mainStreamMessage(f'orientR1:'))
        for sample in samples:
            out18.write(mainStreamMessage(f' {sample}...') +
                        f"vsearch --orient ../{loc_sel}/{sample}_R1.fas --db ../refs/{loc_sel}.fas --fastaout "
                        f"../{loc_sel}/{sample}_R1_orient.fas\n")
        out18.write(mainStreamMessage(f'\n\n'))


def orient_one_sample_6():
    """Orients all the merged and R1 sequences in the same direction (forward) for a selected sample with
    the following script:
    vsearch --orient ../L1/SS_merged.fas --db ../refs/L1.fas --fastaout ../L1/SS_orient.fas
    vsearch --orient ../L2/SS_merged.fas --db ../refs/L2.fas --fastaout ../L1/SS_orient.fas

    L1 = locus name for amplicons < 550 bp
    L2 = locus name for amplicons > 550 bp
    SS = selected sample
    """
    with open("scripts/orient_merged_6." + scriptExt, "w") as out17:
        out17.write(mainStreamMessage(f'loci1:'))
        for locus1b in loci1s:
            out17.write(mainStreamMessage(f' {sam_sel}...') +
                        f"vsearch --orient ../{locus1b}/{sam_sel}_merged.fas --db ../refs/{locus1b}.fas --fastaout "
                        f"../{locus1b}/{sam_sel}_orient.fas\n")
        out17.write(mainStreamMessage(f'\n\n'))

    with open("scripts/orient_R1_6." + scriptExt, "w") as out18:
        out18.write(mainStreamMessage(f'lociR1:'))
        for locus2b in loci2s:
            out18.write(mainStreamMessage(f' {sam_sel}...') +
                        f"vsearch --orient ../{locus2b}/{sam_sel}_merged.fas --db ../refs/{locus2b}.fas --fastaout "
                        f"../{locus2b}/{sam_sel}_R1_orient.fas\n")
        out18.write(mainStreamMessage(f'\n\n'))


# ------------------------- CREATION OF GLOBAL SCRIPTS -----------------------------------
def runall_1():
    """ Defines the run order for option 1
    """
    with open("scripts/runall1." + scriptExt, "w") as out19:
        out19.write("./infor1." + scriptExt + "\n"
                    "./infor2." + scriptExt + "\n"
                    + startLogRedirect('../infiles/res1.log') +
                    "./merging." + scriptExt + "\n"
                    "./fqtofas." + scriptExt + "\n"
                    "./derep." + scriptExt + "\n"
                    "./derep_r1." + scriptExt + "\n"
                    "./cluster." + scriptExt + "\n"
                    "./cluster_r1." + scriptExt + "\n"
                    "./chimera." + scriptExt + "\n"
                    "./chimera_r1." + scriptExt + "\n"
                    "./locimerged." + scriptExt + "\n"
                    "./locir1." + scriptExt + "\n"
                    "./orientloc1." + scriptExt + "\n"
                    "./orientloc2." + scriptExt + "\n"
                    + endLogRedirect())
    os.chdir('scripts')
    sys.stdout.write("The files are being processed by 'vsearch', be patient !")
    subprocess.run([shellCmd, "./runall1." + scriptExt])


def runall_2():
    """ Defines the run order for option 2
    """
    with open("scripts/runall2." + scriptExt, "w") as out19:
        out19.write(startLogRedirect('../infiles/res2.log') +
                    "./merging." + scriptExt + "\n"
                    "./fqtofas." + scriptExt + "\n"
                    "./derep." + scriptExt + "\n"
                    "./derep_r1." + scriptExt + "\n"
                    "./cluster." + scriptExt + "\n"
                    "./cluster_r1." + scriptExt + "\n"
                    "./chimera." + scriptExt + "\n"
                    "./chimera_r1." + scriptExt + "\n"
                    "./locimerged." + scriptExt + "\n"
                    "./locir1." + scriptExt + "\n"
                    "./orientloc1." + scriptExt + "\n"
                    "./orientloc2." + scriptExt + "\n"
                    + endLogRedirect())
    os.chdir('scripts')
    if not winOS:
            for file in os.listdir("."):
                os.chmod(file, 0o755)
    #sys.stdout.write("The files are being copied, then 'vsearch' starts, be patient !")
  
    subprocess.run([shellCmd, "./runall2." + scriptExt])


def runall_3():
    """ Defines the run order for option 3
    """
    with open("scripts/runall3." + scriptExt, "w") as out19:
        out19.write(startLogRedirect('../infiles/res3.log') +
                    "./cluster." + scriptExt + "\n"
                    "./cluster_r1." + scriptExt + "\n"
                    "./chimera." + scriptExt + "\n"
                    "./chimera_r1." + scriptExt + "\n"
                    "./locimerged." + scriptExt + "\n"
                    "./locir1." + scriptExt + "\n"
                    "./orientloc1." + scriptExt + "\n"
                    "./orientloc2." + scriptExt + "\n"
                    + endLogRedirect())
    os.chdir('scripts')
    sys.stdout.write("The files are being processed, be patient !")
    subprocess.run([shellCmd, "./runall3." + scriptExt])


def runall_4():
    """ Defines the run order for option 4
    """
    with open("scripts/runall4." + scriptExt, "w") as out19:
        out19.write(startLogRedirect('../infiles/res4.log') +
                    "./cluster." + scriptExt + "\n"
                    "./chimera." + scriptExt + "\n"
                    "./loci_sel." + scriptExt + "\n"
                    "./orient_merged_4." + scriptExt + "\n"
                    + endLogRedirect())
    os.chdir('scripts')
    sys.stdout.write("The files are being processed, be patient !")
    subprocess.run([shellCmd, "./runall4." + scriptExt])


def runall_5():
    """ Defines the run order for option 5
    """
    with open("scripts/runall5." + scriptExt, "w") as out19:
        out19.write(startLogRedirect('../infiles/res5.log') +
                    "./cluster_r1." + scriptExt + "\n"
                    "./chimera_r1." + scriptExt + "\n"
                    "./locir1_sel." + scriptExt + "\n"
                    "./orient_R1_5." + scriptExt + "\n"
                    + endLogRedirect())
    os.chdir('scripts')
    sys.stdout.write("The files are being processed, be patient !")
    subprocess.run([shellCmd, "./runall5." + scriptExt])


def runall_6():
    """ Defines the run order for option 6
    """
    with open("scripts/runall6." + scriptExt, "w") as out19:
        out19.write(startLogRedirect('../infiles/res6.log') +
                    "./cluster_one_sample_6." + scriptExt + "\n"
                    "./chimera_one_sample_6." + scriptExt + "\n"
                    "./loci_merged_6." + scriptExt + "\n"
                    "./loci_R1_6." + scriptExt + "\n"
                    "./orient_merged_6." + scriptExt + "\n"
                    "./orient_R1_6." + scriptExt + "\n"
                    + endLogRedirect())
    os.chdir('scripts')
    sys.stdout.write("The files are being processed, be patient !")
    subprocess.run([shellCmd, "./runall6." + scriptExt])


# -------------------ANALYSIS OF NUMBER OF READS CLUSTERS BY LOCI ETC.. -------
def nb_clus():
    """Calculates the number of resulting sequences (reads, merged, dereplicates and clusters) according
    to the selected options 'minsize', minseqlength, alpha parameter' and 'identity'.
    """
    os.chdir("../")  # return to the working directory
    with open(sample_user, "rt") as out32, open("infiles/counts.txt", "w") as out33:
        lignes = out32.readlines()
        for ligne in lignes:
            ligne = ligne.rstrip()
            # Nb READS Calculated on R1.fa
            r1fa = open("infiles/" + ligne + "_R1.fa", "rt")
            reads = r1fa.read()
            nb_reads = reads.count(">")
            # Nb merged .fa
            merged = open("infiles/" + ligne + ".fa", "rt")
            mgd = merged.read()
            nb_merged = mgd.count(">")
            # Nb derepmerged Calculated on _derep.fas
            derepfas = open("infiles/" + ligne + "_derep.fas", "rt")
            df = derepfas.read()
            nb_derpm = df.count("sample")
            # Nb clustermerged Calculated on _cluster.fas
            clusmerged = open("infiles/" + ligne + "_cluster.fas", "rt")
            clusmer = clusmerged.read()
            nb_clusm = clusmer.count("sample")
            # Nb clustermergedK Calculated on _cluster_OK.fas
            clusmergedok = open("infiles/" + ligne + "_cluster_OK.fas", "rt")
            clusmerok = clusmergedok.read()
            nb_clusmok = clusmerok.count("sample")
            # Nb derepR1 Calculated on _derep_R1.fas
            derepfasr1 = open("infiles/" + ligne + "_derep_R1.fas", "rt")
            dfr1 = derepfasr1.read()
            nb_derpr1 = dfr1.count("sample")
            # Nb clusterR1 Calculated on _cluter_R1.fas
            clusfasr1 = open("infiles/" + ligne + "_cluster_R1.fas", "rt")
            clusr1 = clusfasr1.read()
            nb_clusr1 = clusr1.count("sample")
            # Nb clusterR1OK Calculated on _cluter_R1_OK.fas
            clusfasr1ok = open("infiles/" + ligne + "_cluster_R1_OK.fas", "rt")
            clusr1ok = clusfasr1ok.read()
            nb_clusr1ok = clusr1ok.count("sample")
            out33.writelines(f"The sample {ligne} has:\n"
                             f"\t{nb_reads} reads\n"
                             f"\t{nb_merged} merged sequences\n"
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
    sequences (amplicons < 550 bp).
    """
    with open("infiles/nbseq_bylocmerged.txt", "w") as out34:
        for sample in samples:
            sample = sample.rstrip()
            for locus1 in loci1s:
                locus1 = locus1.rstrip()
                os.chdir(locus1)
                refs = open(sample + "_merged.fas")
                refsloc = refs.read()
                nb_ref = refsloc.count("sample")
                out34.writelines(f"The sample {sample} showed:\n"
                                 f"\t{nb_ref} clusters of merged sequences for the locus {locus1}\n")
                os.chdir("../")


def nb_seqbyloc_r1():
    """Calculates the number of relevant cluster sequences by locus according
    to the selected options 'minsize', minseqlength, alpha parameter' and 'identity' for R1
    sequences (amplicons > 550 bp).
    """
    with open("infiles/nbseq_bylocR1.txt", "w") as out34:
        for sample in samples:
            sample = sample.rstrip()
            for locus2 in loci2s:
                locus2 = locus2.rstrip()
                os.chdir(locus2)
                refs = open(sample + "_R1.fas")
                refsloc = refs.read()
                nb_ref = refsloc.count("sample")
                out34.writelines(f"The sample {sample} showed:\n"
                                 f"\t{nb_ref} clusters of R1 sequences for the locus {locus2}\n")
                os.chdir("../")


# --------------------- CLEANING FOLDERS AFTER ANALYSES AND OTHER IMPORTANT TOOLS-----------------
def clean():
    """Cleans the heaviest .fastq (if exist) and .fa files in the folder 'infiles'
    """
    os.chdir(current_dir)
    fastq_files = glob.glob('infiles/*.fastq')
    for fastq_file in fastq_files:
        try:
            os.remove(fastq_file)
        except OSError as e:
            sys.stdout.write(f"\nError:{e.strerror}")
    fa_files = glob.glob('infiles/*.fa')
    for fa_file in fa_files:
        try:
            os.remove(fa_file)
        except OSError as e:
            sys.stdout.write(f"\nError:{e.strerror}")


def blast():
    """BLAST the file of clusters of your choice against the 'nt' database of NCBI according to the command:
    NCBIWWW.qblast("blastn", "nt", 'yourfile', megablast=True, hitlist_size=3)

    NB: Needs Biopython and returns the three best hits for each sequence, usually enough to identify
    the sequence and creates a file-result named 'yourFN'_res_blast.txt.
    """
    from Bio.Blast import NCBIWWW
    from Bio import SearchIO
    clus = input("\nWhat cluster file would you blast to NCBI e.g. infiles/TUN27_cluster_OK.fas\n"
                 "or 16S/all_16S.fasta ? : ")
    #  BLAST TO NCBI
    fasta_string = open(clus).read()
    result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string, megablast=True, hitlist_size=3)
    with open(clus + "_.xml", "w") as out_handle:
        out_handle.write(result_handle.read())
        result_handle.close()
    #  PARSING BLAST RESULTS
    blast_qresults = SearchIO.parse(clus + "_.xml", "blast-xml")
    with open(clus + "_res_blast.txt", "w") as outb:
        for blast_qresult in blast_qresults:
            blast_hit = str(blast_qresult[0])
            blast_hit1 = str(blast_qresult[1])
            blast_hit2 = str(blast_qresult[2])
            outb.writelines(f"\n{blast_hit}\n{blast_hit1}\n{blast_hit2}\n")
    outb.close()


def importseq():
    """Imports sequences from NCBI according to their accession numbers. returns two files, one in 'fasta' format,
    the other in GenBank format

    NB: needs Biopython
    """
    from Bio import Entrez
    mail = input("\nEnter you email, mandotory for free accesse at NCBI Entrez : ")
    seq_to_import = input("\nEnter the genbank code (gb accession number) of the sequence you need : e.g. AF359039:  ")
    Entrez.email = mail
    filename_fas = seq_to_import + ".fas"
    filename_gb = seq_to_import + ".gb"
    handle_fas = Entrez.efetch(db="nucleotide", id=seq_to_import, rettype="fasta", retmode="text")
    handle_gb = Entrez.efetch(db="nucleotide", id=seq_to_import, rettype="gb", retmode="text")
    out_handle_fas = open("infiles/" + filename_fas, "w")
    out_handle_fas.write(handle_fas.read())
    out_handle_gb = open("infiles/" + filename_gb, "w")
    out_handle_gb.write(handle_gb.read())
    out_handle_fas.close()
    out_handle_gb.close()

    sys.stdout.write("\nThe sequence is saved as *.fas (fasta format) and *.gb (genbank format) into the folder 'infiles'")


def trim_select_alls():
    """Select the final clusters of merged sequences by locus and by sample according to the size threshold retained
    for each of them.
    """
    dirloc = input("--------------------------------------------------------\n"
                   "Which LOCUS < 550bp do you want to analyze? \n"
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
                out20.writelines(f"" + startLogRedirect('./' + sample + '.log') +
                                 f"#Sum of sizes for {sample} = {a}'\n"
                                 f"#Threshold set to: {ts}'\n"
                                 f"#The sizes > {b} were conserved'\n"
                                 f"#Trimming left and right primers:'\n"
                                 f"#---------------------------------'\n"
                                 f"vsearch --fastx_filter {sample}_orient.fas --fastq_stripleft {trim_left} "
                                 f" --fastq_stripright {trim_right} --fastaout tmp --minsize {b}\n"
                                 f"#Dereplication after trimming:'\n"
                                 f"#---------------------------------'\n"
                                 f"vsearch --derep_fulllength  ./tmp --output {sample}_select.fas --sizein --sizeout\n"
                                 f"" + endLogRedirect())
            sys.stdout.write("----------------------")
            sys.stdout.write(f"\nSum of sizes for {sample} = {a}")
            sys.stdout.write(f"\nThe sizes > {b} were conserved")
            sys.stdout.write("----------------------")
            subprocess.run([shellCmd, "./trim-select." + scriptExt])
            os.remove("tmp")
            os.remove("trim-select." + scriptExt)
        os.chdir('../')
        dirloc = input("--------------------------------------------------------\n"
                       "Which new LOCUS < 550bp do you want to analyze? \n"
                       "Check that the corresponding folder does exists, e.g. nd4\n"
                       "Enter the NEW locus name or if you are done enter: 'end': ")


def trim_select_alls_r1():
    """Select the final clusters of merged sequences by locus and by sample according to the size threshold retained
    for each of them.
    """
    dirloc = input("--------------------------------------------------------\n"
                   "Which LOCUS > 550bp do you want to analyze? \n"
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
                out20.writelines(f"" + startLogRedirect('./{sample}.log') +
                                 f"#Sum of sizes for {sample} = {a}'\n"
                                 f"#Threshold set to: {ts}'\n"
                                 f"#The sizes > {b} were conserved'\n"
                                 f"#Trimming left and right primers:'\n"
                                 f"#---------------------------------'\n"
                                 f"vsearch --fastx_filter {sample}_R1_orient.fas --fastq_stripleft {trim_left} "
                                 f" --fastq_stripright {trim_right} --fastaout tmp --minsize {b}\n"
                                 f"#Dereplication after trimming:'\n"
                                 f"#---------------------------------'\n"
                                 f"vsearch --derep_fulllength  ./tmp --output {sample}_R1_select.fas --sizein --sizeout\n"
                                 f"" + endLogRedirect())
            sys.stdout.write("----------------------")
            sys.stdout.write(f"\nSum of sizes for {sample}= {a}")
            sys.stdout.write(f"\nThe sizes > {b} were conserved")
            sys.stdout.write("----------------------")
            subprocess.run([shellCmd, "./trim-select." + scriptExt])
            os.remove("tmp")
            os.remove("trim-select." + scriptExt)
        os.chdir('../')
        dirloc = input("--------------------------------------------------------\n"
                       "Which new LOCUS > 550bp do you want to analyze? \n"
                       "Check that the corresponding folder does exists, e.g. nd4\n"
                       "Enter the NEW locus name or if you are done enter: 'end': ")


def trim_select():
    """Select the final clusters of merged sequences by locus and by sample according to the size threshold retained
    for each of them.
    """
    dirloc = input("--------------------------------------------------------\n"
                   "Which LOCUS < 550bp do you want to analyze? \n"
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
                out20.writelines(f'' + startLogRedirect('./' + orient + '.log') +
                                 f'#Sum of sizes for {orient} = {a}\n'
                                 f'#Threshold set to: {ts}\n'
                                 f'#The sizes > {b} were conserved\n'
                                 f'#Trimming left and right primers:\n'
                                 f'#---------------------------------\n'
                                 f'vsearch --fastx_filter {orient}_orient.fas --fastq_stripleft {trim_left} '
                                 f' --fastq_stripright {trim_right} --fastaout tmp --minsize {b}\n'
                                 f'#Dereplication after trimming:\n'
                                 f'#---------------------------------\n'
                                 f'vsearch --derep_fulllength  ./tmp --output {orient}_select.fas --sizein --sizeout\n'
                                 f'' + endLogRedirect())
            sys.stdout.write("----------------------")
            sys.stdout.write(f"\nSum of sizes for {orient} = {a}")
            sys.stdout.write(f"\nThe sizes > {b} were conserved")
            sys.stdout.write("----------------------")
            subprocess.run([shellCmd, "./trim-select." + scriptExt])
            os.remove("tmp")
            os.remove("trim-select." + scriptExt)
            orient = input(f"--------------------------------------------------------\n"
                           f"Which NEW FILE do you want to trim for the locus {dirloc}?\n"
                           f"enter the sample name or, if you  terminated with the locus {dirloc} enter 'end': ")
        os.chdir('../')
        dirloc = input("--------------------------------------------------------\n"
                       "Which new LOCUS < 550bp do you want to analyze? \n"
                       "Check that the corresponding folder does exists, e.g. nd4\n"
                       "Enter the NEW locus name or if you are done enter: 'end': ")


def trim_select_r1():
    """Select the final clusters of merged sequences by locus and by sample according to the size threshold retained
    for each of them.
    """
    dirloc = input("--------------------------------------------------------\n"
                   "Which LOCUS > 550bp do you want to analyze? \n"
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
                out20.writelines(f'' + startLogRedirect('./' + orient + '.log') +
                                 f'#Sum of sizes for {orient} = {a}\n'
                                 f'#Threshold set to: {ts}\n'
                                 f'#The sizes > {b} were conserved\n'
                                 f'#Trimming left and right primers:\n'
                                 f'#---------------------------------\n'
                                 f'vsearch --fastx_filter {orient}_R1_orient.fas --fastq_stripleft {trim_left} '
                                 f' --fastq_stripright {trim_right} --fastaout tmp --minsize {b}\n'
                                 f'#Dereplication after trimming:\n'
                                 f'#---------------------------------\n'
                                 f'vsearch --derep_fulllength  ./tmp --output {orient}_R1_select.fas --sizein --sizeout\n'
                                 f'' + endLogRedirect())
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
                       "Which new LOCUS > 550bp do you want to analyze? \n"
                       "Check that the corresponding folder does exists, e.g. nd4\n"
                       "Enter the NEW locus name or if you are done enter: 'end': ")


def concat():
    """Concatenation of selected sequences according size threshold, locus by locus
    """

    loc2cat = input("----------------------------------------------------------\n"
                    "In which LOCUS do you want to concatenate the selected sequences? : ")
    while loc2cat != "end":
        while Path(loc2cat).is_dir() is False:
            loc2cat = input("\nYour locus name is not valid. Please enter a valid name for the locus : ")
        os.chdir(loc2cat)
        files2cat = glob.glob('*_select.fas')
        with open("./" + loc2cat + "_allseq_select.fasta", "w") as out:
            for file in files2cat:
                with open(file, "r") as out2:
                    out.write(out2.read())
        sys.stdout.write(f"\nThe files {files2cat} have been concatenated for locus {loc2cat}")
        os.chdir('../')
        loc2cat = input("----------------------------------------------------\n"
                        "Which NEW locus do you want concatenate?\n"
                        "Enter the new locus name or 'end' if you terminated: ")
    else:
        sys.stdout.write("\n**** YOUR CONCATENATION SESSION is COMPLETE ****\n")
        exit()


if __name__ == "__main__":
    date = datetime.datetime.now()
    current_dir = os.getcwd()
    rmenu = input("\n******************************  MAIN MENU  *****************************\n\n"
                  "1- NEW complete analyze (slow - be patient!)\n"
                  "2- NEW analyze WITHOUT the unnecessary statistical step (faster)\n"
                  "3- Re-analyze all loci, from the clustering step and with other parameters (faster)\n"
                  "4- Re-analyze only one locus of amplicon < 550bp, with other parameters\n"
                  "5- Re-analyze only one locus of amplicon > 550bp, with other parameters\n"
                  "6- Re-analyse only one sample, with other parameters\n"
                  "7- Clean unnecessary files '*.fastq' and '*.fa' files to free up space\n"
                  "8- Blast at NCBI a cluster of sequences for a sample\n"
                  "9- Retrieve a particular sequence from GenBank knowing its gb code (accession number)\n"
                  "10- Sub-menu for SELECTION OF MINIMUM SIZES according to THRESHOLDS\n"
                  "11- Concatenate all the files '*_select.fas' for a final analysing in MEGA7\n"
                  "12- Exit the program\n\n"
                  " *********  Type 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 or 12:  ")
    # NEW complete analyze
    if rmenu == "1":
        inputs_1_2()
        param_1_2()
        folder_1_2()
        stat()
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
        runall_1()
        nb_clus()
        nb_seqbyloc_merged()
        nb_seqbyloc_r1()

        sys.stdout.write("\n**** YOUR RUN OPTION 1 IS COMPLETE ****\n")

    # NEW analyze WITHOUT the unnecessary statistical step
    if rmenu == "2":
        inputs_1_2()
        param_1_2()
        folder_1_2()
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
        runall_2()
        nb_clus()
        nb_seqbyloc_merged()
        nb_seqbyloc_r1()

        sys.stdout.write("\n**** YOUR RUN OPTION 2 IS COMPLETE ****\n")

    # Re-analyze all loci, from the clustering step and with other parameters
    if rmenu == "3":
        prev_param()
        inputs_3()
        param_3()
        cluster_merged()
        cluster_r1()
        chimera_merged()
        chimera_r1()
        runloc_merged()
        runloc_r1()
        orient_merged()
        orient_r1()
        runall_3()
        nb_clus()
        nb_seqbyloc_merged()
        nb_seqbyloc_r1()

        sys.stdout.write("\n**** YOUR RUN OPTION 3 IS COMPLETE ****\n")

    # Re-analyze only one locus of amplicon < 450bp, with other parameters
    if rmenu == "4":
        prev_param()
        inputs_4_5()
        param_4()
        cluster_merged()
        chimera_merged()
        runlocsel_merged()
        orient_one_loc_4()
        runall_4()

        sys.stdout.write("\n**** YOUR RUN OPTION 4 IS COMPLETE ****\n")

    # Re-analyze only one locus of amplicon > 450bp
    if rmenu == "5":
        prev_param()
        inputs_4_5()
        param_5()
        cluster_r1()
        chimera_r1()
        runlocsel_r1()
        orient_one_loc_5()
        runall_5()

        sys.stdout.write("\n**** YOUR RUN OPTION 5 IS COMPLETE ****\n")

    # Re-analyse only one sample with other parameters
    if rmenu == "6":
        prev_param()
        inputs_6()
        param_6()
        cluster_one_sample_6()
        chimera_one_sample_6()
        runloc_one_sample_6()
        orient_one_sample_6()
        runall_6()

        sys.stdout.write("\n**** YOUR RUN OPTION 6 IS COMPLETE ****\n")

    # Clean unnecessary files '*.fastq' and '*.fa' files to free up space
    if rmenu == "7":
        clean()
        sys.stdout.write("\n**** YOU HAVE FREED UP SPACE IN YOUR HARD DRIVE ****\n")

    # Blast at NCBI a cluster of sequences for a sample
    if rmenu == "8":
        sys.stdout.write("BE PATIENT, the NCBI reply can be slow !! Maybe an online search should be faster")
        blast()
        sys.stdout.write("\nResults of your BLAST was returned in a file named 'yourfile'_res_blas.txt\n")

    # Retrieve a particular sequence from GenBank knowing its gb code (accession number)
    if rmenu == "9":
        importseq()
        sys.stdout.write("\nResults of your request were returned in two files 'yourGB'.fas and 'yourGB.gb\n")

    # Trim primers and select the size threshold
    if rmenu == '10':
        submenu = input("\n***************** MININIMUM SIZE SELECTION MENU *****************\n\n"
                        "1- Apply the same size threshold for all samples for locus < 550 bp (merged)\n"
                        "2- Apply the same size threshold for all samples for locus > 550 bp (R1)\n"
                        "3- Apply a specific size threshold for each sample, for locus < 550 bp (merged)\n"
                        "4- Apply a specific size threshold for each sample, for locus > 550 bp (R1)\n\n"
                        " ************************************** Type 1, 2, 3, 4 or end:  ")
        if submenu == "1":
            prev_param()
            trim_select_alls()
            sys.stdout.write("\n**** YOUR SIZE SELECTION for loci < 550bp is COMPLETE ****\n")

        if submenu == "2":
            prev_param()
            trim_select_alls_r1()
            sys.stdout.write("\n**** YOUR SIZE SELECTION for loci > 550bp is COMPLETE ****\n")

        if submenu == "3":
            prev_param()
            trim_select()
            sys.stdout.write("\n**** YOUR SIZE SELECTION sample by sample for loci < 550bp is COMPLETE ****\n")

        if submenu == "4":
            prev_param()
            trim_select_r1()
            sys.stdout.write("\n**** YOUR SIZE SELECTION sample by sample for loci > 550bp is COMPLETE ****\n")

        if submenu == "end":
            exit()

    # Concatenation of sample results
    if rmenu == "11":
        concat()

    # Exit the program
    if rmenu == "12":
        sys.stdout.write("\n**** THANK YOU FOR USING MBCTOOLS3 ****\n")
        exit()
