#!/usr/bin/python3

"""mbctools a cross-platform toolkit to make use of VSEARCH easier and interactive, thus helping analyze
metabarcoding data in the best conditions. It assumes VSEARCH is pre-installed and consists in the following MAIN MENU:

1 -> BASIC ANALYZES
2 -> REMOVAL OF PRIMERS AND SELECTION OF MINIMUM SEQUENCE ABUNDANCES ACCORDING TO USER-DEFINED THRESHOLDS
3 -> CONCATENATION OF ALL SAMPLE SELECTED SEQUENCES BY LOCUS FOR PHYLOGENETIC ANALYZES
4 -> CONVERSION OF ANALYSIS RESULTS INTO metaXplor IMPORT FORMAT

Option 1 offers the following submenu:

1  -> NEW COMPLETE ANALYSIS (mandatory)
1a -> Re-analyze all loci, from the clustering step, modifying parameters
1b -> Re-analyze only one locus of paired-end amplicon (merged reads), modifying parameters
1c -> Re-analyze only one locus of single-end amplicon (R1 only), modifying parameters
1d -> Re-analyse a given sample, modifying parameters
1e -> Optional quality checking of fastq files (slow)

Option 2 offers the following submenu:

2a -> Apply the SAME size threshold for ALL SAMPLES for the loci based on PAIRED-END reads (R1/R2 merged)
2b -> Apply the SAME size threshold for ALL SAMPLES for the loci based on SINGLE-END reads (R1 only)
2c -> Apply a SPECIFIC size threshold  for a given sample, for the loci based on PAIRED-END reads (R1/R2 merged)
2d -> Apply a SPECIFIC size threshold  for a given sample, for the loci based on SINGLE-END reads (R1 only)

Option 4 offers the following submenu:

4a -> Generate sequence files
4b -> Generate assignment file
4c -> Builds metaXplor-format sample metadata file from provided tabulated file
4d -> Compresses all metaXplor files into a final, ready to import, zip archive


VSEARCH reference:
Rognes T, Flouri T, Nichols B, Quince C, Mahe F (2016). VSEARCH: a versatile open source tool for metagenomics
PeerJ 4:e2584 doi: 10.7717/peerj.2584 https://doi.org/10.7717/peerj.2584

Usage:
=====
    python3 mbctools.py
"""

__authors__ = "Christian Barnabe, Guilhem Sempere"
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
import configparser
import traceback
import dateutil.parser
import io
import zipfile

try:
    p = subprocess.run(["vsearch"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
except FileNotFoundError:
    print("Please install vsearch, then try again...")
    exit(1)

winOS = "Windows" in platform.uname()
shellCmd = "powershell" if winOS else "bash"
scriptExt = "ps1" if winOS else "sh"
fileSep = "\\" if winOS else "/"
globalErrorOnStopCmd = "" if winOS else "set -e"
localErrorOnStopCmd = "; If ($LASTEXITCODE -gt 0) { exit $LASTEXITCODE }" if winOS else ""
date = datetime.datetime.now()
current_dir = os.getcwd()

global loci1s, loci2s, loci1, loci2, dir_fastq, fastq_R1, fastqr1s, fastq_R2, fastqr2s, Samples, \
    samples, alpha, identity, loc_sel1, loc_sel2, rmenu, sam_sel, minsize, minseqlength, loc2trim, trim_left, \
    trim_right, ts, ts1, sam2trim2c, sam2trim2d, alloci, next_run, menu, ts2, loc2cat, loc2trim2a, loc2trim2b, \
    loc2trim2c, loc2trim2d

errorStyle = "\033[91m"
warningStyle = "\033[93m"
normalStyle = "\033[0m"
titleStyle = "\033[94m\033[1m"
promptStyle = "\033[96m"
successStyle = "\033[92m"
citationStyle = "\033[1;33m\033[1m\033[3m"

metaXplorFasta = "metaXplor_sequences.fasta"
metaXplorSequenceComposition = "metaXplor_sequences.tsv"
metaXplorAssignments = "metaXplor_assignments.tsv"
metaXplorSamples = "metaXplor_samples.tsv"


def start_log_redirect(filepath):
    """Redirects log files if WinOs or not WinOS
    """
    if not winOS:
        return "exec 3>&1 4>&2 >" + filepath + " 2>&1\n"
    else:
        return "&{\n"


def end_log_redirect(filepath):
    """Redirects log files if WinOs or not WinOS
    """
    if not winOS:
        return "exec 1>&3 2>&4\n"
    else:
        return "} 2> " + filepath + "\n"


def main_stream_message(message):
    """Displays a main stream message on the console if WinOS or not WinOS
    """
    if not winOS:
        return "printf \"" + message + "\" >&3\n"
    else:
        return "Write-Host -NoNewline \"" + message + "\"\n"


def folders():
    """ Creates folders useful for the metabarcoding analyze
    """
    alreadyExisting = []

    global loci1s, loci2s, loci1, loci2
    sys.stdout.write("")
    folder = "scripts"
    path = os.path.join(current_dir, folder)
    if Path(path).is_dir():
        alreadyExisting.append(folder)
    else:
        os.mkdir(path)

    sys.stdout.write("")
    folder = "outputs"
    path = os.path.join(current_dir, folder)
    if Path(path).is_dir():
        alreadyExisting.append(folder)
    else:
        os.mkdir(path)

    sys.stdout.write("")
    folder = "refs"
    path = os.path.join(current_dir, folder)
    if Path(path).is_dir():
        alreadyExisting.append(folder)
    else:
        os.mkdir(path)

    sys.stdout.write("")
    folder = "tmp_files"
    path = os.path.join(current_dir, folder)
    if Path(path).is_dir():
        alreadyExisting.append(folder)
    else:
        os.mkdir(path)

    sys.stdout.write("")
    folder = "loci"
    path = os.path.join(current_dir, folder)
    if Path(path).is_dir():
        alreadyExisting.append(folder)
    else:
        os.mkdir(path)

    sys.stdout.write("")
    for locus in list(set(loci1s) | set(loci2s)):
        path = os.path.join(current_dir + "/loci", locus)
        if Path(path).is_dir():
            alreadyExisting.append("loci/" + locus)
        else:
            os.chdir(f"{current_dir}/loci")
            os.mkdir(path)

    if len(alreadyExisting) > 0:
        sys.stdout.write(
            f"\nThe following folders already exist and will be used in the current analysis: " + ", ".join(
                alreadyExisting))


def in_dir_fastq():
    """Input of the path containing the fastq files, option 1
    """
    global dir_fastq
    dir_fastq = input(
        "\n\n" + promptStyle + "Enter the FULL PATH of the folder where fastq files are located" + normalStyle +
        f"\nDefault is here:\n"
        f"{current_dir}/fastq\n"
        f"enter path: ")
    while Path(dir_fastq).is_dir() is False and dir_fastq not in ['', 'end', 'home', 'exit']:
        dir_fastq = input(errorStyle + "\n--> WRONG INPUT: " + normalStyle + " path not valid, enter a valid path\n"
                                                                             "default is {current_dir}/fastq\n"
                                                                             "OR 'end' 'home' 'exit': ")
    else:
        if dir_fastq == "":
            dir_fastq = f"{current_dir}/fastq"
        elif dir_fastq == "end":
            del dir_fastq
            main_menu1()
        elif dir_fastq == "home":
            del dir_fastq
            main()
        elif dir_fastq == "exit":
            quit_mbctools()


def in_fastq_R1():
    """Input of the file name containing the R1 fastq file names, option 1
    """
    global fastq_R1, fastqr1s
    fastq_R1 = input("\n" + promptStyle + "Enter R1 fastq file name\n" + normalStyle +
                         "default = fastqR1.txt: ")
    while os.path.isfile(fastq_R1) is False and fastq_R1 not in ["", "end", "home", "exit"]:
        fastq_R1 = input(
            errorStyle + "\n--> WRONG INPUT: " + normalStyle + " file name is not valid, please enter a valid file name\n"
                                                               "default = fastqR1.txt\n"
                                                               "OR 'end' 'home' 'exit': ")
    else:
        if fastq_R1 == "":
            fastq_R1 = 'fastqR1.txt'
        elif fastq_R1 == "end":
            del fastq_R1
            main_menu1()
        elif fastq_R1 == "home":
            del fastq_R1
            main()
        elif fastq_R1 == "exit":
            quit_mbctools()
    with open(fastq_R1, "r") as out1:
        fastqr1s = out1.read().splitlines()


def in_fastq_R2():
    """Input of the file name containing the R2 fastq file names, option 1
    """
    global fastq_R2, fastqr2s
    fastq_R2 = input("\n" + promptStyle + "Enter R2 fastq file name\n" + normalStyle +
                         "default = fastqR2.txt: ")
    while os.path.isfile(fastq_R2) is False and fastq_R2 not in ["", "end", "home", "exit"]:
        fastq_R2 = input(
            errorStyle + "\n--> WRONG INPUT: " + normalStyle + " file name is not valid, please enter a valid file name\n"
                                                               "default = fastqR2.txt\n"
                                                               "OR 'end' 'home' 'exit': ")
    else:
        if fastq_R2 == '':
            fastq_R2 = 'fastqR2.txt'
        elif fastq_R2 == "end":
            del fastq_R2
            main_menu1()
        elif fastq_R2 == "home":
            del fastq_R2
            main()
        elif fastq_R2 == "exit":
            quit_mbctools()
    with open(fastq_R2, "r") as out2:
        fastqr2s = out2.read().splitlines()


def in_loci1():
    """Input of the file name containing the list of paired-end based loci, option 1
    """
    global loci1, loci1s
    loci1 = input(
        "\n" + promptStyle + "Enter the name of the file containing loci based on paired-end reads\n" + normalStyle +
        "default = locus1.txt: ")
    while os.path.isfile(loci1) is False and loci1 not in ["", "end", "home", "exit"]:
        loci1 = input(
            errorStyle + "\n--> WRONG INPUT: " + normalStyle + " file name is not valid, please enter a valid file name\n"
                                                               "default = locus1.txt\n"
                                                               "OR 'end' 'home' 'exit': ")
    else:
        if loci1 == '':
            loci1 = 'locus1.txt'
        elif loci1 == "end":
            del loci1
            main_menu1()
        elif loci1 == "home":
            del loci1
            main()
        elif loci1 == "exit":
            quit_mbctools()
    with open(loci1, "r") as out:
        loci1s = out.read().splitlines()
    return loci1s


def in_loci2():
    """Input of the file name containing the list of single-end based loci, option 1
    """
    global loci2, loci2s
    loci2 = input(
        "\n" + promptStyle + "Enter the name of the file containing loci based on single-end reads only (R1)\n" + normalStyle +
        "default = locus2.txt: ")
    while os.path.isfile(loci2) is False and loci2 not in ["", "end", "home", "exit"]:
        loci2 = input(
            errorStyle + "\n--> WRONG INPUT: " + normalStyle + " file name is not valid, please enter a valid file name\n"
                                                               "default = locus2.txt: "
                                                               "OR 'end' 'home' 'exit': ")
    else:
        if loci2 == "":
            loci2 = 'locus2.txt'
        elif loci2 == "end":
            del loci2
            main_menu1()
        elif loci2 == "home":
            del loci2
            main()
        elif loci2 == "exit":
            quit_mbctools()
    with open(loci2, "r") as out:
        loci2s = out.read().splitlines()
    return loci2s


def in_Samples():
    """Input of the file name containing the list of samples, option 1
    """
    global Samples, samples
    Samples = input("\n" + promptStyle + "Enter the name of the file containing sample names\n" + normalStyle +
                        "default = samples.txt: ")
    while os.path.isfile(Samples) is False and Samples not in ["", "end", "home", "exit"]:
        Samples = input(
            errorStyle + "\n--> WRONG INPUT: " + normalStyle + " file name is not valid, please enter a valid file name\n"
                                                               "default = samples.txt\n"
                                                               "OR 'end' 'home' 'exit': ")
    else:
        if Samples == "":
            Samples = 'samples.txt'
        elif Samples == "end":
            del Samples
            main_menu1()
        elif Samples == "home":
            del Samples
            main()
        elif Samples == "exit":
            quit_mbctools()
    with open(Samples, "r") as out5:
        samples = out5.read().splitlines()


def in_minsize():
    """Input of the minimum abundance of sequences to retain for denoising/clustering, options 1, 1a, 1b, 1c and 1d
    """
    global minsize
    minsize = input("\n" + promptStyle + "Enter the minsize option value for clusters,\n" + normalStyle +
                         "i.e. the minimum sequence abundance of the retained clusters\n"
                         "default = 8: ")
    while minsize.isnumeric() is False and minsize not in ["", "end", "home", "exit"]:
        minsize = input(errorStyle + "\n--> WRONG INPUT: " + normalStyle + " minsize option must be an integer\n"
                                                                                "enter an integer (e.g. 2, 10, 50...)\n"
                                                                                "default = 8\n"
                                                                                "OR 'end' 'home' 'exit': ")
    else:
        if minsize == '':
            minsize = '8'
        elif minsize == "end":
            del minsize
            main_menu1()
        elif minsize == "home":
            del minsize
            main()
        elif minsize == "exit":
            quit_mbctools()


def in_minseqlength():
    """Input of the minimum length of sequences to keep for any locus, options 1, 1a, 1b, 1c and 1d
    """
    global minseqlength
    minseqlength = input(
        "\n" + promptStyle + "Enter the minimum length of sequences to keep for any locus\n" + normalStyle +
        "default = 100: ")
    while minseqlength.isnumeric() is False and minseqlength not in ["", "end", "home", "exit"]:
        minseqlength = input(errorStyle + "\n--> WRONG INPUT: " + normalStyle + " minimum length must be an integer\n"
                                                                                "enter an integer (e.g. 100, 150, 180...)\n"
                                                                                "default = 100\n"
                                                                                "OR 'end' 'home' 'exit': ")
    else:
        if minseqlength == '':
            minseqlength = '100'
        elif minseqlength == "end":
            del minseqlength
            main_menu1()
        elif minseqlength == "home":
            del minseqlength
            main()
        elif minseqlength == "exit":
            quit_mbctools()


def in_alpha():
    """Input of the alpha parameter for denoising/clustering, options 1, 1a, 1b, 1c and 1d
    """
    global alpha
    alpha = input("\n" + promptStyle + "Enter alpha parameter for the clustering\n" + normalStyle +
                  "default = 2: ")
    while alpha.isnumeric() is False and alpha not in ["", "end", "home", "exit"]:
        alpha = input(errorStyle + "\n--> WRONG INPUT: " + normalStyle + " alpha parameter must be an integer\n"
                                                                         "enter an integer (e.g. 1, 2, 3...)\n"
                                                                         "default = 2\n"
                                                                         "OR 'end' 'home' 'exit': ")
    else:
        if alpha == '':
            alpha = '2'
        elif alpha == "end":
            del alpha
            main_menu1()
        elif alpha == "home":
            del alpha
            main()
        elif alpha == "exit":
            quit_mbctools()


def in_identity():
    """Input of the identity parameter (0, 1.0) to BLAST the clusters against reference sequences, for affiliating
    clusters to the different loci, options 1, 1a, 1b, 1c and 1d
    """
    global identity
    identity = input(
        "\n" + promptStyle + "Enter identity parameter to BLAST the clusters against references\n" + normalStyle +
        "i.e. the identity percentage, enter an integer from 0 to 100\n"
        "default = 70: ")
    while identity not in ["end", "home", "exit", ""] and identity.isnumeric() is False:
        identity = input(
            errorStyle + "\n--> WRONG INPUT: " + normalStyle + " identity parameter must be an integer from 0 to 100 \n"
                                                               "default = 70\n"
                                                               "OR 'end' 'home' 'exit': ")
    else:
        if identity == '':
            identity = '0.7'
        elif identity == "end":
            del identity
            main_menu1()
        elif identity == "home":
            del identity
            main()
        elif identity.isnumeric() and int(identity) < 100:
            identity = int(identity) / 100
        elif identity == "exit":
            quit_mbctools()


def in_loc_sel_merged():
    """Input of a selected locus based on paired-end reads to rerun for option 1b
    """
    global loc_sel1
    loc_sel1 = input(
        "\n" + promptStyle + "Enter the name of the locus analysed by paired-end reads you want to rerun\n" + normalStyle +
        f"among {loci1s}"
        f"no default: ")
    while loc_sel1 not in loci1s and loc_sel1 not in ["end", "home", "exit"]:
        loc_sel1 = input(
            errorStyle + "\n--> WRONG INPUT: " + normalStyle + " locus name is not valid, please enter valid name of locus\n"
                                                               f"among {loci1s}\n"
                                                               "OR 'end' 'home' 'exit': ")
    else:
        if loc_sel1 == "end":
            main_menu2()
        elif loc_sel1 == "home":
            main()
        elif loc_sel1 == "exit":
            quit_mbctools()


def in_loc_sel_r1():
    """Input of a selected locus based on single-end read (R1) to rerun for option 1c
    """
    global loc_sel2, rmenu
    loc_sel2 = input(
        "\n" + promptStyle + "Enter the name of the locus analysed by only single-end (R1) reads\n" + normalStyle +
        f"among {loci2s} you want to rerun\n"
        f"no default: ")
    while loc_sel2 not in loci2s and loc_sel2 not in ["end", "home", "exit"]:
        loc_sel2 = input(
            errorStyle + "\n--> WRONG INPUT: " + normalStyle + " locus name is not valid, please enter valid name of locus\n"
                                                               f"among {loci2s}\n"
                                                               "OR 'end' 'home' 'exit': ")
    else:
        if loc_sel2 == "end":
            main_menu2()
        elif loc_sel2 == "home":
            main()
        elif loc_sel2 == "exit":
            quit_mbctools()


def in_sam_sel():
    """Input of the sample name to rerun for option 1d
    """
    global sam_sel, samples
    sam_sel = input("\n" + promptStyle + f"Enter the sample name you want to rerun\n" + normalStyle +
                    f"among {samples}\n"
                    f"no default: ")
    while sam_sel not in samples and sam_sel not in ["end", "home", "exit"]:
        sam_sel = input(
            errorStyle + "\n--> WRONG INPUT: " + normalStyle + " sample name is not valid, please enter a valid sample name\n"
                                                               f"among {samples}\n"
                                                               f"no default\n"
                                                               "OR 'end' 'home' 'exit': ")
    else:
        if sam_sel == "end" and rmenu == "1d":
            main_menu1()
        elif sam_sel == "home":
            main()
        elif sam_sel == "exit:":
            quit_mbctools()


def in_loc2trim_2x():
    """Input of the loci names for selection of minium sequence abundances according to user-defined thresholds,
    options 2a, 2b, 2c and 2d
    """
    global loc2trim2a, loc2trim2b, loc2trim2c, loc2trim2d, loci1s, loci2s
    if rmenu == "2a":
        loc2trim2a = input(
            "\n" + promptStyle + "Enter a LOCUS name based on paired-end mergeable reads you want to analyze among:\n" + normalStyle +
            f"{loci1s}\n"
            "OR\n"
            "'end' to run another option 2x\n"
            "'home' to return to main menu\n"
            "'exit' to terminate current session: ")
        while loc2trim2a not in loci1s and loc2trim2a not in ["end", "home", "exit"]:
            loc2trim2a = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " locus name is not valid, please enter a valid name\n"
                                                                   f"among {loci1s}\n"
                                                                   "OR\n"
                                                                   "'end' to run an option 2a\n"
                                                                   "'home' to return to main menu\n"
                                                                   "'exit' to terminate current session: ")
        if loc2trim2a == "end":
            main_menu2()
        if loc2trim2a == "home":
            main()
        if loc2trim2a == "exit":
            sys.stdout.write("\nSelection of minimum sizes according to user-defined thresholds is achieved\n"
                             "Statistical summary for paired-end based loci:\n"
                             f"Results in ---> {current_dir}/outputs/Stats_option_2a.txt\n")
            quit_mbctools()
        return loc2trim2a

    if rmenu == "2b":
        loc2trim2b = input(
            "\n" + promptStyle + "Enter a LOCUS name based on single-end R1 reads you want to analyze among:\n" + normalStyle +
            f"{loci2s}\n"
            "OR\n"
            "'end' to run another option 2x\n"
            "'home' to return to main menu\n"
            "'exit' to terminate current session: ")
        while loc2trim2b not in loci2s and loc2trim2b not in ["end", "home", "exit"]:
            loc2trim2b = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " locus name is not valid, please enter a valid name\n"
                                                                   f"among {loci2s}\n"
                                                                   "OR\n"
                                                                   "'end' to run an option 2x\n"
                                                                   "'home' to return to main menu\n"
                                                                   "'exit' to terminate current session: ")
        if loc2trim2b == "end":
            main_menu2()
        if loc2trim2b == "home":
            main()
        elif loc2trim2b == "exit":
            sys.stdout.write("\nSelection of minimum sizes according to user-defined thresholds is achieved\n"
                             "Statistical summary for paired-end based loci:\n"
                             f"Results in ---> {current_dir}/outputs/Stats_option_2a.txt\n")
            quit_mbctools()
        return loc2trim2b

    if rmenu == "2c":
        loc2trim2c = input(
            "\n" + promptStyle + "Enter a LOCUS name based on paired-end mergeable reads you want to analyze among:\n" + normalStyle +
            f"{loci1s}\n"
            "OR\n"
            "'end' to run another option 2x\n"
            "'home' to return to main menu\n"
            "'exit' to terminate current session: ")
        while loc2trim2c not in loci1s and loc2trim2c not in ["end", "home", "exit"]:
            loc2trim2c = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " locus name is not valid, please enter a valid name\n"
                                                                   f"among {loci1s}\n"
                                                                   "OR\n"
                                                                   "'end' to run an option 2x\n"
                                                                   "'home' to return to main menu\n"
                                                                   "'exit' to terminate current session: ")
        if loc2trim2c == "end":
            main_menu2()
        if loc2trim2c == "home":
            main()
        if loc2trim2c == "exit":
            sys.stdout.write("\nSelection of minimum sizes according to user-defined thresholds is achieved\n"
                             "Statistical summary for paired-end based loci:\n"
                             f"Results in ---> {current_dir}/outputs/Stats_option_2c.txt\n")
            quit_mbctools()
        return loc2trim2c

    if rmenu == "2d":
        loc2trim2d = input(
            "\n" + promptStyle + "Enter a LOCUS name based on single-end (R1) reads you want to analyze among:\n" + normalStyle +
            f"{loci2s}\n"
            "OR\n"
            "'end' to run another option 2x\n"
            "'home' to return to main menu\n"
            "'exit' to terminate current session: ")
        while loc2trim2d not in loci2s and loc2trim2d not in ["end", "home", "exit"]:
            loc2trim2d = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " locus name is not valid, please enter a valid name\n"
                                                                   f"among {loci2s}\n"
                                                                   "OR\n"
                                                                   "'end' to run an option 2x\n"
                                                                   "'home' to return to main menu\n"
                                                                   "'exit' to terminate current session: ")
        if loc2trim2d == "end":
            main_menu2()
        if loc2trim2d == "home":
            main()
        if loc2trim2d == "exit":
            sys.stdout.write("\nSelection of minimum sizes according to user-defined thresholds is achieved\n"
                             "Statistical summary for paired-end based loci:\n"
                             f"Results in ---> {current_dir}/outputs/Stats_option_2d.txt\n")
            quit_mbctools()
        return loc2trim2d


def in_trim_sample2c():
    """Input of sample names for option 2c
    """
    global samples, sam2trim2c
    sam2trim2c = input("\n" + promptStyle + f"Enter SAMPLE do you want to trim among {samples}\n" + normalStyle +
                       f"among {samples}?\n"
                       f"If you finished with locus {loc2trim2c} enter 'end': ")
    while sam2trim2c not in samples and sam2trim2c not in ["end", "home", "exit"]:
        sam2trim2c = input(
            errorStyle + f"\n--> WRONG INPUT: " + normalStyle + " sample name '{sam2trim2c}' is not valid, please enter a valid name\n"
                                                                f"among {samples}\n"
                                                                "OR 'end' 'home' 'exit': ")
    else:
        if sam2trim2c == "end":
            trim_2x()
        if sam2trim2c == "home":
            main_menu2()
        if sam2trim2c == "exit":
            quit_mbctools()
    return sam2trim2c


def in_trim_sample2d():
    """Input of sample names for option 2d
    """
    global samples, sam2trim2d
    sam2trim2d = input("\n" + promptStyle + f"Enter SAMPLE do you want to trim among {samples}\n" + normalStyle +
                       f"among {samples}?\n"
                       f"If you finished with locus {loc2trim2d} enter 'end': ")
    while sam2trim2d not in samples and sam2trim2d not in ["end", "home", "exit"]:
        sam2trim2d = input(
            errorStyle + f"\n--> WRONG INPUT: " + normalStyle + " sample name '{sam2trim2d}' is not valid, please enter a valid name\n"
                                                                f"among {samples}\n"
                                                                "OR 'end' 'home' 'exit': ")
    else:
        if sam2trim2d == "end":
            trim_2x()
        if sam2trim2d == "home":
            main_menu2()
        if sam2trim2d == "exit":
            quit_mbctools()
    return sam2trim2d


def in_trim_left():
    """Input of the number of bp corresponding to the left primer to remove from the clusters,
    options 2a, 2b, 2c and 2d
    """
    global trim_left, loc2trim2a, loc2trim2b, loc2trim2c
    if rmenu == "2a":
        trim_left = input(
            "\n" + promptStyle + f"Enter the number of bp of the left primer for {loc2trim2a}? (e.g. 20): " + normalStyle)
        while trim_left.isnumeric() is False and trim_left not in ["end", "home", "exit"]:
            trim_left = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " enter an integer corresponding to the length of the left primer (e.g. 20)\n"
                                                                   "OR 'end' 'home' 'exit': ")
        else:
            if trim_left == "end":
                main_menu2()
            elif trim_left == "home":
                main()
            elif trim_left == "exit":
                quit_mbctools()

    if rmenu == "2b":
        trim_left = input(
            "\n" + promptStyle + f"Enter the number of bp of the left primer for {loc2trim2b}? (e.g. 20): " + normalStyle)
        while trim_left.isnumeric() is False and trim_left not in ["end", "home", "exit"]:
            trim_left = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " enter an integer corresponding to the length of the left primer (e.g. 20)\n"
                                                                   "OR 'end' 'home' 'exit': ")
        else:
            if trim_left == "end":
                main_menu2()
            elif trim_left == "home":
                main()
            elif trim_left == "exit":
                quit_mbctools()

    if rmenu == "2c":
        trim_left = input(
            "\n" + promptStyle + f"Enter the number of bp of the left primer for {loc2trim2c}? (e.g. 20): " + normalStyle)
        while trim_left.isnumeric() is False and trim_left not in ["end", "home", "exit"]:
            trim_left = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " enter an integer corresponding to the length of the left primer (e.g. 20)\n"
                                                                   "OR 'end' 'home' 'exit': ")
        else:
            if trim_left == "end":
                main_menu2()
            elif trim_left == "home":
                main()
            elif trim_left == "exit":
                quit_mbctools()

    if rmenu == "2d":
        trim_left = input(
            "\n" + promptStyle + f"Enter the number of bp of the left primer for {loc2trim2d}? (e.g. 20): " + normalStyle)
        while trim_left.isnumeric() is False and trim_left not in ["end", "home", "exit"]:
            trim_left = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " enter an integer corresponding to the length of the left primer (e.g. 20)\n"
                                                                   "OR 'end' 'home' 'exit': ")
        else:
            if trim_left == "end":
                main_menu2()
            elif trim_left == "home":
                main()
            elif trim_left == "exit":
                quit_mbctools()
    return trim_left


def in_trim_right():
    """Input of the number of bp corresponding to the right primer to remove from the clusters,
    options 2a, 2b, 2c and 2d
    """
    global loc2trim2a, loc2trim2b, loc2trim2c, loc2trim2d, trim_right
    if rmenu == "2a":
        trim_right = input(
            "\n" + promptStyle + f"Enter the number of bp of the right primer for {loc2trim2a}? (e.g. 22): " + normalStyle)
        while trim_right.isnumeric() is False and trim_right not in ["end", "home", "exit"]:
            trim_right = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " enter an integer corresponding to the length of the "
                                                                   "right primer (e.g. 22)\n"
                                                                   "OR 'end' 'home' 'exit': ")
        else:
            if trim_right == "end":
                main_menu2()
            elif trim_right == "home":
                main()
            elif trim_right == "exit":
                quit_mbctools()

    if rmenu == "2b":
        trim_right = input(
            "\n" + promptStyle + f"Enter the number of bp of the right primer for {loc2trim2b}?\n" + normalStyle +
            f"NB : may be 0 for single-end reads! (e.g. 22): ")
        while trim_right.isnumeric() is False and trim_right not in ["end", "home", "exit"]:
            trim_right = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " enter an integer corresponding to the length of the "
                                                                   "right primer (e.g. 22)\n"
                                                                   "!! May be 0 with single-end based loci"
                                                                   "OR 'end' 'home' 'exit': ")
        else:
            if trim_right == "end":
                main_menu2()
            elif trim_right == "home":
                main()
            elif trim_right == "exit":
                quit_mbctools()

    if rmenu == "2c":
        trim_right = input(
            "\n" + promptStyle + f"Enter the number of bp of the right primer for {loc2trim2c}? (e.g. 22): " + normalStyle)
        while trim_right.isnumeric() is False and trim_right not in ["end", "home", "exit"]:
            trim_right = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " enter an integer corresponding to the length of the "
                                                                   "right primer (e.g. 22)\n"
                                                                   "OR 'end' 'home' 'exit': ")
        else:
            if trim_right == "end":
                main_menu2()
            elif trim_right == "home":
                main()
            elif trim_right == "exit":
                quit_mbctools()

    if rmenu == "2d":
        trim_right = input(
            "\n" + promptStyle + f"Enter the number of bp of the right primer for {loc2trim2d}?\n" + normalStyle +
            f"NB: may be 0 for single-end reads (e.g. 22): ")
        while trim_right.isnumeric() is False and trim_right not in ["end", "home", "exit"]:
            trim_right = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " enter an integer corresponding to the length of the "
                                                                   "right primer (e.g. 22)\n"
                                                                   "!! May be 0 with single-end based loci"
                                                                   "OR 'end' 'home' 'exit': ")
        else:
            if trim_right == "end":
                main_menu2()
            elif trim_right == "home":
                main()
            elif trim_right == "exit":
                quit_mbctools()
    return trim_right


def in_ts():
    """Input of the user-defined threshold of minimum abundance of clusters to retain, options 2a, 2b, 2c and 2d
    """
    global ts, ts1, loc2trim2a, loc2trim2b, loc2trim2c, loc2trim2d
    if rmenu == "2a":
        ts = input(
            "\n" + promptStyle + f"Enter the THRESHOLD (integer between 0 and 100) you want to use for this locus {loc2trim2a}\n" + normalStyle +
            f"Example: if you want to keep only the clusters whose abundance (size) is greater than 5%\n"
            f"of the sum of sizes  for a given sample with {loc2trim2a}, enter 5\n"
            f"no default: ")
        while ts not in ["end", "home", "exit"] and ts.isnumeric() is False or ts == "":
            ts = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " THRESHOLD format is not valid, it must be an integer between 0 and 100\n"
                                                                   f"Enter a new THRESHOLD for {loc2trim2a}, (e.g. 5)\n"
                                                                   f"no default\n"
                                                                   "OR 'end' 'home' 'exit': ")
        else:
            if ts == "end":
                main_menu2()
            if ts == "home":
                main()
            if ts.isnumeric() and int(ts) < 100:
                ts1 = int(ts) / 100
            if ts == "exit":
                quit_mbctools()

    if rmenu == "2b":
        ts = input(
            "\n" + promptStyle + f"Enter the THRESHOLD (integer between 0 and 100) you want to use for this locus {loc2trim2b}\n" + normalStyle +
            f"Example: if you want to keep only the clusters whose abundance (size) is greater than 5%\n"
            f"of the sum of sizes  for a given sample with {loc2trim2b}, enter 5\n"
            f"no default: ")
        while ts not in ["end", "home", "exit"] and ts.isnumeric() is False or ts == "":
            ts = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " THRESHOLD format is not valid, it must be an integer between 0 and 100\n"
                                                                   f"Enter a new THRESHOLD for {loc2trim2b}, (e.g. 5)\n"
                                                                   f"no default\n"
                                                                   "OR 'end' 'home' 'exit': ")
        else:
            if ts == "end":
                main_menu2()
            if ts == "home":
                main()
            if ts.isnumeric() and int(ts) < 100:
                ts1 = int(ts) / 100
            if ts == "exit":
                quit_mbctools()

    if rmenu == "2c":
        ts = input(
            "\n" + promptStyle + f"Enter the THRESHOLD (integer between 0 and 100) you want to use for this sample {sam2trim2c}\n" + normalStyle +
            f"Example: if you want to keep only the clusters whose abundance (size) is greater than 5%\n"
            f"of the sum of sizes for {sam2trim2c}, enter 5\n"
            f"no default: ")
        while ts not in ["end", "home", "exit"] and ts.isnumeric() is False or ts == "":
            ts = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " THRESHOLD format is not valid, it must be an integer between 0 and 100\n"
                                                                   f"Enter a new THRESHOLD for {sam2trim2c}, (e.g. 5)\n"
                                                                   f"no default\n"
                                                                   "OR 'end' 'home' 'exit': ")
        else:
            if ts == "end":
                main_menu2()
            if ts == "home":
                main()
            if ts.isnumeric() and int(ts) < 100:
                ts1 = int(ts) / 100
            if ts == "exit":
                quit_mbctools()

    if rmenu == "2d":
        ts = input(
            "\n" + promptStyle + f"Enter the THRESHOLD (integer between 0 and 100) you want to use for this locus {loc2trim2d}\n" + normalStyle +
            f"Example: if you want to keep only the clusters whose abundance (size) is greater than 5%\n"
            f"of the sum of sizes  for a given sample with {loc2trim2d}, enter 5\n"
            f"no default: ")
        while ts not in ["end", "home", "exit"] and ts.isnumeric() is False or ts == "":
            ts = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " THRESHOLD format is not valid, it must be an integer between 0 and 100\n"
                                                                   f"Enter a new THRESHOLD for {loc2trim2d}, (e.g. 5)\n"
                                                                   f"no default\n"
                                                                   "OR 'end' 'home' 'exit': ")
        else:
            if ts == "end":
                main_menu2()
            if ts == "home":
                main()
            if ts.isnumeric() and int(ts) < 100:
                ts1 = int(ts) / 100
            if ts == "exit":
                quit_mbctools()
    return ts, ts1


def param_1x():
    """ Creates a file with one parameter by line, options 1, 1a, 1b, 1c and 1d
    """
    global dir_fastq, Samples, minsize, rmenu
    os.chdir(current_dir)

    config = configparser.ConfigParser()
    contents = {"date": date, "dir_fastq": dir_fastq, "fastq_R1": fastq_R1, "fastq_R2": fastq_R2,
                "minsize": minsize, "minseqlength": minseqlength, "alpha": alpha, "identity": identity}
    if rmenu == '1' or rmenu == '1a':
        contents["loci1"] = loci1
        contents["loci2"] = loci2
        contents["Samples"] = Samples
        contents["alpha"] = alpha
        contents["identity"] = identity
    elif rmenu == '1b':
        contents["loc_sel1"] = loci1
        contents["Samples"] = Samples
    elif rmenu == '1c':
        contents["loc_sel2"] = loc_sel2
        contents["Samples"] = Samples
    elif rmenu == '1d':
        contents["sam_sel"] = sam_sel

    config['mbctools'] = contents
    with open("outputs/parameters_option_" + rmenu + ".cfg", 'w') as configfile:
        config.write(configfile)


def prev_param(paramConfigFile):
    """ Recalls global variables for different options
    """
    try:
        global fastqr1s, fastqr2s, loci1s, loci2s, samples
        os.chdir(current_dir)

        fileToParse = paramConfigFile if paramConfigFile is not None else "outputs/parameters_option_1.cfg"
        config = configparser.ConfigParser()
        config.read(fileToParse)
        contents = config['mbctools']

        if paramConfigFile is not None:
            global minsize
            minsize = contents["minsize"]
            global minseqlength
            minseqlength = contents["minseqlength"]
            global alpha
            alpha = contents["alpha"]
            global identity
            identity = contents["identity"]

        global dir_fastq
        dir_fastq = contents["dir_fastq"]
        global fastq_R1
        fastq_R1 = contents["fastq_R1"]
        global fastq_R2
        fastq_R2 = contents["fastq_R2"]
        global loci1
        loci1 = contents["loci1"]
        global loci2
        loci2 = contents["loci2"]
        global Samples
        Samples = contents["Samples"]

        with open(fastq_R1, "r") as out1:
            fastqr1s = out1.read().splitlines()

        with open(fastq_R2, "r") as out2:
            fastqr2s = out2.read().splitlines()

        with open(loci1, "r") as out3:
            loci1s = out3.read().splitlines()

        with open(loci2, "r") as out4:
            loci2s = out4.read().splitlines()

        with open(Samples, "r") as out5:
            samples = out5.read().splitlines()

        return fastqr1s, fastqr2s, loci1s, loci2s, samples
    except KeyError:
        print(errorStyle + "\nMissing parameter in configuration file " + fileToParse + " - " + normalStyle)
        traceback.print_exc(limit=0)
        exit(1)


def quality():
    """Tests the quality of each 'fastq' file (option 1e) by the VSEARCH command:

    vsearch --fastq_eestats2 dir_fastq/fastqF-R1 --output ../outputs/SN_R1_quality.txt
    vsearch --fastq_eestats2 dir_fastq/fastqF-R2 --output ../outputs/SN_R2_quality.txt

    dir_fastq = directory containing the fastq files
    fastqF-R1 and  fastqF-R2 = fastq file names for the sample
    SN = sample name
    """
    global fastqr2s, fastqr1s
    with open("scripts/infor1." + scriptExt, "w") as out1:
        i = 0
        out1.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Quality statistical tests on R1 reads '
                                                                     f'for samples:\n'))
        while i < len(samples):
            sample = samples[i]
            fastqr1 = fastqr1s[i]
            i = i + 1
            out1.write(main_stream_message(f' {sample}...'))
            out1.write(f"vsearch --fastq_eestats2 {dir_fastq}{fileSep}{fastqr1} --output "
                       f"../outputs/{sample}_R1_quality.txt" + localErrorOnStopCmd + "\n")
        out1.write(main_stream_message(f'\n\n'))

    with open("scripts/infor2." + scriptExt, "w") as out2:
        i = 0
        out2.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Quality statistical tests on R2 reads '
                                                                     f'for samples:\n'))
        while i < len(samples):
            sample = samples[i]
            fastqr2 = fastqr2s[i]
            i = i + 1
            out2.write(main_stream_message(f' {sample}...'))
            out2.write(f"vsearch --fastq_eestats2 {dir_fastq}{fileSep}{fastqr2} --output "
                       f"../outputs/{sample}_R2_quality.txt" + localErrorOnStopCmd + "\n")
        out2.write(main_stream_message(f'\n\n'))


def merging():
    """Merges paired-end reads into one sequence, when the length of the expected amplicon allows it (option 1)
    according the VSEARCH command:

    vsearch --fastq_mergepairs dir_fastq/fastqR1 --reverse dir_fastq/fastqR2 --fastaout
    ../tmp_files/SN_pairedEnd.fa   --fastq_allowmergestagger --relabel sample=SN_merged.

    dir_fastq =  directory containing the fastq files
    fastqR1 = complete file name of fastq R1 read
    fastqR2 = fastqR2 complete name
    option  --fastq_allowmergestagger allows the merging of short fragments
    """
    with open("scripts/merging." + scriptExt, "w") as out:
        i = 0
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Merging paired-end reads  for a given sample:\n'))
        while i < len(samples):
            sample = samples[i]
            fastqr1 = fastqr1s[i]
            fastqr2 = fastqr2s[i]
            out.write(main_stream_message(
                f" {sample}...") +
                      f"vsearch --fastq_mergepairs {dir_fastq}{fileSep}{fastqr1} --reverse {dir_fastq}{fileSep}{fastqr2} "
                      f"--fastaout ../tmp_files/{sample}_pairedEnd.fa --fastq_allowmergestagger --relabel "
                      f"sample={sample}_merged." + localErrorOnStopCmd + "\n")
            i = i + 1
        out.write(main_stream_message(f'\n\n'))


def fastq2fas():
    """When the merging R1/R2 is impossible because of an unadapted size of amplicon, the reads R1 of 301 bp
    (better than R2) are used to search the relevant sequences.
    First, all R1 'fastq' files have to be transformed into 'fasta' files by the VSEARCH command:

    vsearch --fastq_filter dir_fastq/fastaqR1 –fastaout ../tmp_files/SN_singleEnd.fa

    Where : dir_fastq = directory containing the fastq files ; SN = sample name ; fastqR1 = name of the 'fastq' file
    containing the R1 reads
    """
    with open("scripts/fqtofas." + scriptExt, "w") as out:
        i = 0
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Converting FASTQ files into FASTA format '
                                                                    f'for samples:\n'))
        while i < len(samples):
            sample = samples[i]
            fastqr1 = fastqr1s[i]
            out.write(main_stream_message(
                f' {sample}...') +
                      f"vsearch --fastq_filter {dir_fastq}{fileSep}{fastqr1} --fastaout ../tmp_files/{sample}_singleEnd.fa --relabel "
                      f" sample={sample}_R1." + localErrorOnStopCmd + "\n")
            i = i + 1
        out.write(main_stream_message(f'\n\n'))


def derep_1():
    """Dereplicates merged sequences in a given 'fasta' file with the VSEARCH command:

    vsearch --fastx_uniques ../tmp_files/SN_pairedEnd.fa  --fastaout ../tmp_files/SN_pairedEnd_derep.fas --sizeout
    --strand both

    And dereplicates the R1 sequences in a given 'fasta' file with the VSEARCH command:

    vsearch --fastx_uniques ../tmp_files/SN_singleEnd.fa  --fastaout ../tmp_files/SN_singleEnd_derep.fas --sizeout
    --strand both

    Both commands dereplicate in both strands (option --strand both) and write abundance annotation (frequency)
    to output (option --sizeout).

    SN = sample name
    """
    with open("scripts/derep." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Dereplicating merged reads  for a given sample:\n'))
        for sample in samples:
            out.write(main_stream_message(
                f' {sample}...') +
                      f"vsearch --fastx_uniques ../tmp_files/{sample}_pairedEnd.fa  --fastaout "
                      f"../tmp_files/{sample}_pairedEnd_derep.fas --sizeout --strand both" + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))

    with open("scripts/derep_r1." + scriptExt, "w") as out1:
        out1.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Dereplicating single-end reads '
                                                                     f'for samples:\n'))
        for sample in samples:
            out1.write(main_stream_message(f' {sample}...') +
                       f"vsearch --fastx_uniques ../tmp_files/{sample}_singleEnd.fa --fastaout "
                       f"../tmp_files/{sample}_singleEnd_derep.fas --sizeout --strand both" + localErrorOnStopCmd
                       + "\n")
        out1.write(main_stream_message(f'\n\n'))


def cluster_1x():
    """Denoises and clusters Illumina dereplicated merged sequences and gives in output the centroids sequences
    to 'fasta' files (options 1, 1a, 1b, 1c and 1d) with the following VSEARCH commands:

    For loci based in paired-end reads:

    vsearch --cluster_unoise ../tmp_files/SN_pairedEnd_derep.fas --sizein --centroids
    ../tmp_files/SN(or SS)_pairedEnd_cluster.fas --strand both --minsize int --sizeout --unoise_alph int
    --minseqlength int

    For loci based in single-end reads:

    vsearch --cluster_unoise ../tmp_files/SN_singleEnd_derep.fas --sizein --centroids
    ../tmp_files/SN(or SS)_singleEnd_cluster.fas --strand both --minsize int --sizeout --unoise_alph int
    --minseqlength int

    SN = sample name; int = integer; SS = selected sample
    """
    global rmenu
    if rmenu in ["1", "1a", "1b"]:
        with open("scripts/cluster." + scriptExt, "w") as out:
            i = 0
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Clustering merged reads  for a given sample:\n'))
            while i < len(samples):
                sample = samples[i]
                out.write(main_stream_message(
                    f' {sample}...') +
                          f"vsearch --cluster_unoise ../tmp_files/{sample}_pairedEnd_derep.fas --sizein --centroids "
                          f"../tmp_files/{sample}_pairedEnd_cluster.fas --strand both --minsize {minsize} --sizeout "
                          f"--unoise_alph {alpha} --minseqlength {minseqlength}" + localErrorOnStopCmd + "\n")
                i = i + 1
            out.write(main_stream_message(f'\n\n'))

    if rmenu in ["1", "1a", "1c"]:
        with open("scripts/cluster_r1." + scriptExt, "w") as out:
            i = 0
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Clustering single-end reads '
                                                                        f' for a given sample:\n'))
            while i < len(samples):
                sample = samples[i]
                out.write(main_stream_message(
                    f' {sample}...') +
                          f"vsearch --cluster_unoise ../tmp_files/{sample}_singleEnd_derep.fas --sizein --centroids "
                          f"../tmp_files/{sample}_singleEnd_cluster.fas --strand both --minsize {minsize} --sizeout "
                          f"--unoise_alpha {alpha} --minseqlength {minseqlength}" + localErrorOnStopCmd
                          + "\n")
                i = i + 1
            out.write(main_stream_message(f'\n\n'))

    elif rmenu == "1d":
        with open("scripts/cluster_one_sample_1d." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(
                f'\nClustering reads for '
                f'selected sample {sam_sel}:') +
                      f"vsearch --cluster_unoise ../tmp_files/{sam_sel}_pairedEnd_derep.fas --sizein --centroids "
                      f"../tmp_files/{sam_sel}_pairedEnd_cluster.fas --strand both --minsize {minsize} --sizeout "
                      f"--unoise_alph {alpha} --minseqlength {minseqlength}" + localErrorOnStopCmd + "\n"
                                                                                                     f"vsearch --cluster_unoise ../tmp_files/{sam_sel}_singleEnd_derep.fas --sizein --centroids "
                                                                                                     f"../tmp_files/{sam_sel}_singleEnd_cluster.fas --strand both --minsize {minsize} "
                                                                                                     f"--sizeout --unoise_alph {alpha} --minseqlength {minseqlength}"
                      + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))


def chimera_remove():
    """Detects and removes potential chimeras in denoised merged sequences or single-end or selected sample
    (options 1, 1a, 1b, 1c and 1d), by the VSEARCH commands:

    vsearch --uchime3_denovo ../tmp_files/SN(or SS)_pairedEnd_cluster.fas --nonchimeras
    ../tmp_files/SN(or SS)_pairedEnd_cluster_OK.fas
    vsearch --uchime3_denovo ../tmp_files/SN(or SS)_singleEnd_cluster.fas --nonchimeras
    ../tmp_files/SN(or SS)_singleEnd_cluster_OK.fas

    SN = sample name; SS = selected sample
    """
    if rmenu in ["1", "1a", "1b"]:
        with open("scripts/chimera." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd +
                      "\n" + main_stream_message(f'Detecting and removing chimeras within merged reads of samples:\n'))
            for sample in samples:
                out.write(main_stream_message(
                    f' {sample}...') +
                          f"vsearch --uchime3_denovo ../tmp_files/{sample}_pairedEnd_cluster.fas --nonchimeras "
                          f"../tmp_files/{sample}_pairedEnd_cluster_OK.fas" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

    if rmenu in ["1", "1a", "1c"]:
        with open("scripts/chimera_r1." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Detecting and removing chimeras within '
                                                                        f'single-end reads of samples:\n'))
            for sample in samples:
                out.write(main_stream_message(
                    f' {sample}...') +
                          f"vsearch --uchime3_denovo ../tmp_files/{sample}_singleEnd_cluster.fas --nonchimeras"
                          f" ../tmp_files/{sample}_singleEnd_cluster_OK.fas" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

    if rmenu == "1d":
        with open("scripts/chimera_one_sample_1d." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(
                f'Detecting and removing chimeras for selected sample {sam_sel}') +
                      f"vsearch --uchime3_denovo ../tmp_files/{sam_sel}_pairedEnd_cluster.fas --nonchimeras "
                      f"../tmp_files/{sam_sel}_pairedEnd_cluster_OK.fas" + localErrorOnStopCmd +
                      "\n"
                      f"vsearch --uchime3_denovo ../tmp_files/{sam_sel}_singleEnd_cluster.fas --nonchimeras "
                      f"../tmp_files/{sam_sel}_singleEnd_cluster_OK.fas" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))


def runloc_merged():
    """Searches similarities between merged sequences, denoised and non-chimera sequences and the local reference
    database (-db), options 1 and 1a, by the VSEARCH command:

    vsearch --usearch_global ../tmp_files/SN_pairedEnd_cluster_OK.fas --db ../refs/L1.fas --matched
    ../loci/L1/SN_pairedEnd.fas --id int --strand both

    L1 = locus name for amplicons based on paired-end reads
    SN = sample name
    id = minimum identity accepted (0-1.0)
    """
    with open("scripts/locimerged." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd +
                  "\n" + main_stream_message(f'Affiliating clusters to loci for merged reads of samples:\n'))
        for loci1b in loci1s:
            for sample in samples:
                out.write(main_stream_message(
                    f' {sample} vs {loci1b}...') +
                          f"vsearch --usearch_global ../tmp_files/{sample}_pairedEnd_cluster_OK.fas --db ../refs/{loci1b}.fas"
                          f" --matched ../loci/{loci1b}/{sample}_pairedEnd.fas --id {identity} --strand both"
                          + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))


def runloc_r1():
    """Searches similarities between single-end sequences (in case of amplicons where the merging R1/R2 is impossible)
    denoised and non-chimera sequences and the local reference database (-db), options 1 and 1a,
    by the VSEARCH command:

     vsearch --usearch_global ../tmp_files/SN_singleEnd_cluster_OK.fas --db ../refs/L2.fas --matched
     ../loci/L2/SN_singledEnd.fa --id int --strand both

    L2 = locus for amplicons with no mergeable R1/R2
    SN = sample name
    id = minimum identity accepted (0-1.0)
    """
    with open("scripts/locir1." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd +
                  "\n" + main_stream_message(f'Affiliating clusters to loci for single-end reads of samples:\n'))
        for locus2b in loci2s:
            for sample in samples:
                out.write(main_stream_message(f' {sample} vs {locus2b}...') +
                          f"vsearch --usearch_global ../tmp_files/{sample}_singleEnd_cluster_OK.fas --db "
                          f"../refs/{locus2b}.fas --matched ../loci/{locus2b}/{sample}_singleEnd.fas --id {identity} "
                          f"--strand both" + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))


def runlocsel_merged():
    """Searches similarities between merged sequences, denoised and non-chimera sequences and the local reference
    database (-db), for a selected paired-end based locus, option 1b, by the VSEARCH command:

    vsearch --usearch_global ../tmp_files/SN_pairedEnd_cluster_OK.fas --db ../refs/SL.fas --matched
    ../loci/SL/SN_merged.fas --id real --strand both

    SN = sample name
    SL = Selected locus based on paired-end reads
    id = minimum identity accepted (0-1.0)
    """
    with open("scripts/loci_sel." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters to selected '
                                                                    f'locus {loc_sel1}  for a given sample:\n'))
        for sample in samples:
            out.write(main_stream_message(f' {sample}...') +
                      f"vsearch --usearch_global ../tmp_files/{sample}_pairedEnd_cluster_OK.fas --db "
                      f"../refs/{loc_sel1}.fas --matched ../loci/{loc_sel1}/{sample}_pairedEnd.fas --id {identity} "
                      f"--strand both"
                      + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))


def runlocsel_r1():
    """Searches similarities between single-end sequences (R1), denoised and non-chimera sequences and the local
    reference database (-db), for a selected single-end based locus, option 1c, by the VSEARCH command:

    vsearch --usearch_global ../tmp_files/SN_singleEnd_cluster_OK.fas --db ../refs/SL.fas --matched
    ../loci/SL/SN_singleEnd.fas --id real --strand both

    SN = sample name
    SL = Selected locus based on single-end read (R1)
    id = minimum identity accepted (0-1.0)
    """
    with open("scripts/locir1_sel." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters to selected '
                                                                    f'locus {loc_sel2}  for a given sample:\n'))
        for sample in samples:
            out.write(main_stream_message(f' {sample}...') +
                      f"vsearch --usearch_global ../tmp_files/{sample}_singleEnd_cluster_OK.fas --db "
                      f"../refs/{loc_sel2}.fas --matched ../loci/{loc_sel2}/{sample}_singleEnd.fas --id {identity} "
                      f"--strand both"
                      + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))


def runloc_one_sample_1d():
    """ Searches similarities between denoised and non-chimera sequences and local reference database (db)
    only for a selected sample (option 1d) by the VSEARCH command:

    vsearch --usearch_global ../tmp_files/SS_pairedEnd_cluster_OK.fas --db ../refs/L1.fas --matched
    ../loci/L1/SS_pairedEnd.fas --id real --strand both

    vsearch --usearch_global ../tmp_files/SS_singleEnd_cluster_OK.fas --db ../refs/L2.fas --matched
    ../loci/L2/SS_singleEnd.fas --id real --strand both

    SS = selected sample
    L1 = locus name for amplicons based on paired-end reads
    L2 = locus name for amplicons with no mergeable R1/R2
    id = minimum identity accepted (0-1.0)
    """
    with open("scripts/loci_merged_1d." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters of selected '
                                                                    f'sample {sam_sel} to the loci:\n'))
        for loci1b in loci1s:
            out.write(main_stream_message(f' {loci1b}...') +
                      f"vsearch --usearch_global ../tmp_files/{sam_sel}_pairedEnd_cluster_OK.fas --db "
                      f"../refs/{loci1b}.fas --matched ../loci/{loci1b}/{sam_sel}_pairedEnd.fas --id {identity} "
                      f"--strand both" + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))

    with open("scripts/loci_R1_1d." + scriptExt, "w") as out1:
        out1.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters fo selected '
                                                                     f'sample {sam_sel} to the loci:\n'))
        for locus2b in loci2s:
            out1.write(main_stream_message(f' {locus2b}...') +
                       f"vsearch --usearch_global ../tmp_files/{sam_sel}_singleEnd_cluster_OK.fas --db "
                       f"../refs/{locus2b}.fas --matched ../loci/{locus2b}/{sam_sel}_singleEnd.fas "
                       f"--id {identity} --strand both"
                       + localErrorOnStopCmd + "\n")
            out1.write(main_stream_message(f'\n\n'))


def orient_1x():
    """Orients all the sequences in the same direction (forward) than references, options 1, 1a, 1c and 1d,
    with the following script:

    For paired-end merged sequences:
    vsearch --orient ../loci/L1/SN_pairedEnd.fas --db ../refs/L1.fas --fastaout ../loci/L1/SN_pairedEnd_orient.fas

    For single-end R1 sequences:
    vsearch --orient ../loci/L2/SN_singleEnd.fas --db ../refs/L2.fas --fastaout ../loci/L2/SN_singleEnd_orient.fas

    For selected loci:
    vsearch --orient ../loci/SL/SN_pairedEnd.fas --db ../refs/SL.fas --fastaout ../loci/SL/SN_pairedEnd_orient.fas
    or
    vsearch --orient ../loci/SL/SN_singleEnd.fas --db ../refs/SL.fas --fastaout ../loci/SL/SN_singleEnd_orient.fas

    For selected sample:
    vsearch --orient ../loci/L1/SS_pairedEnd.fas --db ../refs/L1.fas --fastaout ../loci/L1/SS_pairedEnd_orient.fas
    and
    vsearch --orient ../loci/L2/SS_singleEnd.fas --db ../refs/L2.fas --fastaout ../loci/L1/SS_singleEnd_orient.fas

    SN = locus name; SS = selected sample; SL = selected locus
    L1 = locus name for amplicons based on paired-end reads
    L2 = locus name for amplicons with no mergeable R1/R2
    """
    if rmenu in ["1", "1a"]:
        with open("scripts/orientloc1." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" +
                      main_stream_message(f"Orienting all merged reads in the same direction  for a given sample:\n"))
            for locus1b in loci1s:
                for sample in samples:
                    out.write(main_stream_message(f' {sample} vs {locus1b}...') +
                              f"vsearch --orient ../loci/{locus1b}/{sample}_pairedEnd.fas --db ../refs/{locus1b}.fas "
                              f"--fastaout ../loci/{locus1b}/{sample}_pairedEnd_orient.fas"
                              + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

        with open("scripts/orientloc2." + scriptExt, "w") as out1:
            out1.write(globalErrorOnStopCmd + "\n" +
                       main_stream_message(f"Orienting all single-end reads in "
                                           f"the same direction  for a given sample:\n"))
            for locus2b in loci2s:
                for sample in samples:
                    out1.write(main_stream_message(f' {sample} vs {locus2b}...') +
                               f"vsearch --orient ../loci/{locus2b}/{sample}_singleEnd.fas --db ../refs/{locus2b}.fas "
                               f"--fastaout ../loci/{locus2b}/{sample}_singleEnd_orient.fas" + localErrorOnStopCmd
                               + "\n")
            out1.write(main_stream_message(f'\n\n'))

    elif rmenu == "1b":
        with open("scripts/orient_merged_1b." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" +
                      main_stream_message(f'Orientation of all merged reads in the same '
                                          f'direction for selected locus {loc_sel1} '
                                          f'and  for a given sample:\n'))
            for sample in samples:
                out.write(main_stream_message(f' {sample}...') +
                          f"vsearch --orient ../loci/{loc_sel1}/{sample}_pairedEnd.fas --db "
                          f"../refs/{loc_sel1}.fas --fastaout "
                          f"../loci/{loc_sel1}/{sample}_pairedEnd_orient.fas" +
                          localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

    elif rmenu == "1c":
        with open("scripts/orient_R1_1c." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" +
                      main_stream_message(f'Orientation of all single-end reads in the '
                                          f'same direction for selected locus '
                                          f'{loc_sel2} and  for a given sample:\n'))
            for sample in samples:
                out.write(main_stream_message(f' {sample}...') +
                          f"vsearch --orient ../loci/{loc_sel2}/{sample}_singleEnd.fas --db ../refs/{loc_sel2}.fas "
                          f"--fastaout ../loci/{loc_sel2}/{sample}_singleEnd_orient.fas" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

    if rmenu == "1d":
        with open("scripts/orient_merged_1d." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" +
                      main_stream_message(f'Orientation of all clusters of selected '
                                          f'sample {sam_sel} for the loci:\n'))
            for locus1b in loci1s:
                out.write(main_stream_message(f' {locus1b}...') +
                          f"vsearch --orient ../loci/{locus1b}/{sam_sel}_pairedEnd.fas --db ../refs/{locus1b}.fas --fastaout "
                          f"../loci/{locus1b}/{sam_sel}_pairedEnd_orient.fas" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

        with open("scripts/orient_R1_1d." + scriptExt, "w") as out18:
            out18.write(main_stream_message(f'Orientation of all clusters of selected sample {sam_sel} for '
                                            f'the loci:\n'))
            for locus2b in loci2s:
                out18.write(main_stream_message(f' {locus2b}...') +
                            f"vsearch --orient ../loci/{locus2b}/{sam_sel}_singleEnd.fas --db ../refs/{locus2b}.fas "
                            f"--fastaout ../loci/{locus2b}/{sam_sel}_singleEnd_orient.fas" + localErrorOnStopCmd + "\n")
            out18.write(main_stream_message(f'\n\n'))


def runs_1x():
    """Creation of global scripts for the options 1, 1a, 1b, 1c, Ad and 1e
    """
    global rmenu, p
    if rmenu == "1":
        with open(f"scripts/runall1.{scriptExt}", "w") as out:
            if not winOS:
                out.write(start_log_redirect('../outputs/res1.log') +
                          "scriptArray=('./merging." + scriptExt + "' './fqtofas." + scriptExt + "' './derep."
                          + scriptExt + "' './derep_r1." + scriptExt + "' './cluster." + scriptExt + "' './cluster_r1."
                          + scriptExt + "' './chimera." + scriptExt + "' './chimera_r1." + scriptExt + "' './locimerged."
                          + scriptExt + "' './locir1." + scriptExt + "' './orientloc1." + scriptExt + "' './orientloc2."
                          + scriptExt + "')\n" +
                          'for script in "${scriptArray[@]}"\n' +
                          '    do\n' +
                          '    if ! ${script}; then\n' +
                          '        printf "Error executing ${script}\\n\\n" >&3\n' +
                          '        exit 1\n' +
                          '    fi\n' +
                          'done\n' +
                          end_log_redirect('../outputs/res1.log'))
            else:
                out.write(start_log_redirect('../outputs/res1.log') +
                          '   $scriptArray = @("./merging.' + scriptExt + '", "./fqtofas.'
                          + scriptExt + '", "./derep.' + scriptExt + '", "./derep_r1.'
                          + scriptExt + '", "./cluster.' + scriptExt + '", "./cluster_r1.'
                          + scriptExt + '", "./chimera.' + scriptExt + '", "./chimera_r1.'
                          + scriptExt + '", "./locimerged.' + scriptExt + '", "./locir1.'
                          + scriptExt + '", "./orientloc1.' + scriptExt + '", "./orientloc2.'
                          + scriptExt + '")\n' +
                          '   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
                          '        $script = $scriptArray[$i]\n' +
                          '        & "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; '
                          'exit $LASTEXITCODE }\n' +
                          '   }\n' +
                          end_log_redirect('../outputs/res1.log'))
        os.chdir('scripts')
        if not winOS:
            for file in os.listdir("."):
                os.chmod(file, 0o755)
        p = subprocess.run([shellCmd, "./runall1." + scriptExt])
        if p.returncode > 0:
            print(f"\nMain analysis execution failed, please check {current_dir}/outputs/res1.log")
            exit(1)
        print(successStyle + f"Step {rmenu} ended successfully" + normalStyle)

    if rmenu == "1a":
        with open(f"scripts/runall1a.{scriptExt}", "w") as out:
            if not winOS:
                out.write(start_log_redirect('../outputs/res1a.log') + "scriptArray=('"
                                                                       "./cluster." + scriptExt + "' '"
                                                                                                  "./cluster_r1." + scriptExt + "' '"
                                                                                                                                "./chimera." + scriptExt + "' '"
                                                                                                                                                           "./chimera_r1." + scriptExt + "' '"
                                                                                                                                                                                         "./locimerged." + scriptExt + "' '"
                                                                                                                                                                                                                       "./locir1." + scriptExt + "' '"
                                                                                                                                                                                                                                                 "./orientloc1." + scriptExt + "' '"
                                                                                                                                                                                                                                                                               "./orientloc2." + scriptExt + "')\n" +
                          'for script in "${scriptArray[@]}"\n' +
                          '    do\n' +
                          '    if ! ${script}; then\n' +
                          '        printf "Error executing ${script}\\n\\n" >&3\n' +
                          '        exit 1\n' +
                          '    fi\n' +
                          'done\n' +
                          end_log_redirect('../outputs/res1a.log'))
            else:
                out.write(start_log_redirect('../outputs/res1a.log') +
                          '   $scriptArray = @("'
                          './cluster.' + scriptExt + '", "'
                                                     './cluster_r1.' + scriptExt + '", "'
                                                                                   './chimera.' + scriptExt + '", "'
                                                                                                              './chimera_r1.' + scriptExt + '", "'
                                                                                                                                            './locimerged.' + scriptExt + '", "'
                                                                                                                                                                          './locir1.' + scriptExt + '", "'
                                                                                                                                                                                                    './orientloc1.' + scriptExt + '", "'
                                                                                                                                                                                                                                  './orientloc2.' + scriptExt + '")\n' +
                          '   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
                          '        $script = $scriptArray[$i]\n' +
                          '        & "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; '
                          'exit $LASTEXITCODE }\n' +
                          '   }\n' +
                          end_log_redirect('../outputs/res1a.log'))
        os.chdir('scripts')
        if not winOS:
            for file in os.listdir("."):
                os.chmod(file, 0o755)
        p = subprocess.run([shellCmd, "./runall1a." + scriptExt])
        if p.returncode > 0:
            print(f"\nMain analysis execution failed, please check {current_dir}/outputs/res1a.log")
            exit(1)
        print(successStyle + f"Step {rmenu} ended successfully" + normalStyle)

    if rmenu == "1b":
        with open(f"scripts/runall1b.{scriptExt}", "w") as out:
            if not winOS:
                out.write(start_log_redirect('../outputs/res1b.log') + "scriptArray=('"
                                                                       "./cluster." + scriptExt + "' '"
                                                                                                  "./chimera." + scriptExt + "' '"
                                                                                                                             "./loci_sel." + scriptExt + "' '"
                                                                                                                                                         "./orient_merged_1b." + scriptExt + "')\n" +
                          'for script in "${scriptArray[@]}"\n' +
                          '    do\n' +
                          '    if ! ${script}; then\n' +
                          '        printf "Error executing ${script}\\n\\n" >&3\n' +
                          '        exit 1\n' +
                          '    fi\n' +
                          'done\n' +
                          end_log_redirect('../outputs/res1b.log'))
            else:
                out.write(start_log_redirect('../outputs/res1b.log') +
                          '   $scriptArray = @("'
                          './cluster.' + scriptExt + '", "'
                                                     './chimera.' + scriptExt + '", "'
                                                                                './loci_sel' + scriptExt + '", "'
                                                                                                           './orient_merged_1b.' + scriptExt + '")\n' +
                          '   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
                          '        $script = $scriptArray[$i]\n' +
                          '        & "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; '
                          'exit $LASTEXITCODE }\n' +
                          '   }\n' +
                          end_log_redirect('../outputs/res1b.log'))
        os.chdir('scripts')
        if not winOS:
            for file in os.listdir("."):
                os.chmod(file, 0o755)
        p = subprocess.run([shellCmd, "./runall1b." + scriptExt])
        if p.returncode > 0:
            print(f"\nMain analysis execution failed, please check {current_dir}/outputs/res1b.log")
            exit(1)
        print(successStyle + f"Step {rmenu} ended successfully" + normalStyle)

    if rmenu == "1c":
        with open(f"scripts/runall1c.{scriptExt}", "w") as out:
            if not winOS:
                out.write(start_log_redirect('../outputs/res1c.log') + "scriptArray=('"
                                                                       "./cluster_r1." + scriptExt + "' '"
                                                                                                     "./chimera_r1." + scriptExt + "' '"
                                                                                                                                   "./locir1_sel." + scriptExt + "' '"
                                                                                                                                                                 "./orient_R1_1c." + scriptExt + "')\n" +
                          'for script in "${scriptArray[@]}"\n' +
                          '    do\n' +
                          '    if ! ${script}; then\n' +
                          '        printf "Error executing ${script}\\n\\n" >&3\n' +
                          '        exit 1\n' +
                          '    fi\n' +
                          'done\n' +
                          end_log_redirect('../outputs/res1c.log'))
            else:
                out.write(start_log_redirect('../outputs/res1c.log') +
                          '   $scriptArray = @("'
                          './cluster_r1.' + scriptExt + '", "'
                                                        './chimera_r1.' + scriptExt + '", "'
                                                                                      './locir1_sel.' + scriptExt + '", "'
                                                                                                                    './orient_R1_1c.' + scriptExt + '")\n' +
                          '   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
                          '        $script = $scriptArray[$i]\n' +
                          '        & "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; '
                          'exit $LASTEXITCODE }\n' +
                          '   }\n' +
                          end_log_redirect('../outputs/res1c.log'))
        os.chdir('scripts')
        if not winOS:
            for file in os.listdir("."):
                os.chmod(file, 0o755)
        p = subprocess.run([shellCmd, "./runall1c." + scriptExt])
        if p.returncode > 0:
            print(f"\nMain analysis execution failed, please check {current_dir}/outputs/res1c.log")
            exit(1)
        print(successStyle + f"Step {rmenu} ended successfully" + normalStyle)

    if rmenu == "1d":
        with open(f"scripts/runall1d.{scriptExt}", "w") as out:
            if not winOS:
                out.write(start_log_redirect('../outputs/res1d.log') +
                          "scriptArray=('"
                          "./cluster_one_sample_1d." + scriptExt + "' '"
                                                                   "./chimera_one_sample_1d." + scriptExt + "' '"
                                                                                                            "./loci_merged_1d." + scriptExt + "' '"
                                                                                                                                              "./loci_R1_1d." + scriptExt + "' '"
                                                                                                                                                                            "./orient_merged_1d." + scriptExt + "' '"
                                                                                                                                                                                                                "./orient_R1_1d." + scriptExt + "')\n" +
                          'for script in "${scriptArray[@]}"\n' +
                          '    do\n' +
                          '    if ! ${script}; then\n' +
                          '        printf "Error executing ${script}\\n\\n" >&3\n' +
                          '        exit 1\n' +
                          '    fi\n' +
                          'done\n' +
                          end_log_redirect('../outputs/res1d.log'))
            else:
                out.write(start_log_redirect('../outputs/res1d.log') +
                          '   $scriptArray = @("'
                          './cluster_one_sample_1d.' + scriptExt + '", "'
                                                                   './chimera_one_sample_1d.' + scriptExt + '", "'
                                                                                                            './loci_merged_1d.' + scriptExt + '", "'
                                                                                                                                              './loci_R1_1d.' + scriptExt + '", "'
                                                                                                                                                                            './orient_merged_1d.' + scriptExt + '", "'
                                                                                                                                                                                                                './orient_R1_1d.' + scriptExt + '")\n' +
                          '   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
                          '        $script = $scriptArray[$i]\n' +
                          '        & "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; '
                          'exit $LASTEXITCODE }\n' +
                          '   }\n' +
                          end_log_redirect('../outputs/res1d.log'))
        os.chdir('scripts')
        if not winOS:
            for file in os.listdir("."):
                os.chmod(file, 0o755)
        p = subprocess.run([shellCmd, "./runall1d." + scriptExt])
        if p.returncode > 0:
            print(f"\nMain analysis execution failed, please check {current_dir}/outputs/res1d.log")
            exit(1)
        print(successStyle + f"Step {rmenu} ended successfully" + normalStyle)

    if rmenu == "1e":
        with open("scripts/runall1e." + scriptExt, "w") as out:
            if not winOS:
                out.write(start_log_redirect('../outputs/res1e.log') +
                          "scriptArray=('"
                          "./infor1." + scriptExt + "' '"
                                                    "./infor2." + scriptExt + "')\n" +
                          'for script in "${scriptArray[@]}"\n' +
                          '    do\n' +
                          '    if ! ${script}; then\n' +
                          '        printf "Error executing ${script}\\n\\n" >&3\n' +
                          '        exit 1\n' +
                          '    fi\n' +
                          'done\n' + end_log_redirect('../outputs/res1e.log'))
            else:
                out.write(start_log_redirect('../outputs/res1e.log') +
                          '   $scriptArray = @("'
                          './infor1.' + scriptExt + '", "'
                                                    './infor2.' + scriptExt + '")\n' +
                          '   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
                          '        $script = $scriptArray[$i]\n' +
                          '        & "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; exit $LASTEXITCODE }\n' +
                          '   }\n' + end_log_redirect('../outputs/res1e.log'))
        os.chdir('scripts')
        if not winOS:
            for file in os.listdir("."):
                os.chmod(file, 0o755)
        sys.stdout.write("Quality checking is being processed, it is slow, be patient!\n\n")
        p = subprocess.run([shellCmd, "./runall1e." + scriptExt])
        if p.returncode > 0:
            print(
                errorStyle + f"\nMain analysis execution failed, please check {current_dir}/outputs/res1e.log" + normalStyle)
            exit(1)
        print(successStyle + f"Step {rmenu} ended successfully" + normalStyle)
        sys.stdout.write(successStyle + "\nExecution of option 1e is complete\n" + normalStyle +
                         f"Statistical test files (*_quality.txt) are located at ---> {current_dir}/outputs\n\n")
        os.chdir(current_dir)


def stats_1x():
    """Calculates the number of resulting sequences (reads, merged, dereplicates and clusters) according
    to the selected options 'minsize', 'minseqlength', 'alpha parameter' and 'identity', options 1, 1a, 1b, 1c and 1d
    """
    if rmenu == "1":
        os.chdir(current_dir)
        print("\nComputing statistics after using option 1, wait...\n"
              f"Results in ---> {current_dir}/outputs/Stats_option_1.txt:")
        with open("outputs/Stats_option_1.txt", "w") as out:
            out.write("With option 1, parameters set to:\n\n"
                      f"Directory = {current_dir}\n"
                      f"Fastq R1 file name = {fastq_R1}\n"
                      f"Fastq R2 file name = {fastq_R2}\n"
                      f"Paired-end based loci = {loci1s}\n"
                      f"Single-end based (R1) loci = {loci2s}\n"
                      f"Samples = {samples}\n"
                      f"Minimum abundance for clusters = {minsize}\n"
                      f"Minimum length for sequences = {minseqlength}\n"
                      f"Alpha clustering parameter = {alpha}\n"
                      f"Identity for allocating clusters = {identity}\n\n")
            for sample in samples:
                sample = sample.rstrip()
                # Nb READS Calculated on singleEnd.fa
                r1fa = open(f"tmp_files/{sample}_singleEnd.fa", "rt")
                reads = r1fa.read()
                nb_reads = reads.count(">")
                # Nb merged _pairedEnd.fa
                merged = open(f"tmp_files/{sample}_pairedEnd.fa", "rt")
                mgd = merged.read()
                nb_merged = mgd.count(">")
                percent_merging = (nb_merged / nb_reads) * 100
                percent = '{:.2f}%'.format(percent_merging)
                # Number of dereplicated sequences calculated on _pairedEnd_derep.fas
                derepfas = open(f"tmp_files/{sample}_pairedEnd_derep.fas", "rt")
                df = derepfas.read()
                nb_derpm = df.count("sample")
                # Nb clustermerged Calculated on _cluster_pairedEnd.fas
                clusmerged = open(f"tmp_files/{sample}_pairedEnd_cluster.fas", "rt")
                clusmer = clusmerged.read()
                nb_clusm = clusmer.count("sample")
                # Number of merged clusters without chimeras calculated on _pairedEnd_cluster_OK.fas
                clusmergedok = open(f"tmp_files/{sample}_pairedEnd_cluster_OK.fas", "rt")
                clusmerok = clusmergedok.read()
                nb_clusmok = clusmerok.count("sample")
                # Number of dereplicated sequences (R1) calculated on _singleEnd_derep.fas
                derepfasr1 = open(f"tmp_files/{sample}_singleEnd_derep.fas", "rt")
                dfr1 = derepfasr1.read()
                nb_derpr1 = dfr1.count("sample")
                # Nb clusterR1 Calculated on _singleEnd_cluster.fas
                clusfasr1 = open(f"tmp_files/{sample}_singleEnd_cluster.fas", "rt")
                clusr1 = clusfasr1.read()
                nb_clusr1 = clusr1.count("sample")
                # Nb clusterR1OK Calculated on _singleEnd_cluster_OK.fas
                clusfasr1ok = open(f"tmp_files/{sample}_singleEnd_cluster_OK.fas", "rt")
                clusr1ok = clusfasr1ok.read()
                nb_clusr1ok = clusr1ok.count("sample")
                out.writelines(f"\nThe sample {sample} has:\n"
                               f"\t{nb_reads} reads\n"
                               f"\t{nb_merged} merged sequences\n"
                               f"\tThe percentage of merging is {percent}\n"
                               f"\t{nb_derpm} dereplicated merged sequences\n"
                               f"\t{nb_clusm} merged clusters\n"
                               f"\t{nb_clusmok} merged clusters without chimera (OK)\n\n"
                               f"\t{nb_derpr1} dereplicated R1 sequences\n"
                               f"\t{nb_clusr1} R1 clusters\n"
                               f"\t{nb_clusr1ok} R1 clusters without chimera (R1_OK)\n\n")
                for locus1 in loci1s:
                    os.chdir(f"loci/{locus1}")
                    refs = open(f"{sample}_pairedEnd.fas")
                    refsloc = refs.read()
                    nb_ref = refsloc.count("sample")
                    out.writelines(f"\t{nb_ref} clusters of merged sequences of {sample} affiliated "
                                   f"to locus {locus1}\n")
                    os.chdir(current_dir)
                for locus2 in loci2s:
                    os.chdir(f"loci/{locus2}")
                    refs2 = open(sample + "_singleEnd.fas")
                    refsloc2 = refs2.read()
                    nb_ref2 = refsloc2.count("sample")
                    out.writelines(f"\t{nb_ref2} clusters of single-end (R1) sequences of {sample} affiliated "
                                   f"to locus {locus2}\n")
                    os.chdir(current_dir)
            print(successStyle + "\nThe MAIN MANDATORY ANALYSIS option 1 is complete\n" + normalStyle)

    elif rmenu == "1a":
        os.chdir(current_dir)
        print("\nComputing statistics using option 1a, wait...\n"
              f"Results in ---> {current_dir}/outputs/Stats_option_1a.txt:")
        with open("outputs/Stats_option_1a.txt", "w") as out:
            out.write("With option 1a, parameters set to:\n\n"
                      f"Directory = {current_dir}\n"
                      f"Fastq R1 file name = {fastq_R1}\n"
                      f"Fastq R2 file name = {fastq_R2}\n"
                      f"Paired-end based loci = {loci1s}\n"
                      f"Single-end based (R1) loci = {loci2s}\n"
                      f"Samples = {samples}\n"
                      f"Minimum abundance for clusters = {minsize}\n"
                      f"Minimum length for sequences = {minseqlength}\n"
                      f"Alpha clustering parameter = {alpha}\n"
                      f"Identity for allocating clusters = {identity}\n\n")
            for sample in samples:
                sample = sample.rstrip()
                # Nb READS Calculated on singleEnd.fa
                r1fa = open(f"tmp_files/{sample}_singleEnd.fa", "rt")
                reads = r1fa.read()
                nb_reads = reads.count(">")
                # Nb merged _pairedEnd.fa
                merged = open(f"tmp_files/{sample}_pairedEnd.fa", "rt")
                mgd = merged.read()
                nb_merged = mgd.count(">")
                percent_merging = (nb_merged / nb_reads) * 100
                percent = '{:.2f}%'.format(percent_merging)
                # Number of dereplicated sequences calculated on _pairedEnd_derep.fas
                derepfas = open(f"tmp_files/{sample}_pairedEnd_derep.fas", "rt")
                df = derepfas.read()
                nb_derpm = df.count("sample")
                # Nb clustermerged Calculated on _cluster_pairedEnd.fas
                clusmerged = open(f"tmp_files/{sample}_pairedEnd_cluster.fas", "rt")
                clusmer = clusmerged.read()
                nb_clusm = clusmer.count("sample")
                # Number of merged clusters without chimeras calculated on _pairedEnd_cluster_OK.fas
                clusmergedok = open(f"tmp_files/{sample}_pairedEnd_cluster_OK.fas", "rt")
                clusmerok = clusmergedok.read()
                nb_clusmok = clusmerok.count("sample")
                # Number of dereplicated sequences (R1) calculated on _singleEnd_derep.fas
                derepfasr1 = open(f"tmp_files/{sample}_singleEnd_derep.fas", "rt")
                dfr1 = derepfasr1.read()
                nb_derpr1 = dfr1.count("sample")
                # Nb clusterR1 Calculated on _singleEnd_cluster.fas
                clusfasr1 = open(f"tmp_files/{sample}_singleEnd_cluster.fas", "rt")
                clusr1 = clusfasr1.read()
                nb_clusr1 = clusr1.count("sample")
                # Nb clusterR1OK Calculated on _singleEnd_cluster_OK.fas
                clusfasr1ok = open(f"tmp_files/{sample}_singleEnd_cluster_OK.fas", "rt")
                clusr1ok = clusfasr1ok.read()
                nb_clusr1ok = clusr1ok.count("sample")
                out.writelines(f"\nThe sample {sample} has:\n"
                               f"\t{nb_reads} reads\n"
                               f"\t{nb_merged} merged sequences\n"
                               f"\tThe percentage of merging is {percent}\n"
                               f"\t{nb_derpm} dereplicated merged sequences\n"
                               f"\t{nb_clusm} merged clusters\n"
                               f"\t{nb_clusmok} merged clusters without chimera (OK)\n\n"
                               f"\t{nb_derpr1} dereplicated R1 sequences\n"
                               f"\t{nb_clusr1} R1 clusters\n"
                               f"\t{nb_clusr1ok} R1 clusters without chimera (R1_OK)\n\n")
                for locus1 in loci1s:
                    os.chdir(f"./loci/{locus1}")
                    refs = open(f"{sample}_pairedEnd.fas")
                    refsloc = refs.read()
                    nb_ref = refsloc.count("sample")
                    out.writelines(f"\t{nb_ref} clusters of merged sequences of {sample} affiliated "
                                   f"to locus {locus1}\n")
                    os.chdir(current_dir)
                for locus2 in loci2s:
                    os.chdir(f"./loci/{locus2}")
                    refs2 = open(sample + "_singleEnd.fas")
                    refsloc2 = refs2.read()
                    nb_ref2 = refsloc2.count("sample")
                    out.writelines(f"\t{nb_ref2} clusters of single-end (R1) sequences of {sample} affiliated "
                                   f"to locus {locus2}\n")
                    os.chdir(current_dir)
            print(successStyle + "\nExecution of option 1a is complete\n" + normalStyle)

    elif rmenu == "1b":
        os.chdir(current_dir)
        print("\nComputing statistics after using option 1b, wait...\n"
              f"Results in ---> {current_dir}/outputs/Stats_option_1b.txt:")
        with open("outputs/Stats_option_1b.txt", "w") as out:
            out.write("With option 1b, parameters set to:\n\n"
                      f"Directory = {current_dir}\n"
                      f"Fastq R1 file name = {fastq_R1}\n"
                      f"Fastq R2 file name = {fastq_R2}\n"
                      f"Paired-end based loci = {loci1s}\n"
                      f"Single-end based (R1) loci = {loci2s}\n"
                      f"Samples = {samples}\n"
                      f"Minimum abundance for clusters = {minsize}\n"
                      f"Minimum length for sequences = {minseqlength}\n"
                      f"Alpha clustering parameter = {alpha}\n"
                      f"Identity for allocating clusters = {identity}\n\n"
                      f"The selected locus is {loc_sel1}\n"
                      f"For {loc_sel1}:\n")
            os.chdir(f"./loci/{loc_sel1}")
            for sample in samples:
                a = open(sample + "_pairedEnd.fas")
                b = a.read()
                c = b.count("sample")
                out.write(f"\t{sample} has {c} clusters\n")
            os.chdir(current_dir)
        print(successStyle + "\nExecution of option 1b is complete\n" + normalStyle)

    elif rmenu == "1c":
        os.chdir(current_dir)
        print("\nComputing statistics after using option 1c, wait...\n"
              f"Results in ---> {current_dir}/outputs/Stats_option_1c.txt")
        with open("outputs/Stats_option_1c.txt", "w") as out:
            out.write("With option 1c, parameters set to:\n\n"
                      f"Directory = {current_dir}\n"
                      f"Fastq R1 file name = {fastq_R1}\n"
                      f"Fastq R2 file name = {fastq_R2}\n"
                      f"Paired-end based loci = {loci1s}\n"
                      f"Single-end based (R1) loci = {loci2s}\n"
                      f"Samples = {samples}\n"
                      f"Minimum abundance for clusters = {minsize}\n"
                      f"Minimum length for sequences = {minseqlength}\n"
                      f"Alpha clustering parameter = {alpha}\n"
                      f"Identity for allocating clusters = {identity}\n\n"
                      f"The selected locus is {loc_sel2}\n"
                      f"For {loc_sel2}:\n")
            os.chdir(f"./loci/{loc_sel2}")
            for sample in samples:
                a = open(sample + "_singleEnd.fas")
                b = a.read()
                c = b.count("sample")
                out.write(f"\t{sample} has {c} clusters\n")
            os.chdir(current_dir)
        print(successStyle + "\nExecution of option 1c is complete\n" + normalStyle)

    elif rmenu == "1d":
        os.chdir(current_dir)
        print("\nComputing statistics after using option 1d, wait\n"
              f"Results in ---> {current_dir}/outputs/Stats_option_1d.txt:")
        with open("outputs/Stats_option_1d.txt", "w") as out:
            out.write("With option 1d, parameters set to:\n\n"
                      f"Directory = {current_dir}\n"
                      f"Fastq R1 file name = {fastq_R1}\n"
                      f"Fastq R2 file name = {fastq_R2}\n"
                      f"Paired-end based loci = {loci1s}\n"
                      f"Single-end based (R1) loci = {loci2s}\n"
                      f"Samples = {samples}\n"
                      f"Minimum abundance for clusters = {minsize}\n"
                      f"Minimum length for sequences = {minseqlength}\n"
                      f"Alpha clustering parameter = {alpha}\n"
                      f"Identity for allocating clusters = {identity}\n\n"
                      f"The selected sample is {sam_sel}\n")

            for locus1 in loci1s:
                os.chdir(f"./loci/{locus1}")
                a = open(sam_sel + "_pairedEnd.fas")
                b = a.read()
                c = b.count("sample")
                out.write(f"\t{sam_sel} has {c} clusters for locus {locus1} based on paired-end reads\n")
                os.chdir(current_dir)

            for locus2 in loci2s:
                os.chdir(f"./loci/{locus2}")
                a = open(sam_sel + "_singleEnd.fas")
                b = a.read()
                c = b.count("sample")
                out.write(f"\t{sam_sel} has {c} clusters for locus {locus2} based on single-end reads\n")
                os.chdir(current_dir)
        print(successStyle + "\nExecution of option 1d is complete\n" + normalStyle)


def trim_2x():
    """Removes primers and selects clusters according minimum abundance thresholds fixed by the user
    """
    global loc2trim2a, loc2trim2b, loc2trim2c, loc2trim2d, ts, ts1, trim_left, trim_right, sam2trim2c, sam2trim2d
    if rmenu == "2a":
        os.chdir(current_dir)
        stat_2a = open(f"outputs/Stats_option_2a.txt", "w")
        while True:
            loc2trim2a = in_loc2trim_2x()
            os.chdir(f"./loci/{loc2trim2a}")
            trim_left = in_trim_left()
            trim_right = in_trim_right()
            ts = in_ts()
            stat_2a.write(f"Locus {loc2trim2a} trimmed {trim_left} bp (left) and {trim_right} bp (right) "
                          f"with threshold set at {ts1}\n")
            for sample in samples:
                with open(f"{sample}_pairedEnd_orient.fas", 'r') as filin, \
                        open("trim-select." + scriptExt, "w") as out:
                    targets = [line for line in filin if "size" in line]
                    a = 0
                    for target in targets:
                        size = re.search('size=(.+?)$', target).group(1)
                        a = a + int(size)
                    b = int(a * float(ts1) + 1)
                    stat_2a.writelines(f"\tSum of sizes for {sample} = {a}\n"
                                       f"\tThe sizes > {b} for {sample} were retained\n")
                    out.write(f"" + start_log_redirect('./' + sample + '.log') +
                              f"vsearch --fastx_filter {sample}_pairedEnd_orient.fas --fastq_stripleft {trim_left} "
                              f"  --fastq_stripright {trim_right} --fastaout tmp --minsize {b}"
                              + localErrorOnStopCmd + "\n"
                                                      f"vsearch --derep_fulllength ./tmp --output {sample}_pairedEnd_select.fas "
                                                      f"--sizein --sizeout"
                              + localErrorOnStopCmd + "\n"
                                                      f"" + end_log_redirect('./' + sample + '.log'))
                subprocess.run([shellCmd, "./trim-select." + scriptExt])
                selected = open(f'./{sample}_pairedEnd_select.fas', 'r')
                nb_selected = selected.read().count('>')
                stat_2a.writelines(f"\tNumber of selected clusters for sample {sample}: {nb_selected}\n\n")
                sys.stdout.write(f"\nSum of sizes for {sample} at locus {loc2trim2a} = {a}\n"
                                 f"With threshold set at {ts}, sizes > {b} were retained\n"
                                 f"Number of selected clusters for sample {sample}: {nb_selected}\n")
            os.remove("trim-select." + scriptExt)
            os.remove("tmp")
            os.chdir(current_dir)

    if rmenu == "2b":
        os.chdir(current_dir)
        stat_2b = open(f"outputs/Stats_option_2b.txt", "w")
        while True:
            loc2trim2b = in_loc2trim_2x()
            os.chdir(f"./loci/{loc2trim2b}")
            trim_left = in_trim_left()
            trim_right = in_trim_right()
            ts = in_ts()
            stat_2b.write(f"Locus {loc2trim2b} trimmed {trim_left} bp (left) and {trim_right} bp (right) "
                          f"with threshold set at {ts1}\n")
            for sample in samples:
                with open(sample + "_singleEnd_orient.fas", 'r') as filin, \
                        open("trim-select." + scriptExt, "w") as out:
                    targets = [line for line in filin if "size" in line]
                    a = 0
                    for target in targets:
                        size = re.search('size=(.+?)$', target).group(1)
                        a = a + int(size)
                    b = int(a * float(ts1) + 1)
                    stat_2b.writelines(f"\tSum of sizes for {sample} = {a}\n"
                                       f"\tThe sizes > {b} for {sample} were retained\n")
                    out.write(f"" + start_log_redirect('./' + sample + '.log') +
                              f"vsearch --fastx_filter {sample}_singleEnd_orient.fas --fastq_stripleft {trim_left} "
                              f"  --fastq_stripright {trim_right} --fastaout tmp --minsize {b}"
                              + localErrorOnStopCmd + "\n"
                                                      f"vsearch --derep_fulllength ./tmp --output {sample}_singleEnd_select.fas "
                                                      f"--sizein --sizeout" + localErrorOnStopCmd + "\n"
                                                                                                    f"" + end_log_redirect(
                        './' + sample + '.log'))
                subprocess.run([shellCmd, "./trim-select." + scriptExt])
                selected = open('./' + sample + '_singleEnd_select.fas', 'r')
                nb_selected = selected.read().count('>')
                stat_2b.writelines(f"\tNumber of selected clusters for sample {sample}: {nb_selected}\n\n")
                sys.stdout.write(f"\nSum of sizes for {sample} at locus {loc2trim2b} = {a}\n"
                                 f"With threshold set at {ts}, sizes > {b} were retained\n"
                                 f"Number of selected clusters for sample {sample}: {nb_selected}\n")
            os.remove("trim-select." + scriptExt)
            os.remove("tmp")
            os.chdir(current_dir)

    if rmenu == "2c":
        os.chdir(current_dir)
        stat_2c = open(f"{current_dir}/outputs/Stats_option_2c.txt", "a")
        while True:
            loc2trim2c = in_loc2trim_2x()
            os.chdir(f"{current_dir}/loci/{loc2trim2c}")
            trim_left = in_trim_left()
            trim_right = in_trim_right()
            stat_2c.write(f"Locus {loc2trim2c} trimmed at {trim_left} bp (left) and {trim_right} bp (right)\n")
            while True:
                sam2trim2c = in_trim_sample2c()
                ts = in_ts()
                with open(sam2trim2c + "_pairedEnd_orient.fas", "r") as filin, \
                        open("trim-select." + scriptExt, "w") as filout:
                    targets = [line for line in filin if "size" in line]
                    a = 0
                    for target in targets:
                        size = re.search('size=(.+?)$', target).group(1)
                        a = a + int(size)
                    b = int(a * float(ts1) + 1)
                    stat_2c.write(f"\tSum of sizes for {sam2trim2c} at locus {loc2trim2c} = {a}\n"
                                  f"\tAt threshold {ts} sizes > {b} for {sam2trim2c} were retained\n")
                    filout.writelines(f'' + start_log_redirect('./' + loc2trim2c + '.log') +
                                      f' vsearch --fastx_filter {sam2trim2c}_pairedEnd_orient.fas --fastq_stripleft {trim_left} '
                                      f'  --fastq_stripright {trim_right} --fastaout tmp --minsize {b}\n'
                                      + localErrorOnStopCmd + "\n"
                                                              f' vsearch --derep_fulllength ./tmp --output {sam2trim2c}_pairedEnd_select.fas --sizein '
                                                              f'--sizeout\n' + localErrorOnStopCmd + "\n"
                                                                                                     f'' + end_log_redirect(
                        './' + sam2trim2c + '.log'))
                subprocess.run([shellCmd, "./trim-select." + scriptExt])
                selected = open("./" + sam2trim2c + "_pairedEnd_select.fas", "r")
                nb_selected = selected.read().count(">")
                stat_2c.write(f"\tNumber of selected clusters for {sam2trim2c} is: {nb_selected}\n\n")
                sys.stdout.write(f"\n{sam2trim2c}: sum of sizes = {a}\n"
                                 f"The sizes > {b} were retained\n"
                                 f"The number of selected clusters for {sam2trim2c} at locus {loc2trim2c} "
                                 f"= {nb_selected}\n")
                os.remove("trim-select." + scriptExt)
                os.remove("tmp")

    if rmenu == "2d":
        os.chdir(current_dir)
        stat_2d = open(f"{current_dir}/outputs/Stats_option_2d.txt", "a")
        while True:
            loc2trim2d = in_loc2trim_2x()
            os.chdir(f"{current_dir}/loci/{loc2trim2d}")
            trim_left = in_trim_left()
            trim_right = in_trim_right()
            stat_2d.write(f"Locus {loc2trim2d} trimmed at {trim_left} bp (left) and {trim_right} bp (right)\n")
            while True:
                sam2trim2d = in_trim_sample2d()
                ts = in_ts()
                with open(sam2trim2d + "_singleEnd_orient.fas", "r") as filin, \
                        open("trim-select." + scriptExt, "w") as filout:
                    targets = [line for line in filin if "size" in line]
                    a = 0
                    for target in targets:
                        size = re.search('size=(.+?)$', target).group(1)
                        a = a + int(size)
                    b = int(a * float(ts1) + 1)
                    stat_2d.write(f"\tSum of sizes for {sam2trim2d} at locus {loc2trim2d} = {a}\n"
                                  f"\tAt threshold {ts} sizes > {b} for {sam2trim2d} were retained\n")
                    filout.writelines(f'' + start_log_redirect('./' + loc2trim2d + '.log') +
                                      f' vsearch --fastx_filter {sam2trim2d}_singleEnd_orient.fas --fastq_stripleft {trim_left} '
                                      f'  --fastq_stripright {trim_right} --fastaout tmp --minsize {b}\n'
                                      + localErrorOnStopCmd + "\n"
                                                              f' vsearch --derep_fulllength ./tmp --output {sam2trim2d}_singleEnd_select.fas --sizein '
                                                              f'--sizeout\n' + localErrorOnStopCmd + "\n"
                                                                                                     f'' + end_log_redirect(
                        './' + sam2trim2d + '.log'))
                subprocess.run([shellCmd, "./trim-select." + scriptExt])
                selected = open("./" + sam2trim2d + "_singleEnd_select.fas", "r")
                nb_selected = selected.read().count(">")
                stat_2d.write(f"\tNumber of selected clusters for {sam2trim2d} is: {nb_selected}\n\n")
                sys.stdout.write(f"\n{sam2trim2d}: sum of sizes = {a}\n"
                                 f"The sizes > {b} were retained\n"
                                 f"Number of selected clusters for {sam2trim2d} at locus {loc2trim2d} "
                                 f"= {nb_selected}\n")
                os.remove("trim-select." + scriptExt)
                os.remove("tmp")


def concat_3():
    """Clustering of all sample sequences by locus in a unique file (option 3)
    """
    os.system("cls" if winOS else "clear")
    global alloci, samples, loc2cat
    nb_samples = len(samples)
    alloci = loci1s + list(set(loci2s) - set(loci1s))
    if os.path.exists("outputs/Stats_option_3.txt"):
        os.remove("outputs/Stats_option_3.txt")
    while True:
        loc2cat = input(
            titleStyle + "\n----- MENU 3 - CONCATENATION OF ALL SAMPLE SELECTED SEQUENCES BY LOCUS FOR PHYLOGENETIC ANALYZES -----" + promptStyle + "\n"
                                                                                                                                  f"\nFor which LOCUS do you want to concatenate all sample sequences?\n" + normalStyle +
            f"among {alloci}?\n"
            f"OR 'end' 'home' 'exit': ")
        while loc2cat not in alloci and loc2cat not in ["end", "home", "exit"]:
            loc2cat = input(
                errorStyle + "\n--> WRONG INPUT: " + normalStyle + " locus name is not valid, please enter a valid name \n"
                                                                   f"among {alloci}: ")
        else:
            if loc2cat in ["end", "home"]:
                main()
            elif loc2cat == "exit":
                sys.stdout.write("\nSee the results of concatenation session\n"
                                 f"Results in ---> {current_dir}/outputs/Stats_option_3.txt\n" +
                                 successStyle + "\n\nCLUSTERING SESSION FOR PHYLOGENETIC PAIRED-END LOCI is complete\n" + normalStyle)
                quit_mbctools()
        stat_3 = open('./outputs/Stats_option_3.txt', 'a')
        os.chdir(f"./loci/{loc2cat}")
        files2cat = glob.glob('*_select.fas')
        if len(files2cat) == 0:
            print(
                errorStyle + f"\nSequences for locus {loc2cat} have not been filtered using an abundance threshold. Please run step 2 on all loci for which you want to run step 3\n" + normalStyle)
        else:
            with open(f"./{loc2cat}_allseq_select.fasta", "w") as out:
                for file in files2cat:
                    if os.path.exists(file):
                        with open(file, "r") as out2:
                            out.write(out2.read())
            tot = open("./" + loc2cat + "_allseq_select.fasta")
            nb_tot = tot.read().count(">")
            stat_3.writelines(f"The locus {loc2cat} has {nb_tot} sequences\n")
            sys.stdout.write(
                successStyle + f"\nLocus {loc2cat}: {nb_tot} sequences from {nb_samples} samples have been concatenated\n" + normalStyle +
                f"Results in ---> {current_dir}/{loc2cat}/{loc2cat}_allseq_select.fasta\n")
        os.chdir(current_dir)
        stat_3.close()


def prevent():
    """Forces the user to run the mandatory option 1 before any other option
    """
    global current_dir
    if os.path.isfile(f"{current_dir}/outputs/parameters_option_1.cfg") is False:
        sys.stdout.write(
            "\nYou have to run mandatory OPTION 1 " + warningStyle + "before" + normalStyle + " running this option \n")
        q = input(promptStyle + "\nDo you want to run OPTION 1? Reply 'yes', or anything else to exit: " + normalStyle)
        if q == 'yes':
            main_menu1()
            quit_mbctools()
        else:
            quit_mbctools()


def quit_mbctools():
    print(successStyle + "\n\nThanks for using mbctools!" + normalStyle)
    printHowToCite()
    quit()


def printHowToCite():
    print("\nPlease cite this software as follows:" +
          citationStyle + "\nmbctools: An open source tool to facilitate DNA sequence data analysis using VSEARCH in metabarcoding studies"
                          "\nChristian Barnabé, Guilhem Sempéré and Etienne Waleckx. https://github.com/GuilhemSempere/mbctools" + normalStyle + "\n\n")


def main_menu1():
    """Displays submenu 1
    """
    os.system("cls" if winOS else "clear")
    global rmenu
    rmenu = input(
        titleStyle + "\n----- MENU 1 - BASIC ANALYSIS - only option 1 is strictly mandatory -----" + normalStyle + "\n\n"
                                                                                                                   "1  -> NEW COMPLETE ANALYSIS (mandatory)\n"
                                                                                                                   "1a -> Re-analyze all loci, from the clustering step, modifying parameters\n"
                                                                                                                   "1b -> Re-analyze only one locus of paired-end amplicon (merged reads), modifying parameters\n"
                                                                                                                   "1c -> Re-analyze only one locus of single-end amplicon (R1 only), modifying parameters\n"
                                                                                                                   "1d -> Re-analyse a given sample, modifying parameters\n"
                                                                                                                   "1e -> Optional quality checking of fastq files (slow)\n\n"
                                                                                                                   "\n" + promptStyle + "Enter '1' '1a' '1b' '1c' '1d' '1e' to run analysis\n" + normalStyle +
        "OR 'end' 'home' 'exit': ")
    while rmenu not in ['1', '1a', '1b', '1c', '1d', '1e', 'end', 'home', 'exit']:
        rmenu = input(errorStyle + f"\n--> WRONG INPUT: " + normalStyle + " Please enter a CORRECT NAME of an option\n"
                                                                          f" among '1' '1a' '1b' '1c' '1d' '1e' 'home' 'exit': ")
    else:
        if rmenu in ["home", "end"]:
            main()
        elif rmenu == 'exit':
            quit_mbctools()
        elif rmenu == "1e":
            menu1e()
        elif rmenu == "1":
            menu1()
        elif rmenu == "1a":
            menu1a()
        elif rmenu == "1b":
            menu1b()
        elif rmenu == "1c":
            menu1c()
        elif rmenu == "1d":
            menu1d()
    return rmenu


def main_menu2():
    """Displays submenu 2
    """
    os.system("cls" if winOS else "clear")
    global rmenu
    rmenu = input(
        titleStyle + "\n----- MENU 2 - REMOVAL OF PRIMERS AND SELECTION OF MINIMUM SEQUENCE ABUNDANCES ACCORDING TO USER-DEFINED THRESHOLDS -----" + normalStyle + "\n\n"
                                                                                                                                            "2a -> Apply the SAME size threshold for ALL SAMPLES for the loci based on PAIRED-END reads "
                                                                                                                                            "(R1/R2 merged)\n"
                                                                                                                                            "\ti.e. you want to keep only sequences whose abundance (size)\n"
                                                                                                                                            "\tis greater than x% of the total number of sequences for a given sample.\n"
                                                                                                                                            "\tThis threshold of x% can be chosen for each locus.\n\n"
                                                                                                                                            "2b -> Apply the SAME size threshold for ALL SAMPLES for the loci based on SINGLE-END reads "
                                                                                                                                            "(R1 only)\n"
                                                                                                                                            "\tsame as option 1 but only using the R1 reads instead of merged ones.\n\n"
                                                                                                                                            "2c -> Apply a SPECIFIC size threshold  for a given sample, for the loci based on PAIRED-END reads "
                                                                                                                                            "(R1/R2 merged)\n"
                                                                                                                                            "\ti.e. you want to modulate the threshold of x% by locus but also by sample\n"
                                                                                                                                            "\twithin a particular locus.\n\n"
                                                                                                                                            "2d -> Apply a SPECIFIC size threshold  for a given sample, for the loci based on SINGLE-END reads "
                                                                                                                                            "(R1 only)\n"
                                                                                                                                            "\tsame as option 2c but only using the R1 sequences instead of merged ones.\n\n"
                                                                                                                                            "\n" + promptStyle + "Enter '2a' '2b' '2c' '2d' to run analysis\n" + normalStyle +
        "OR 'end' 'home' 'exit': ")
    while rmenu not in ['2a', '2b', '2c', '2d', 'end', 'home', 'exit']:
        rmenu = input(errorStyle + f"\n--> WRONG INPUT: " + normalStyle + " Please enter a CORRECT option number\n"
                                                                          f" among '2a' '2b' '2c' '2d' 'end' 'home' 'exit': ")
    else:
        if rmenu in ["home", "end"]:
            main()
        elif rmenu == 'exit':
            quit_mbctools()
        elif rmenu == "2a":
            menu2a()
        elif rmenu == "2b":
            menu2b()
        elif rmenu == "2c":
            menu2c()
        elif rmenu == "2d":
            menu2d()


def main_menu3():
    """Displays submenu 3
    """
    prevent()
    prev_param(None)
    concat_3()
    rerun()


def main_menu4():
    """Displays submenu 4
    """
    os.system("cls" if winOS else "clear")
    global rmenu
    rmenu = input(
        titleStyle + "\n----- MENU 4 - CONVERSION OF ANALYSIS RESULTS INTO metaXplor IMPORT FORMAT -----" + normalStyle + "\n\n"
                                                                                                                          "4a -> Generate sequence files\n"
                                                                                                                          "\tCompiles all sequences selected for all loci into a single fasta\n"
                                                                                                                          "\tOutputs a .tsv file indicating samples weights for each sequence\n\n"
                                                                                                                          "4b -> Generate assignment file\n"
                                                                                                                          "\tConverts blastn results (obtained from blasting above-mentioned fasta file) from 'Hit table (text)' (format #7) into metaXplor format\n\n"
                                                                                                                          "4c -> Builds metaXplor-format sample metadata file from provided tabulated file\n\n"
                                                                                                                          "4d -> Compresses all metaXplor files into a final, ready to import, zip archive\n\n"
                                                                                                                          "\n" + promptStyle + "Enter '4a' '4b' '4c' '4d' to generate required types of files\n" + normalStyle +
        "OR 'end' 'home' 'exit': ")
    while rmenu not in ['4a', '4b', '4c', '4d', 'end', 'home', 'exit']:
        rmenu = input(errorStyle + f"\n--> WRONG INPUT: " + normalStyle + " Please enter a CORRECT option number\n"
                                                                          f" among '4a' '4b' '4c' '4d' 'end' 'home' 'exit': ")
    else:
        if rmenu in ["home", "end"]:
            main()
        elif rmenu == 'exit':
            quit_mbctools()
        elif rmenu == "4a":
            menu4a()
        elif rmenu == "4b":
            menu4b()
        elif rmenu == "4c":
            menu4c()
        elif rmenu == "4d":
            menu4d(True)


def menu1():
    """Runs option 1
    """
    try:
        dir_fastq
    except NameError:
        in_dir_fastq()
    try:
        fastq_R1
    except NameError:
        in_fastq_R1()
    try:
        fastq_R2
    except NameError:
        in_fastq_R2()
    try:
        loci1
    except NameError:
        in_loci1()
    try:
        loci2
    except NameError:
        in_loci2()
    try:
        Samples
    except NameError:
        in_Samples()
    try:
        minsize
    except NameError:
        in_minsize()
    try:
        minseqlength
    except NameError:
        in_minseqlength()
    try:
        alpha
    except NameError:
        in_alpha()
    try:
        identity
    except NameError:
        in_identity()
    folders()
    param_1x()
    merging()
    fastq2fas()
    derep_1()
    cluster_1x()
    chimera_remove()
    runloc_merged()
    runloc_r1()
    orient_1x()
    sys.stdout.write("\n\n")
    runs_1x()
    stats_1x()
    rerun()
    if len(sys.argv) == 0:
        rerun()


def menu1a():
    """Runs option 1a
    """
    prevent()
    prev_param(None)
    in_minsize()
    in_minseqlength()
    in_alpha()
    in_identity()
    param_1x()
    cluster_1x()
    chimera_remove()
    runloc_merged()
    runloc_r1()
    orient_1x()
    runs_1x()
    stats_1x()
    rerun()


def menu1b():
    """Runs option 1b
    """
    prevent()
    prev_param(None)
    in_loc_sel_merged()
    in_minsize()
    in_minseqlength()
    in_alpha()
    in_identity()
    param_1x()
    cluster_1x()
    chimera_remove()
    runlocsel_merged()
    orient_1x()
    runs_1x()
    stats_1x()
    rerun()


def menu1c():
    """Runs option 1c
    """
    prevent()
    prev_param(None)
    in_loc_sel_r1()
    in_minsize()
    in_minseqlength()
    in_alpha()
    in_identity()
    param_1x()
    cluster_1x()
    runlocsel_r1()
    orient_1x()
    runs_1x()
    stats_1x()
    rerun()


def menu1d():
    """Runs option 1d
    """
    prevent()
    prev_param(None)
    in_sam_sel()
    in_minsize()
    in_minseqlength()
    in_alpha()
    in_identity()
    param_1x()
    cluster_1x()
    chimera_remove()
    runloc_one_sample_1d()
    orient_1x()
    runs_1x()
    stats_1x()
    rerun()


def menu1e():
    """Runs option 1e
    """
    prevent()
    prev_param(None)
    quality()
    runs_1x()
    rerun()


def menu2a():
    """Runs option 2a
    """
    prevent()
    prev_param(None)
    trim_2x()
    rerun()


def menu2b():
    """Runs option 2b
    """
    prevent()
    prev_param(None)
    trim_2x()
    rerun()


def menu2c():
    """Runs option 2c
    """
    prevent()
    prev_param(None)
    if os.path.exists("outputs/Stats_option_2c.txt"):
        os.remove("outputs/Stats_option_2c.txt")
    trim_2x()
    rerun()


def menu2d():
    """Runs option 2d
    """
    prevent()
    prev_param(None)
    if os.path.exists("outputs/Stats_option_2d.txt"):
        os.remove("outputs/Stats_option_2d.txt")
    trim_2x()
    rerun()


def menu3():
    """Runs option 3
    """
    prevent()
    prev_param(None)
    concat_3()
    rerun()


def menu4a():
    """Runs option 4a
    """
    prevent()
    prev_param(None)
    uniqueLoci = list(set(loci1s + loci2s))

    with open(metaXplorFasta, "w") as fastaFile, open(metaXplorSequenceComposition, "w") as seqCompositionFile:
        print("\nCompiling sequences for samples " + ", ".join(samples) + "\n")

        i = 0
        seqCompositionFile.write("qseqid")
        while i < len(samples):
            seqCompositionFile.write("\t" + samples[i])
            i += 1

        skippedLoci = []
        i = 0
        totalSeqCount = 0
        while i < len(uniqueLoci):
            try:
                file1 = open("loci/" + uniqueLoci[i] + fileSep + uniqueLoci[i] + '_allseq_select.fasta', 'r')
                lines = file1.readlines()

                j = 0
                seqCount = 0
                while j < len(lines):
                    line = lines[j].replace("\r\n", "").replace("\n", "")
                    if line.startswith(">sample="):
                        seqCount += 1
                        qseqid = re.sub(r';size=.*', '', line.replace(">sample=", ""))
                        fastaFile.write(">" + qseqid + "\n")
                        seqCompositionFile.write("\n" + qseqid)

                        k = 0
                        while k < len(samples):
                            seqCompositionFile.write(
                                "\t" + (line.split(";size=")[1] if qseqid.startswith(samples[k]) else "0"))
                            k += 1
                    else:
                        fastaFile.write(lines[j])
                    j += 1

                print("Processed locus " + uniqueLoci[i] + " with " + str(seqCount) + " sequences")
                totalSeqCount += seqCount
            except FileNotFoundError:
                print(warningStyle + "File not found: loci/" + uniqueLoci[i] + fileSep + uniqueLoci[
                    i] + '_allseq_select.fasta: skipping locus ' + uniqueLoci[i] + normalStyle)
                skippedLoci.append(uniqueLoci[i])
            i += 1

    if totalSeqCount == 0:
        print(
            errorStyle + "\nNo concatenated sequences found for any loci. Please run step 3 before retrying" + normalStyle)
        os.remove(metaXplorFasta)
        os.remove(metaXplorSequenceComposition)
    else:
        print(successStyle + "\n" + str(totalSeqCount) + " project sequences were compiled into .fasta and .tsv files")
        if len(skippedLoci) > 0:
            print(
                warningStyle + "Warning: not all sequences could be included because concatenation step (#3) was not run on some loci: " + ", ".join(
                    skippedLoci) + successStyle)
        print(
            "You may now run blastn on " + metaXplorFasta + ", download all results as 'Hit table (text)' (format #7), then come back and launch step 4b" + normalStyle)
    rerun()


def menu4b():
    """Runs option 4b
    """
    prevent()
    prev_param(None)

    blastTextHitTable = None
    print()
    while blastTextHitTable == None or blastTextHitTable.strip() == "" or (
            blastTextHitTable.strip() not in ['end', 'home', 'exit'] and not os.path.isfile(blastTextHitTable.strip())):
        if blastTextHitTable != None and blastTextHitTable.strip() != "":
            print(errorStyle + "\n--> UNEXISTING FILE: " + blastTextHitTable + normalStyle)
        blastTextHitTable = input(promptStyle + "Enter path to blastn hit-table (text format #7): " + normalStyle)
    if blastTextHitTable == "end":
        main_menu4()
    elif blastTextHitTable == "home":
        main()
    elif blastTextHitTable == "exit":
        quit_mbctools()

    with open(blastTextHitTable.strip(), "r") as infile:
        lines = re.sub('\s\s+', '\t', infile.read()).splitlines()

    if len(lines) == 0:
        print(errorStyle + "Provided hit-table file is empty!" + normalStyle)
        menu4b()

    database = None
    previousQseqId = None
    blastType = lines[0].split(" ")[1].strip()

    maxHits = None
    print()
    while maxHits == None or not maxHits.isnumeric() or int(maxHits) < 1:
        if maxHits != None:
            print(errorStyle + "\n--> WRONG INPUT: " + maxHits + normalStyle)
        maxHits = input(
            promptStyle + "Enter maximum number of retained hits per query." + normalStyle + " Default is 5: ")
        if maxHits == "":
            maxHits = "5"
    maxHits = int(maxHits)

    print()
    with open(metaXplorAssignments, "w") as outfile:
        i = 0
        nHitsForQseqId = 0
        while i < len(lines):
            if lines[i].startswith("#"):
                if database == None and "Database:" in lines[i]:
                    try:
                        database = lines[i].split(" ")[2]
                        if database == "nt":
                            database = "n"
                        elif database == "nr":
                            database = "p"
                        else:
                            raise Exception("Unsupported Database type: " + database)
                    except:
                        print(errorStyle + "Unable to parse accession type prefix in '" + lines[i] + "'" + normalStyle)
                        os.remove(metaXplorAssignments)
                        menu4b()
                elif previousQseqId == None and "Fields:" in lines[i]:
                    outfile.write(
                        lines[i].split(":")[1].strip().replace(", ", "\t").replace("query acc.ver", "qseqid").replace(
                            "subject acc.ver", "sseqid") + "\tassignment_method\tbest_hit\n")
            else:
                if database == None:
                    print(errorStyle + "Unable to determine accession type prefix" + normalStyle)
                    os.remove(metaXplorAssignments)
                    menu4b()

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

    print(successStyle + "File " + metaXplorAssignments + " was successfully written" + normalStyle)
    menu4d(False)
    rerun()


def menu4c():
    """Runs option 4c
    """
    prevent()
    prev_param(None)

    sampleMetadataFile = None
    print(
        "\nYou must now provide a tabulated metadata file for your samples. A header field named 'Sample' is expected for the column featuring sample names")
    print("Any field with 'date' in its header name will be considered to be the sample collection date")
    print("Collection location may be specified:")
    print(
        "\t- either as commma-separated decimal-format values in a single field named 'LatLon' (e.g. -17.7127, -67.9905)")
    print(
        "\t- or in separate columns named 'Latitude' and 'Longitude', in decimal format (e.g. -17.7127) or DMS format (e.g. 16°42'45.6\"S)")
    print("Any additional columns will remain named as provided")
    while sampleMetadataFile == None or sampleMetadataFile.strip() == "" or (
            sampleMetadataFile.strip() not in ['end', 'home', 'exit'] and not os.path.isfile(
            sampleMetadataFile.strip())):
        if sampleMetadataFile != None and sampleMetadataFile.strip() != "":
            print(errorStyle + "\n--> UNEXISTING FILE: " + sampleMetadataFile + normalStyle)
        sampleMetadataFile = input(promptStyle + "Enter path to tabulated sample metadata file: " + normalStyle)
    if sampleMetadataFile == "end":
        main_menu4()
    elif sampleMetadataFile == "home":
        main()
    elif sampleMetadataFile == "exit":
        quit_mbctools()

    print()
    with open(sampleMetadataFile.strip(), "r") as infile:
        lines = infile.read().splitlines()

    if len(lines) == 0:
        print(errorStyle + "Provided sample metadata file is empty!" + normalStyle)
        menu4c()

    headerCols = re.sub('\s\s+', '\t', lines[0]).split("\t")
    sampleIndex = None
    collectionDateIndex = None
    latitudeIndex = None
    longitudeIndex = None
    latLonIndex = None

    i = 0
    while i < len(headerCols):
        headerCol = headerCols[i].lower()
        # print(str(i) + " : " + headerCol)
        if "sample" in headerCol:
            if sampleIndex != None:
                print(errorStyle + "Ambiguity identifying sample name column between '" + headerCols[
                    sampleIndex] + "'' and '" + headerCols[i] + "'" + normalStyle)
                menu4c()
            sampleIndex = i
        elif "date" in headerCol:
            if collectionDateIndex != None:
                print(errorStyle + "Ambiguity identifying collection date column between '" + headerCols[
                    collectionDateIndex] + "'' and '" + headerCols[i] + "'" + normalStyle)
                menu4c()
            collectionDateIndex = i
        elif headerCol.startswith("lat"):
            if "lon" in headerCol:
                if latLonIndex != None:
                    print(errorStyle + "Ambiguity identifying LatLon column between '" + headerCols[
                        latLonIndex] + "'' and '" + headerCols[i] + "'" + normalStyle)
                    menu4c()
                latLonIndex = i
            else:
                if latitudeIndex != None:
                    print(errorStyle + "Ambiguity identifying latitude column between '" + headerCols[
                        latitudeIndex] + "'' and '" + headerCols[i] + "'" + normalStyle)
                    menu4c()
                latitudeIndex = i
        elif headerCol.startswith("lon"):
            if "lat" in headerCol:
                if latLonIndex != None:
                    print(errorStyle + "Ambiguity identifying LatLon column between '" + headerCols[
                        latLonIndex] + "'' and '" + headerCols[i] + "'" + normalStyle)
                    menu4c()
                latLonIndex = i
            else:
                if longitudeIndex != None:
                    print(errorStyle + "Ambiguity identifying longitude column between '" + headerCols[
                        longitudeIndex] + "'' and '" + headerCols[i] + "'" + normalStyle)
                    menu4c()
                longitudeIndex = i
        i += 1

    if sampleIndex == None:
        print(
            errorStyle + "Unable to identify sample name column! Please read instructions carefully and submit a corrected file" + normalStyle)
        menu4c()

    if latLonIndex == None and (not None not in [latitudeIndex, longitudeIndex]):
        print(
            warningStyle + "Unable to identify latitude and/or longitude column(s)! Generated file will contain empty values for this field" + normalStyle)
    if collectionDateIndex == None:
        print(
            warningStyle + "Unable to identify collection date column! Generated file will contain empty values for this field" + normalStyle)

    specialIndexes = [sampleIndex, latitudeIndex, longitudeIndex, latLonIndex, collectionDateIndex]
    with open(metaXplorSamples, "w") as outfile:
        j = 0

        outfile.write("sample_name\tlat_lon\tcollection_date")
        while j < len(headerCols):
            if j not in specialIndexes:
                outfile.write("\t" + headerCols[j])
            j += 1
        outfile.write("\n")
        i = 1
        samplesToProcess = samples.copy()
        while i < len(lines):
            splitLine = re.sub('  +', '\t', lines[i]).split("\t")

            if splitLine[sampleIndex] not in samples:
                print(warningStyle + "Skipping unknown sample: " + splitLine[sampleIndex] + normalStyle)
            else:
                while len(splitLine) < len(headerCols):
                    splitLine.append("")

                outfile.write(splitLine[sampleIndex] + "\t" + determineLatLon(splitLine, latLonIndex, latitudeIndex,
                                                                              longitudeIndex, sampleIndex) + "\t")
                collDate = splitLine[collectionDateIndex].strip() if collectionDateIndex != None else None
                if collDate != None:
                    try:
                        outfile.write(str(dateutil.parser.parse(collDate)).replace(" 00:00:00", ""))
                    except Exception as e:
                        print(errorStyle + "Unable to parse date at line " + str(i) + ": " + collDate + normalStyle)
                        menu4c()

                j = 0
                while j < len(splitLine):
                    if j not in specialIndexes:
                        outfile.write("\t" + splitLine[j])
                    j += 1
                outfile.write("\n")
                samplesToProcess.remove(splitLine[sampleIndex])
            i += 1

    if len(samplesToProcess) > 0:
        print(errorStyle + "Provided file lacks lines for the following sample(s): " + ", ".join(
            samplesToProcess) + normalStyle)
        os.remove(metaXplorSamples)
    else:
        print(successStyle + "File " + metaXplorSamples + " was successfully written" + normalStyle)
    menu4d(False)
    rerun()


def menu4d(invokedByUser):
    """Runs option 4d
    """
    if not os.path.isfile(metaXplorFasta) or os.path.getsize(metaXplorFasta) == 0:
        if invokedByUser:
            print(errorStyle + "\nFile " + metaXplorFasta + " is missing or empty. Please run step 4a" + normalStyle)
        rerun()
    if not os.path.isfile(metaXplorSequenceComposition) or os.path.getsize(metaXplorSequenceComposition) == 0:
        if invokedByUser:
            print(
                errorStyle + "\nFile " + metaXplorSequenceComposition + " is missing or empty. Please run step 4a" + normalStyle)
        rerun()
    if not os.path.isfile(metaXplorAssignments) or os.path.getsize(metaXplorAssignments) == 0:
        if invokedByUser:
            print(
                errorStyle + "\nFile " + metaXplorAssignments + " is missing or empty. Please run step 4b" + normalStyle)
        rerun()
    if not os.path.isfile(metaXplorSamples) or os.path.getsize(metaXplorSamples) == 0:
        if invokedByUser:
            print(errorStyle + "\nFile " + metaXplorSamples + " is missing or empty. Please run step 4c" + normalStyle)
        rerun()

    global zipNow
    zipNow = input(
        "All metaXplor files seem to be ready. Zip them now to create the final import file? Please enter 'yes' or 'no': ")
    while zipNow not in ["yes", "no"]:
        zipNow = input(errorStyle + "\n--> WRONG INPUT: " + normalStyle + " Please enter 'yes' or 'no': ")

    if zipNow == "yes":
        b = io.BytesIO()
        zf = zipfile.ZipFile(b, mode='w')
        zf.write(metaXplorSamples, metaXplorSamples)
        zf.write(metaXplorAssignments, metaXplorAssignments)
        zf.write(metaXplorFasta, metaXplorFasta)
        zf.write(metaXplorSequenceComposition, metaXplorSequenceComposition)
        zf.close()
        zipFileName = 'mbctools_metaXplor_export_' + datetime.datetime.strptime("12/10/2020", "%d/%m/%Y").strftime(
            '%Y%m%d') + '.zip'
        open(zipFileName, 'wb').write(b.getbuffer())
        print(successStyle + "\nmetaXplor import archive was successfully created as " + zipFileName + normalStyle)
    rerun()


def determineLatLon(cellArray, latLonIndex, latitudeIndex, longitudeIndex, sampleIndex):
    gotSeparateLatAndLong = None not in [latitudeIndex, longitudeIndex]
    if latLonIndex != None:
        splitCoords = cellArray[latLonIndex].replace(";", ",").split(",")
        try:
            return str(round(float(splitCoords[0]), 6)) + ", " + str(round(float(splitCoords[1]), 6))
        except:
            pass

    latLon = ""
    if gotSeparateLatAndLong:
        try:
            latLon += str(round(float(cellArray[latitudeIndex]), 6))
        except:
            try:
                latLon += str(dmsToDecimal(cellArray[latitudeIndex]))
            except:
                latLon = ""
        if latLon != "":
            try:
                latLon += ", " + str(round(float(cellArray[longitudeIndex]), 6))
            except:
                try:
                    latLon += ", " + str(dmsToDecimal(cellArray[longitudeIndex]))
                except:
                    latLon = ""

    if latLon == "":
        msg = "Sample " + cellArray[sampleIndex] + ":"
        if latLonIndex != None:
            msg += " Unable to parse LatLon field" + (
                "" if gotSeparateLatAndLong else (" '" + cellArray[latLonIndex] + "'")) + "."
        if gotSeparateLatAndLong:
            msg += " Unable to parse latitude / longitude fields."
        print(warningStyle + msg + normalStyle)
    return latLon


def dmsToDecimal(dmsString):
    deg, minutes, seconds, direction = re.split('[°\'"]', dmsString.replace("''", "\"").replace(" ", ""))
    return round((float(deg) + float(minutes) / 60 + float(seconds) / (60 * 60)) * (
        -1 if direction.upper() in ['W', 'S'] else 1), 6)


def rerun():
    """Prompts the user to continue using mbctools or not
    """
    global next_run
    next_run = input("Do you want to continue with mbctools? Please enter 'yes' or 'no': ")
    while next_run not in ["yes", "no", ""]:
        next_run = input(errorStyle + "\n--> WRONG INPUT: " + normalStyle + " Please enter 'yes' or 'no'\n"
                                                                            "Default is 'yes': ")
    if next_run in ["yes", ""]:
        main()
    if next_run == "no":
        quit_mbctools()


def main():
    """Displays the main menu
    """

    if len(sys.argv) > 1:
        if os.path.isfile(sys.argv[1]):
            prev_param(sys.argv[1])
            global rmenu
            rmenu = "1"
            menu1()
            exit(0)
        else:
            print(errorStyle + "\nUnexisting configuration file: " + sys.argv[1] + "\n" + normalStyle)
            exit(1)

    os.system("cls" if winOS else "clear")

    global menu
    sys.stdout.write(titleStyle + "-------------------- mbctools - MAIN MENU --------------------" + normalStyle + "\n")
    printHowToCite()
    sys.stdout.write("Validating without typing anything applies the default value, if any\n"
                     "Entering 'end' returns to the program upper level, if any\n"
                     "Entering 'home' returns to this main menu\n"
                     "Entering 'exit' leaves the program\n\n")
    menu = input("\n1 -> BASIC ANALYZES\n\n"
                 "2 -> REMOVAL OF PRIMERS AND SELECTION OF MINIMUM SEQUENCE ABUNDANCES ACCORDING TO USER-DEFINED THRESHOLDS\n\n"
                 "3 -> CONCATENATION OF ALL SAMPLE SELECTED SEQUENCES BY LOCUS FOR PHYLOGENETIC ANALYZES\n\n"
                 "4 -> CONVERSION OF ANALYSIS RESULTS INTO metaXplor IMPORT FORMAT\n\n"
                 "\n" + promptStyle + "Enter '1', '2', '3', '4'\n" + normalStyle +
                 "OR 'exit' to quit the program: ")
    while menu not in ['1', '2', '3', '4', 'exit']:
        menu = input(
            errorStyle + "\n--> WRONG INPUT: " + normalStyle + " Please enter a CORRECT NAME of an option among '1', '2', '3', '4'\n"
                                                               "OR 'exit' to quit the program: ")
    else:
        if menu == 'exit':
            quit_mbctools()
        if menu == '1':
            main_menu1()
        if menu == '2':
            main_menu2()
        if menu == '3':
            main_menu3()
        if menu == '4':
            main_menu4()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n")
