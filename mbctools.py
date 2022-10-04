#!/usr/bin/python3

"""This Python program is a toolkit to make the use of VSEARCH easier and interactive, to analyze
metabarcoding NGS data in the best conditions. It proposes the following MAIN MENU:

0  -> Initial and optional quality checking of fastq files (slow - be patient!)
1  -> NEW COMPLETE ANALYSIS
1a -> Re-analyze all loci, from the clustering step, modifying parameters
1b -> Re-analyze only one locus of paired-end amplicon (merged reads), modifying parameters
1c -> Re-analyze only one locus of single-end amplicon (R1 only), modifying parameters
1d -> Re-analyse only one sample, modifying parameters
2  -> SELECTION OF MINIMUM SIZES ACCORDING TO USER-DEFINED THRESHOLDS
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
import numpy as np


# BEFORE ANYTHING ELSE: MAKE SURE VSEARCH IS INSTALLED ################################################################
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
date = datetime.date.today()
current_dir = os.getcwd()
liste = str(np.arange(0, 1, 0.001))

global loci1s, loci2s, loci1_user, loci2_user, dir_fastq, fastqr1_user, fastqr1s, fastqr2_user, fastqr2s, sample_user, \
    samples, alpha, identity, loc_sel1, loc_sel2, rmenu, sam_sel, minsize_user, minseqlength, loc2trim, trim_left, \
    trim_right, ts, ts1, orient2c, orient2d, alloci, next_run, menu, ts2, loc2cat, loc2trim2a, loc2trim2b, \
    loc2trim2c, loc2trim2d


# FUNCTIONS FOR MULTIPLATFORM SUPPORT #################################################################################
def start_log_redirect(filepath):
    if not winOS:
        return "exec 3>&1 4>&2 >" + filepath + " 2>&1\n"
    else:
        return "&{\n"


def end_log_redirect(filepath):
    if not winOS:
        return "exec 1>&3 2>&4\n"
    else:
        return "} 2> " + filepath + "\n"


def main_stream_message(message):
    if not winOS:
        return "printf \"" + message + "\" >&3\n"
    else:
        return "Write-Host -NoNewline \"" + message + "\"\n"


# FOLDER CREATION #####################################################################################################
def folders():
    global loci1s, loci2s, loci1_user, loci2_user
    sys.stdout.write("")
    folder = "scripts"
    path = os.path.join(current_dir, folder)
    if Path(path).is_dir():
        sys.stdout.write(f"\nThe folder {folder} already exists and will be used in the current analysis")
    else:
        os.mkdir(path)
    sys.stdout.write("")
    folder = "outputs"
    path = os.path.join(current_dir, folder)
    if Path(path).is_dir():
        sys.stdout.write(f"\nThe folder {folder} already exists and will be used in the current analysis")
    else:
        os.mkdir(path)
    sys.stdout.write("")
    folder = "refs"
    path = os.path.join(current_dir, folder)
    if Path(path).is_dir():
        sys.stdout.write(f"\nThe folder {folder} already exists and will be used in the current analysis")
    else:
        os.mkdir(path)
    sys.stdout.write("")
    folder = "tmp_files"
    path = os.path.join(current_dir, folder)
    if Path(path).is_dir():
        sys.stdout.write(f"\nThe folder {folder} already exists and will be used in the current analysis")
    else:
        os.mkdir(path)
    sys.stdout.write("")  # loci folders
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


# INPUTS AND VARIABLES ################################################################################################
def in_dir_fastq():  # menu1 1
    global dir_fastq
    dir_fastq = input("\n\nEnter the FULL PATH to the folder where fastq files are located"
                      f"\nDefault is here:\n"
                      f"{current_dir}/fastq\n"
                      f"enter path: ")
    while Path(dir_fastq).is_dir() is False and dir_fastq not in ['', 'end', 'home', 'exit']:
        dir_fastq = input("\n-->ERROR: path not valid, enter a valid path\n"
                          "default is {current_dir}/fastq\n"
                          "OR 'end' 'home' 'exit': ")
    else:
        if dir_fastq == "":
            dir_fastq = f"{current_dir}/fastq"
        elif dir_fastq == "end":
            main_menu1()
        elif dir_fastq == "home":
            main()
        elif dir_fastq == "exit":
            sys.stdout.write("\n\tBye!...\n\n")
            quit()


def in_fastqr1_user():
    global fastqr1_user, fastqr1s
    fastqr1_user = input("\nEnter R1 fastq file name\n"
                         "default = fastqR1.txt: ")
    while os.path.isfile(fastqr1_user) is False and fastqr1_user not in ["", "end", "home", "exit"]:
        fastqr1_user = input("\n-->ERROR: file name is not valid, please enter a valid file name\n"
                             "default = fastqR1.txt\n"
                             "OR 'end' 'home' 'exit': ")
    else:
        if fastqr1_user == "":
            fastqr1_user = 'fastqR1.txt'
        elif fastqr1_user == "end":
            main_menu1()
        elif fastqr1_user == "home":
            main()
        elif fastqr1_user == "exit":
            sys.stdout.write("\n\tBye!...\n\n")
            quit()
    with open(fastqr1_user, "r") as out1:
        fastqr1s = out1.read().splitlines()


def in_fastqr2_user():
    global fastqr2_user, fastqr2s
    fastqr2_user = input("\nEnter R2 fastq file name\n"
                         "default = fastqR2.txt: ")
    while os.path.isfile(fastqr2_user) is False and fastqr2_user not in ["", "end", "home", "exit"]:
        fastqr2_user = input("\n-->ERROR: file name is not valid, please enter a valid file name\n"
                             "default = fastqR2.txt\n"
                             "OR 'end' 'home' 'exit': ")
    else:
        if fastqr2_user == '':
            fastqr2_user = 'fastqR2.txt'
        elif fastqr2_user == "end":
            main_menu1()
        elif fastqr2_user == "home":
            main()
        elif fastqr2_user == "exit":
            sys.stdout.write("\n\tBye!...\n\n")
            quit()
    with open(fastqr2_user, "r") as out2:
        fastqr2s = out2.read().splitlines()


def in_loci1_user():  # menu1 4
    global loci1_user, loci1s
    loci1_user = input("\nEnter the name of the file containing loci based on paired-end reads\n"
                       "default = locus1.txt: ")
    while os.path.isfile(loci1_user) is False and loci1_user not in ["", "end", "home", "exit"]:
        loci1_user = input("\n-->ERROR: file name is not valid, please enter a valid file name\n"
                           "default = locus1.txt\n"
                           "OR 'end' 'home' 'exit': ")
    else:
        if loci1_user == '':
            loci1_user = 'locus1.txt'
        elif loci1_user == "end":
            main_menu1()
        elif loci1_user == "home":
            main()
        elif loci1_user == "exit":
            sys.stdout.write("\n\tBye!...\n\n")
            quit()
    with open(loci1_user, "r") as out:
        loci1s = out.read().splitlines()
    return loci1s


def in_loci2_user():
    global loci2_user, loci2s
    loci2_user = input("\nEnter the name of the file containing loci based on single-end reads only (R1)\n"
                       "default = locus2.txt: ")
    while os.path.isfile(loci2_user) is False and loci2_user not in ["", "end", "home", "exit"]:
        loci2_user = input("\n-->ERROR: file name is not valid, please enter a valid file name\n"
                           "default = locus2.txt: "
                           "OR 'end' 'home' 'exit': ")
    else:
        if loci2_user == "":
            loci2_user = 'locus2.txt'
        elif loci2_user == "end":
            main_menu1()
        elif loci2_user == "home":
            main()
        elif loci2_user == "exit":
            sys.stdout.write("\n\tBye!...\n\n")
            quit()
    with open(loci2_user, "r") as out:
        loci2s = out.read().splitlines()
    return loci2s


def in_sample_user():
    global sample_user, samples
    sample_user = input("\nEnter the name of the file containing sample names\n"
                        "default = samples.txt: ")
    while os.path.isfile(sample_user) is False and sample_user not in ["", "end", "home", "exit"]:
        sample_user = input("\n-->ERROR: file name is not valid, please enter a valid file name\n"
                            "default = samples.txt\n"
                            "OR 'end' 'home' 'exit': ")
    else:
        if sample_user == "":
            sample_user = 'samples.txt'
        elif sample_user == "end":
            main_menu1()
        elif sample_user == "home":
            main()
        elif sample_user == "exit":
            sys.stdout.write("\n\tBye!...\n\n")
            quit()
    with open(sample_user, "r") as out5:
        samples = out5.read().splitlines()


def in_minsize_user():  # menu1 7# menu1a 3
    global minsize_user
    minsize_user = input("\nEnter the minsize option value for clusters,\n"
                         "i.e. the minimum sequence abundance of the retained clusters\n"
                         "default = 8: ")
    while minsize_user.isnumeric() is False and minsize_user not in ["", "end", "home", "exit"]:
        minsize_user = input("\n-->ERROR: minsize option must be an integer\n"
                             "enter an integer (e.g. 2, 10, 50...)\n"
                             "default = 8\n"
                             "OR 'end' 'home' 'exit': ")
    else:
        if minsize_user == '':
            minsize_user = '8'
        elif minsize_user == "end":
            main_menu1()
        elif minsize_user == "home":
            main()
        elif minsize_user == "exit":
            sys.stdout.write("\n\tBye!...\n\n")
            quit()


def in_minseqlength():  # menu1 8# menu1a 4
    global minseqlength
    minseqlength = input("\nEnter the minimum length of sequences to keep for any locus\n"
                         "default = 100: ")
    while minseqlength.isnumeric() is False and minseqlength not in ["", "end", "home", "exit"]:
        minseqlength = input("\n-->ERROR: minimum length must be an integer\n"
                             "enter an integer (e.g. 100, 150, 180...)\n"
                             "default = 100\n"
                             "OR 'end' 'home' 'exit': ")
    else:
        if minseqlength == '':
            minseqlength = '100'
        elif minseqlength == "end":
            main_menu1()
        elif minseqlength == "home":
            main()
        elif minseqlength == "exit":
            sys.stdout.write("\n\tBye!...\n\n")
            quit()


def in_alpha():  # menu1 9# menu1a 5
    global alpha
    alpha = input("\nEnter alpha parameter for the clustering\n"
                  "default = 2: ")
    while alpha.isnumeric() is False and alpha not in ["", "end", "home", "exit"]:
        alpha = input("\n-->ERROR: alpha parameter must be an integer\n"
                      "enter an integer (e.g. 1, 2, 3...)\n"
                      "default = 2\n"
                      "OR 'end' 'home' 'exit': ")
    else:
        if alpha == '':
            alpha = '2'
        elif alpha == "end":
            main_menu1()
        elif alpha == "home":
            main()
        elif alpha == "exit":
            sys.stdout.write("\n\tBye!...\n\n")
            quit()


def in_identity():  # menu1 10# menu1a 6
    global identity, liste
    identity = input("\nEnter identity parameter to BLAST the clusters against references\n"
                     "i.e. the identity percentage, enter an integer from 0 to 100\n"
                     "default = 70: ")
    while identity not in ["end", "home", "exit", ""] and identity.isnumeric() is False:
        identity = input("\n-->ERROR: identity parameter must be an integer from 0 to 100 \n"
                         "default = 70\n"
                         "OR 'end' 'home' 'exit': ")
    else:
        if identity == '':
            identity = '0.7'
        if identity == "end":
            main_menu1()
        if identity == "home":
            main()
        if identity.isnumeric() and int(identity) < 100:
            identity = int(identity) / 100            
        if identity == "exit":
            sys.stdout.write("\n\tBye!...\n\n")
            quit()


def in_loc_sel_merged():
    global loc_sel1
    loc_sel1 = input("\nEnter the name of the locus analysed by paired-end reads you want to rerun\n"
                     f"among {loci1s}"
                     f"no default: ")
    while loc_sel1 not in loci1s and loc_sel1 not in ["end", "home", "exit"]:
        loc_sel1 = input("\n-->ERROR: locus name is not valid, please enter valid name of locus\n"
                         f"among {loci1s}\n"
                         "OR 'end' 'home' 'exit': ")
    else:
        if loc_sel1 == "end":
            main_menu2()
        elif loc_sel1 == "home":
            main()
        elif loc_sel1 == "exit":
            sys.stdout.write("\n\tBye!...\n\n")
            quit()


def in_loc_sel_r1():
    global loc_sel2, rmenu
    loc_sel2 = input("\nEnter the name of the locus analysed by only single-end (R1) reads\n"
                     f"among {loci2s} you want to rerun\n"
                     f"no default: ")
    while loc_sel2 not in loci2s and loc_sel2 not in ["end", "home", "exit"]:
        loc_sel2 = input("\n-->ERROR: locus name is not valid, please enter valid name of locus\n"
                         f"among {loci2s}\n"
                         "OR 'end' 'home' 'exit': ")
    else:
        if loc_sel2 == "end":
            main_menu2()
        elif loc_sel2 == "home":
            main()
        elif loc_sel2 == "exit":
            sys.stdout.write("\n\tBye!...\n\n")
            quit()


def in_sam_sel():
    global sam_sel, samples
    sam_sel = input(f"\nEnter the sample name you want to rerun\n"
                    f"among {samples}\n"
                    f"no default: ")
    while sam_sel not in samples and sam_sel not in ["end", "home", "exit"]:
        sam_sel = input("\n-->ERROR: sample name is not valid, please enter a valid sample name\n"
                        f"among {samples}\n"
                        f"no default\n"
                        "OR 'end' 'home' 'exit': ")
    else:
        if sam_sel == "end" and rmenu == "1d":
            main_menu1()
        elif sam_sel == "home":
            main()
        elif sam_sel == "exit:":
            sys.stdout.write("\n\tBye!...\n\n")
            quit()


def in_loc2trim_2x():
    global loc2trim2a, loc2trim2b, loc2trim2c, loc2trim2d, loci1s, loci2s
    if rmenu == "2a":
        loc2trim2a = input("\nEnter a (new) LOCUS name based on paired-end mergeable reads you want to analyze among:\n"
                           f"{loci1s}\n"
                           "OR 'end' 'home' 'exit': ")

        while loc2trim2a not in loci1s and loc2trim2a not in ["end", "home", "exit"]:
            loc2trim2a = input("\n-->ERROR: locus name is not valid, please enter a valid name\n"
                               f"among {loci1s}\n"
                               "OR 'end' 'home' 'exit': ")
        else:
            if loc2trim2a == "end":
                main_menu2()
            if loc2trim2a == "home":
                main()
            elif loc2trim2a == "exit":
                sys.stdout.write("\nSelection of minimum sizes according to user-defined thresholds is achieved\n"
                                 "Statistical summary for paired-end based loci:\n"
                                 f"Results in ---> {current_dir}/outputs/Stats_option_2a.txt\n")
                sys.stdout.write("\n\tBye!...\n\n")
                quit()
        return loc2trim2a

    if rmenu == "2b":
        loc2trim2b = input("\nEnter a (new) LOCUS name based on single-end R1 reads you want to analyze among:\n"
                           f"{loci2s}\n"
                           "OR 'end' 'home' 'exit': ")

        while loc2trim2b not in loci1s and loc2trim2b not in ["end", "home", "exit"]:
            loc2trim2b = input("\n-->ERROR: locus name is not valid, please enter a valid name\n"
                               f"among {loci1s}\n"
                               "OR 'end' 'home' 'exit': ")
        else:
            if loc2trim2b == "end":
                main_menu2()
            if loc2trim2b == "home":
                main()
            elif loc2trim2b == "exit":
                sys.stdout.write("\nSelection of minimum sizes according to user-defined thresholds is achieved\n"
                                 "Statistical summary for paired-end based loci:\n"
                                 f"Results in ---> {current_dir}/outputs/Stats_option_2a.txt\n")
                sys.stdout.write("\n\tBye!...\n\n")
                quit()
        return loc2trim2b

    if rmenu == "2c":
        loc2trim2c = input("\nEnter a (new) LOCUS name based on paired-end mergeable reads you want to analyze among:\n"
                           f"{loci1s}\n"
                           "OR 'end' 'home' 'exit': ")

        while loc2trim2c not in loci1s and loc2trim2c not in ["end", "home", "exit"]:
            loc2trim2c = input("\n-->ERROR: locus name is not valid, please enter a valid name\n"
                               f"among {loci1s}\n"
                               "OR 'end' 'home' 'exit': ")
        else:
            if loc2trim2c == "end":
                main_menu2()
            if loc2trim2c == "home":
                main()
            elif loc2trim2c == "exit":
                sys.stdout.write("\nSelection of minimum sizes according to user-defined thresholds is achieved\n"
                                 "Statistical summary for paired-end based loci:\n"
                                 f"Results in ---> {current_dir}/outputs/Stats_option_2a.txt\n")
                sys.stdout.write("\n\tBye!...\n\n")
                quit()
        return loc2trim2c

    if rmenu == "2d":
        loc2trim2d = input("\nEnter a (new) LOCUS name based on paired-end mergeable reads you want to analyze among:\n"
                           f"{loci2s}\n"
                           "OR 'end' 'home' 'exit': ")

        while loc2trim2d not in loci2s and loc2trim2d not in ["end", "home", "exit"]:
            loc2trim2d = input("\n-->ERROR: locus name is not valid, please enter a valid name\n"
                               f"among {loci1s}\n"
                               "OR 'end' 'home' 'exit': ")
        else:
            if loc2trim2d == "end":
                main_menu2()
            if loc2trim2d == "home":
                main()
            elif loc2trim2d == "exit":
                sys.stdout.write("\nSelection of minimum sizes according to user-defined thresholds is achieved\n"
                                 "Statistical summary for paired-end based loci:\n"
                                 f"Results in ---> {current_dir}/outputs/Stats_option_2a.txt\n")
                sys.stdout.write("\n\tBye!...\n\n")
                quit()
        return loc2trim2d


def in_trim_left():
    global trim_left, loc2trim2a, loc2trim2b,  loc2trim2c
    if rmenu == "2a":
        trim_left = input(f"\nEnter the number of bp of the left primer for {loc2trim2a}? (e.g. 20): ")
        while trim_left.isnumeric() is False and trim_left not in ["end", "home", "exit"]:
            trim_left = input("\n-->ERROR: enter an integer corresponding to the length of the left primer (e.g. 20)\n"
                              "OR 'end' 'home' 'exit': ")
        else:
            if trim_left == "end":
                main_menu2()
            elif trim_left == "home":
                main()
            elif trim_left == "exit":
                sys.stdout.write("\n\tBye!...\n\n")
                quit()

    if rmenu == "2b":
        trim_left = input(f"\nEnter the number of bp of the left primer for {loc2trim2b}? (e.g. 20): ")
        while trim_left.isnumeric() is False and trim_left not in ["end", "home", "exit"]:
            trim_left = input("\n-->ERROR: enter an integer corresponding to the length of the left primer (e.g. 20)\n"
                              "OR 'end' 'home' 'exit': ")
        else:
            if trim_left == "end":
                main_menu2()
            elif trim_left == "home":
                main()
            elif trim_left == "exit":
                sys.stdout.write("\n\tBye!...\n\n")
                quit()

    if rmenu == "2c":
        trim_left = input(f"\nEnter the number of bp of the left primer for {loc2trim2c}? (e.g. 20): ")
        while trim_left.isnumeric() is False and trim_left not in ["end", "home", "exit"]:
            trim_left = input("\n-->ERROR: enter an integer corresponding to the length of the left primer (e.g. 20)\n"
                              "OR 'end' 'home' 'exit': ")
        else:
            if trim_left == "end":
                main_menu2()
            elif trim_left == "home":
                main()
            elif trim_left == "exit":
                sys.stdout.write("\n\tBye!...\n\n")
                quit()

    if rmenu == "2d":
        trim_left = input(f"\nEnter the number of bp of the left primer for {loc2trim2d}? (e.g. 20): ")
        while trim_left.isnumeric() is False and trim_left not in ["end", "home", "exit"]:
            trim_left = input("\n-->ERROR: enter an integer corresponding to the length of the left primer (e.g. 20)\n"
                              "OR 'end' 'home' 'exit': ")
        else:
            if trim_left == "end":
                main_menu2()
            elif trim_left == "home":
                main()
            elif trim_left == "exit":
                sys.stdout.write("\n\tBye!...\n\n")
                quit()
    return trim_left


def in_trim_right():
    global loc2trim2a, loc2trim2b, loc2trim2c, loc2trim2d, trim_right
    if rmenu == "2a":
        trim_right = input(f"\nEnter the number of bp of the right primer for {loc2trim2a}? (e.g. 22): ")
        while trim_right.isnumeric() is False and trim_right not in ["end", "home", "exit"]:
            trim_right = input("\n-->ERROR: enter an integer corresponding to the length of the "
                               "right primer (e.g. 22)\n"
                               "OR 'end' 'home' 'exit': ")
        else:
            if trim_right == "end":
                main_menu2()
            elif trim_right == "home":
                main()
            elif trim_right == "exit":
                sys.stdout.write("\n\tBye!...\n\n")
                quit()

    if rmenu == "2b":
        trim_right = input(f"\nEnter the number of bp of the right primer for {loc2trim2b}? (e.g. 22): ")
        while trim_right.isnumeric() is False and trim_right not in ["end", "home", "exit"]:
            trim_right = input("\n-->ERROR: enter an integer corresponding to the length of the "
                               "right primer (e.g. 22)\n"
                               "!! May be 0 with single-end based loci"
                               "OR 'end' 'home' 'exit': ")
        else:
            if trim_right == "end":
                main_menu2()
            elif trim_right == "home":
                main()
            elif trim_right == "exit":
                sys.stdout.write("\n\tBye!...\n\n")
                quit()

    if rmenu == "2c":
        trim_right = input(f"\nEnter the number of bp of the right primer for {loc2trim2c}? (e.g. 22): ")
        while trim_right.isnumeric() is False and trim_right not in ["end", "home", "exit"]:
            trim_right = input("\n-->ERROR: enter an integer corresponding to the length of the "
                               "right primer (e.g. 22)\n"
                               "OR 'end' 'home' 'exit': ")
        else:
            if trim_right == "end":
                main_menu2()
            elif trim_right == "home":
                main()
            elif trim_right == "exit":
                sys.stdout.write("\n\tBye!...\n\n")
                quit()

    if rmenu == "2d":
        trim_right = input(f"\nEnter the number of bp of the right primer for {loc2trim2d}? (e.g. 22): ")
        while trim_right.isnumeric() is False and trim_right not in ["end", "home", "exit"]:
            trim_right = input("\n-->ERROR: enter an integer corresponding to the length of the "
                               "right primer (e.g. 22)\n"
                               "!! May be 0 with single-end based loci"
                               "OR 'end' 'home' 'exit': ")
        else:
            if trim_right == "end":
                main_menu2()
            elif trim_right == "home":
                main()
            elif trim_right == "exit":
                sys.stdout.write("\n\tBye!...\n\n")
                quit()
    return trim_right


def in_ts():
    global ts, ts1, liste, loc2trim2a, loc2trim2b, loc2trim2c, loc2trim2d
    if rmenu == "2a":
        ts = input(f"\nEnter the THRESHOLD do you want to use for this locus {loc2trim2a}?\n"
                   f"Example: if you want to keep only the clusters whose abundance (size) is greater than 5%\n"
                   f"of the sum of sizes for each sample with {loc2trim2a}, enter 5\n"
                   f"Enter an integer between 0 and 100\n"
                   f"no default: ")
        while ts not in ["end", "home", "exit"] and ts.isnumeric() is False or ts == "":
            ts = input("\n-->ERROR: THRESHOLD format is not valid, it must be an integer between 0 and 100\n"
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
                sys.stdout.write("\n\tBye!...\n\n")
                quit()

    if rmenu == "2b":
        ts = input(f"\nEnter the THRESHOLD do you want to use for this locus {loc2trim2b}?\n"
                   f"Example: if you want to keep only the clusters whose abundance (size) is greater than 5%\n"
                   f"of the sum of sizes for each sample with {loc2trim2b}, enter 5\n"
                   f"Enter an integer between 0 and 100\n"
                   f"no default: ")
        while ts not in ["end", "home", "exit"] and ts.isnumeric() is False or ts == "":
            ts = input("\n-->ERROR: THRESHOLD format is not valid, it must be an integer between 0 and 100\n"
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
                sys.stdout.write("\n\tBye!...\n\n")
                quit()

    if rmenu == "2c":
        ts = input(f"\nEnter the THRESHOLD do you want to use for this locus {loc2trim2c}?\n"
                   f"Example: if you want to keep only the clusters whose abundance (size) is greater than 5%\n"
                   f"of the sum of sizes for each sample with {loc2trim2c}, enter 5\n"
                   f"Enter an integer between 0 and 100\n"
                   f"no default: ")
        while ts not in ["end", "home", "exit"] and ts.isnumeric() is False or ts == "":
            ts = input("\n-->ERROR: THRESHOLD format is not valid, it must be an integer between 0 and 100\n"
                       f"Enter a new THRESHOLD for {loc2trim2c}, (e.g. 5)\n"
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
                sys.stdout.write("\n\tBye!...\n\n")
                quit()

    if rmenu == "2d":
        ts = input(f"\nEnter the THRESHOLD do you want to use for this locus {loc2trim2d}?\n"
                   f"Example: if you want to keep only the clusters whose abundance (size) is greater than 5%\n"
                   f"of the sum of sizes for each sample with {loc2trim2d}, enter 5\n"
                   f"Enter an integer between 0 and 100\n"
                   f"no default: ")
        while ts not in ["end", "home", "exit"] and ts.isnumeric() is False or ts == "":
            ts = input("\n-->ERROR: THRESHOLD format is not valid, it must be an integer between 0 and 100\n"
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
                sys.stdout.write("\n\tBye!...\n\n")
                quit()
    return ts, ts1


def in_trim_sample2c():
    global samples, orient2c
    orient2c = input(f"\nEnter (new) SAMPLE do you want to trim among {samples}\n"
                     f"among {samples}?\n"
                     f"If you finished with locus {loc2trim2c} enter 'end': ")
    while orient2c not in samples and orient2c not in ["end", "home", "exit"]:
        orient2c = input(f"\n-->ERROR: sample name '{orient2c}' is not valid, please enter a valid name\n"
                         f"among {samples}\n"
                         "OR 'end' 'home' 'exit': ")
    else:
        if orient2c == "end":
            trim_2x()
        if orient2c == "home":
            main_menu2()
        if orient2c == "exit":
            sys.stdout.write("\n\tBye!...\n\n")
            quit()
    return orient2c


def in_trim_sample2d():
    global samples, orient2d
    orient2d = input(f"\nEnter (new) SAMPLE do you want to trim among {samples}\n"
                     f"among {samples}?\n"
                     f"If you finished with locus {loc2trim2d} enter 'end': ")
    while orient2d not in samples and orient2d not in ["end", "home", "exit"]:
        orient2d = input(f"\n-->ERROR: sample name '{orient2d}' is not valid, please enter a valid name\n"
                         f"among {samples}\n"
                         "OR 'end' 'home' 'exit': ")
    else:
        if orient2d == "end":
            trim_2x()
        if orient2d == "home":
            main_menu2()
        if orient2d == "exit":
            sys.stdout.write("\n\tBye!...\n\n")
            quit()
    return orient2d


# PARAMETER FILES #####################################################################################################
def param_1x():  # menu1 12
    """ Creates a file with one parameter by line for options 1x
    """
    global dir_fastq, sample_user, minsize_user, rmenu
    if rmenu == '1':
        with open("outputs/parameters_option_1.txt", "w") as out1:
            out1.write(f"Run option 1: {date}\n{dir_fastq}\n{fastqr1_user}\n{fastqr2_user}\n{loci1_user}\n"
                       f"{loci2_user}\n{sample_user}\n{minsize_user}\n{minseqlength}\n{alpha}\n{identity}\n"
                       f"Samples used = {samples}\nLoci paired-end used = {loci1s}\nLoci single-end (R1) "
                       f"used = {loci2s}\n")
    elif rmenu == '1a':
        with open("outputs/parameters_option_1a.txt", "w") as out2:
            out2.write(f"Run option 1a: {date}\n{dir_fastq}\n{fastqr1_user}\n{fastqr2_user}\n{loci1_user}\n"
                       f"{loci2_user}\n{sample_user}\n{minsize_user}\n{minseqlength}\n{alpha}\n{identity}\n"
                       f"Samples used = {samples}\nLoci paired-end used = {loci1s}\nLoci single-end (R1) "
                       f"used = {loci2s}\n")
    elif rmenu == '1b':
        with open("outputs/parameters_option_1b.txt", "w") as out3:
            out3.write(f"Run option 1b: {date}\n{dir_fastq}\n{fastqr1_user}\n{fastqr2_user}\n{loc_sel1}\n"
                       f"{sample_user}\n{minsize_user}\n{minseqlength}\n{alpha}\n{identity}\n")

    elif rmenu == '1c':
        with open("outputs/parameters_option_1c.txt", "w") as out4:
            out4.writelines(f"Run option 1c: {date}\n{dir_fastq}\n{fastqr1_user}\n{fastqr2_user}\n{loc_sel2}\n"
                            f"{sample_user}\n{minsize_user}\n{minseqlength}\n{alpha}\n{identity}\n")

    elif rmenu == '1d':
        with open("outputs/parameters_option_1d.txt", "w") as out5:
            out5.writelines(f"Run option 1d: {date}\n{dir_fastq}\n{fastqr1_user}\n{fastqr2_user}\n{sam_sel}\n"
                            f"{minsize_user}\n{minseqlength}\n{alpha}\n{identity}\n")


def prev_param():  # menu1a 2
    """ Recalls global variables for different options
    """
    global fastqr1s, fastqr2s, loci1s, loci2s, samples
    os.chdir(current_dir)
    with open("outputs/parameters_option_1.txt", "r") as infile:
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

    with open(fastqr1_user, "r") as out1:
        fastqr1s = out1.read().splitlines()

    with open(fastqr2_user, "r") as out2:
        fastqr2s = out2.read().splitlines()

    with open(loci1_user, "r") as out3:
        loci1s = out3.read().splitlines()

    with open(loci2_user, "r") as out4:
        loci2s = out4.read().splitlines()

    with open(sample_user, "r") as out5:
        samples = out5.read().splitlines()

    return fastqr1s, fastqr2s, loci1s, loci2s, samples


# STATISTICS ON FASTQ FILES ###########################################################################################
def quality():
    """Tests the quality of each 'fastq' file by the VSEARCH command:
    vsearch --fastq_eestats2 ../tmp_files/fastqF-R1 --output ../tmp_files/SN_R1_quality.txt
    vsearch --fastq_eestats2 ../tmp_files/fastqF-R2 --output ../tmp_files/SN_R2_quality.txt

    fastqF = fastq file name
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


# STEP 1: MERGING FOR AMPLICONS BASED ON PAIRED-END READS #############################################################
def merging():  # menu1 17
    """Merges paired-end reads into one sequence, when the length of the expected amplicon allows it
    according the VSEARCH command:
    vsearch --fastq_mergepairs ../tmp_files/fastqR1 --reverse ../tmp_files/fastqR2 --fastaout
    ./tmp_files/SN.fa --fastq_allowmergestagger --relabel sample='sample=SN'_merged

    fastqR1 = fastqR1 complete name
    fastqR2 = fastqR2 complete name
    option --fastq_allowmergestagger allows the merging of short fragments
    """
    with open("scripts/merging." + scriptExt, "w") as out:
        i = 0
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Merging paired-end reads for each sample:\n'))
        while i < len(samples):
            sample = samples[i]
            fastqr1 = fastqr1s[i]
            fastqr2 = fastqr2s[i]
            out.write(main_stream_message(
                f" {sample}...") +
                f"vsearch --fastq_mergepairs {dir_fastq}{fileSep}{fastqr1} "
                f"--reverse {dir_fastq}{fileSep}{fastqr2} --fastaout "
                f"../tmp_files/{sample}.fa --fastq_allowmergestagger --fastq_maxee 1 --relabel "
                f"sample={sample}_merged" + localErrorOnStopCmd + "\n")
            i = i + 1
        out.write(main_stream_message(f'\n\n'))


# FASTQ TO FASTA FOR AMPLICONS BASED ON SINGLE-END READS (R1 only) ####################################################
def fastq2fas():  # menu1 18
    """When the merging R1/R2 is impossible because of an unadapted size of amplicon, the reads R1 of 301 bp
    (better than R2) are used to search the relevant sequences.
    First, all R1 'fastq' files have to be transformed into 'fasta' files by the VSEARCH command:
    vsearch --fastq_filter ../tmp_files/fastaqR1 --fastaout ../tmp_files/SN_R1.fa

    SN = sample name
    fastqR1 = fastqR1 file name
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
                f"vsearch --fastq_filter {dir_fastq}{fileSep}{fastqr1} --fastaout ../tmp_files/{sample}_R1.fa" +
                localErrorOnStopCmd + "\n")
            i = i + 1
        out.write(main_stream_message(f'\n\n'))


# DEREPLICATION #######################################################################################################
def derep_1():  # menu1 19
    """Dereplicates merged sequences in a given 'fasta' file with the VSEARCH command:
    vsearch --derep_fulllength ../tmp_files/SN.fa --output ../tmp_files/SN_derep.fas --sizeout
    --strand both

    Dereplicates R1 sequences in a given 'fasta' file with the VSEARCH command:
    vsearch --derep_fulllength ../tmp_files/_R1_SN'.fa --output ../tmp_files/SN_derep.fas --sizeout
    --strand both --relabel sample='SN'_merged.

    Dereplicates in both strands and writes abundance annotation (frequency) to output.
    SN = sample name
    """
    with open("scripts/derep." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Dereplicating merged reads for each sample:\n'))
        for sample in samples:
            out.write(main_stream_message(
                f' {sample}...') +
                f"vsearch --derep_fulllength ../tmp_files/{sample}.fa --output ../tmp_files/{sample}_derep.fas "
                f"--sizeout --strand both" + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))

    with open("scripts/derep_r1." + scriptExt, "w") as out1:
        out1.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Dereplicating single-end reads '
                                                                     f'for samples:\n'))
        for sample in samples:
            out1.write(main_stream_message(f' {sample}...') +
                       f"vsearch --derep_fulllength ../tmp_files/{sample}_R1.fa --output ../tmp_files/{sample}"
                       f"_derep_R1.fas --sizeout --strand both --relabel sample={sample}_R1." + localErrorOnStopCmd
                       + "\n")
        out1.write(main_stream_message(f'\n\n'))


# UNOISING AND CLUSTERING #############################################################################################
def cluster_1x():
    """Denoises and clusters Illumina dereplicated merged sequences and gives in output the centroids sequences
    to 'fasta' files with the following VSEARCH commands:

    For loci based in paired-end reads:
    vsearch --cluster_unoise ../tmp_files/SN_derep.fas --sizein --centroids ../tmp_files/SN_cluster.fas
    --strand both --minsize int --sizeout --sizeorder --unoise_alph int --minseqlength int

    For loci based in single-end reads:
    vsearch --cluster_unoise ../tmp_files/SN_derep.fas --sizein --centroids ../tmp_files/SN_cluster.fas
    --strand both --minsize int --sizeout --sizeorder --unoise_alph int --minseqlength int

    For a single selected sample:
    vsearch --cluster_unoise ../tmp_files/SS_derep.fas --sizein --centroids ../tmp_files/SS_cluster.fas
    --strand both --minsize int --sizeout --sizeorder --unoise_alph int --minseqlength int
    For loci based in single-end reads

    SN = sample name; int = integer; SS = selected sample
    """
    global rmenu
    if rmenu in ["1", "1a", "1b"]:
        with open("scripts/cluster." + scriptExt, "w") as out:
            i = 0
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Clustering merged reads for each sample:\n'))
            while i < len(samples):
                sample = samples[i]
                out.write(main_stream_message(
                    f' {sample}...') +
                    f"vsearch --cluster_unoise ../tmp_files/{sample}_derep.fas --sizein --centroids "
                    f"../tmp_files/{sample}_cluster.fas --strand both --minsize {minsize_user} --sizeout "
                    f"--sizeorder "
                    f"--unoise_alph {alpha} --minseqlength {minseqlength}" + localErrorOnStopCmd + "\n")
                i = i + 1
            out.write(main_stream_message(f'\n\n'))

    if rmenu in ["1", "1a", "1c"]:
        with open("scripts/cluster_r1." + scriptExt, "w") as out:
            i = 0
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Clustering single-end reads '
                                                                        f'for each sample:\n'))
            while i < len(samples):
                sample = samples[i]
                out.write(main_stream_message(
                    f' {sample}...') +
                    f"vsearch --cluster_unoise ../tmp_files/{sample}_derep_R1.fas --sizein --centroids "
                    f"../tmp_files/{sample}_cluster_R1.fas --strand both --minsize {minsize_user} --sizeout "
                    f"--sizeorder --unoise_alpha {alpha} --minseqlength {minseqlength}" + localErrorOnStopCmd
                    + "\n")
                i = i + 1
            out.write(main_stream_message(f'\n\n'))

    elif rmenu == "1d":
        with open("scripts/cluster_one_sample_1d." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(
                f'\nClustering reads for '
                f'selected sample {sam_sel}:') +
                f"vsearch --cluster_unoise ../tmp_files/{sam_sel}_derep.fas --sizein --centroids "
                f"../tmp_files/{sam_sel}_cluster.fas --strand both --minsize {minsize_user} --sizeout "
                f"--sizeorder --unoise_alph {alpha} --minseqlength {minseqlength}" + localErrorOnStopCmd + "\n"
                f"vsearch --cluster_unoise ../tmp_files/{sam_sel}_derep_R1.fas --sizein --centroids "
                f"../tmp_files/{sam_sel}_cluster_R1.fas --strand both --minsize {minsize_user} "
                f"--sizeout --sizeorder --unoise_alph {alpha} --minseqlength {minseqlength}"
                + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))


# REMOVING CHIMERA ####################################################################################################
def chimera_remove():
    """Detects and removes potential chimeras in denoised merged sequences or single-end or selected sample
    by the VSEARCH commands:
    vsearch --uchime3_denovo ../tmp_files/SN_cluster.fas --nonchimeras ../tmp_files/SN_cluster_OK.fas
    vsearch --uchime3_denovo ../tmp_files/SN_cluster-R1.fas --nonchimeras ../tmp_files/SN_cluster_R1_OK.fas
    vsearch --uchime3_denovo ../tmp_files/SS_cluster.fas --nonchimeras ../tmp_files/SS_cluster_OK.fas
    vsearch --uchime3_denovo ../tmp_files/SS_cluster_R1.fas --nonchimeras ../tmp_files/SS_cluster_R1_OK.fas
    SN = sample name
    ss = selected sample
    After denoising, the chimeras are very scarce
    """
    if rmenu in ["1", "1a", "1b"]:
        with open("scripts/chimera." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd +
                    "\n" + main_stream_message(f'Detecting and removing chimeras within merged reads of samples:\n'))
            for sample in samples:
                out.write(main_stream_message(
                    f' {sample}...') +
                    f"vsearch --uchime3_denovo ../tmp_files/{sample}_cluster.fas --nonchimeras "
                    f"../tmp_files/{sample}_cluster_OK.fas" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

    if rmenu in ["1", "1a", "1c"]:
        with open("scripts/chimera_r1." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Detecting and removing chimeras within '
                                                                        f'single-end reads of samples:\n'))
            for sample in samples:
                out.write(main_stream_message(
                    f' {sample}...') +
                    f"vsearch --uchime3_denovo ../tmp_files/{sample}_cluster_R1.fas --nonchimeras"
                    f" ../tmp_files/{sample}_cluster_R1_OK.fas" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

    if rmenu == "1d":
        with open("scripts/chimera_one_sample_1d." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" + main_stream_message(
                f'Detecting and removing chimeras for selected sample {sam_sel}') +
                f"vsearch --uchime3_denovo ../tmp_files/{sam_sel}_cluster.fas --nonchimeras "
                f"../tmp_files/{sam_sel}_cluster_OK.fas" + localErrorOnStopCmd +
                "\n"
                f"vsearch --uchime3_denovo ../tmp_files/{sam_sel}_cluster_R1.fas --nonchimeras "
                f"../tmp_files/{sam_sel}_cluster_R1_OK.fas" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))


# DISTRIBUTION OF CLUSTERS INTO THE DIFFERENT LOCI ####################################################################
def runloc_merged():  # menu1 25
    """Searches similarities between merged, denoised and non-chimera sequences and the local reference
    database (-db) by the VSEARCH command:
    vsearch --usearch_global ../tmp_files/SN_cluster_OK.fas --db ../refs/L1.fas --matched ../L1/SN_merged.fas
    --id int --strand both

    L1 = locus name for amplicons based on paired-end reads
    This is the second crucial step metabarcoding analyze, depending on a required minimal identity setup by
    the option 'identity'.
    """
    with open("scripts/locimerged." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd +
                  "\n" + main_stream_message(f'Affiliating clusters to loci for merged reads of samples:\n'))
        for loci1b in loci1s:
            for sample in samples:
                out.write(main_stream_message(
                    f' {sample} VS {loci1b}...') +
                    f"vsearch --usearch_global ../tmp_files/{sample}_cluster_OK.fas --db ../refs/{loci1b}.fas"
                    f" --matched ../{loci1b}/{sample}_merged.fas --id {identity} --strand both"
                    + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))


def runloc_r1():  # menu1 26
    """Similar to runloc_merged but with the R1 denoised and non-chimera sequences
    in case of amplicons where the merging R1/R2 is impossible, with the VSEARCH command:
    vsearch --usearch_global ../tmp_files/SN_cluster_R1_OK.fas --db ../refs/L2.fas --matched ../L2/SN_R1.fas
     --id real --strand both

     L1 = locus for amplicons with no mergeable R1/R2
     SN = sample name
     real = a real from 0 to 1, generally around 0.7
    """
    with open("scripts/locir1." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd +
                  "\n" + main_stream_message(f'Affiliating clusters to loci for single-end reads of samples:\n'))
        for locus2b in loci2s:
            for sample in samples:
                out.write(main_stream_message(f' {sample} VS {locus2b}...') +
                            f"vsearch --usearch_global ../tmp_files/{sample}_cluster_R1_OK.fas --db "
                            f"../refs/{locus2b}.fas --matched ../{locus2b}/{sample}_R1.fas --id {identity} "
                            f"--strand both" + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))


# DISTRIBUTION OF CLUSTERS FOR ONE SELECTED LOCUS #####################################################################
def runlocsel_merged():
    """Similar to runloc_merged but with only one selected locus with the VSEARCH command:
    vsearch --usearch_global ../tmp_files/SN_cluster_OK.fas --db ../refs/L1.fas --matched
    ../L2/SN_merged.fas --id real --strand both

    SN = sample name
    L1 = locus for amplicons based on paired-end reads
    """
    with open("scripts/loci_sel." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters to selected '
                                                                      f'locus {loc_sel1} for each sample:\n'))
        for sample in samples:
            out.write(main_stream_message(f' {sample}...') +
                        f"vsearch --usearch_global ../tmp_files/{sample}_cluster_OK.fas --db ../refs/{loc_sel1}.fas"
                        f" --matched ../{loc_sel1}/{sample}_merged.fas --id {identity} --strand both"
                        + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))


def runlocsel_r1():
    """identical to runlocsel_merged but with R1 denoised and non-chimera sequences with the VSEARCH command:
    vsearch --usearch_global ../tmp_files/SN_cluster_R1_OK.fas --db ../refs/L2.fas --matched
    ../L2/SN_R1.fas --id real --strand both

    SN = sample name
    L2 = locus name for amplicons with no mergeable R1/R2
    """
    with open("scripts/locir1_sel." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters to selected '
                                                                      f'locus {loc_sel2} for each sample:\n'))
        for sample in samples:
            out.write(main_stream_message(f' {sample}...') +
                        f"vsearch --usearch_global ../tmp_files/{sample}_cluster_R1_OK.fas --db ../refs/{loc_sel2}.fas"
                        f" --matched ../{loc_sel2}/{sample}_R1.fas --id {identity} --strand both"
                        + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))


# DISTRIBUTION OF CLUSTERS FOR ONE SELECTED SAMPLE ####################################################################
def runloc_one_sample_1d():
    """ Searches similarities between denoised and non-chimera sequences and local reference database (db)
    by the VSEARCH command, but only for a selected sample:
    vsearch --usearch_global ../tmp_files/SS_cluster_OK.fas --db ../refs/L1.fas --matched ../L1/SS_merged.fas
     --id real --strand both
    vsearch --usearch_global ../tmp_files/SS_cluster_OK.fas --db ../refs/L2.fas --matched ../L2/SS_merged.fas
     --id real --strand both

    SN = locus name
    L1 = locus name for amplicons based on paired-end reads
    L2 = locus name for amplicons with no mergeable R1/R2
    real = real from 0 to 1
    This is the second crucial step metabarcoding analyze, depending on a required minimal identity setup by
    the option 'identity'
    """
    with open("scripts/loci_merged_1d." + scriptExt, "w") as out:
        out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters of selected '
                                                                      f'sample {sam_sel} to the loci:\n'))
        for loci1b in loci1s:
            out.write(main_stream_message(f' {loci1b}...') +
                        f"vsearch --usearch_global ../tmp_files/{sam_sel}_cluster_OK.fas --db ../refs/{loci1b}.fas"
                        f" --matched ../{loci1b}/{sam_sel}_merged.fas --id {identity} --strand both"
                        + localErrorOnStopCmd + "\n")
        out.write(main_stream_message(f'\n\n'))

    with open("scripts/loci_R1_1d." + scriptExt, "w") as out1:
        out1.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters fo selected '
                                                                      f'sample {sam_sel} to the loci:\n'))
        for locus2b in loci2s:
            out1.write(main_stream_message(f' {locus2b}...') +
                        f"vsearch --usearch_global ../tmp_files/{sam_sel}_cluster_R1_OK.fas --db ../refs/{locus2b}.fas"
                        f" --matched ../{locus2b}/{sam_sel}_R1.fas --id {identity} --strand both"
                        + localErrorOnStopCmd + "\n")
            out1.write(main_stream_message(f'\n\n'))


# ORIENTATION OF SEQUENCES INTO FORWARD DIRECTION ####################################################################
def orient_1x():
    """Orients all the sequences in the same direction (forward) than references with the following script:

    For paired-end merged sequences:
    vsearch --orient ../L1/SN_merged.fas --db ../refs/L1.fas --fastaout ../L1/SN_orient.fas

    For single-end R1 sequences:
    vsearch --orient ../L2/SN_R1.fas --db ../refs/L2.fas --fastaout ../L2/SN_R1_orient.fas

    For selected loci:
    vsearch --orient ../SL/SN_merged.fas --db ../refs/SL.fas --fastaout ../SL/SN_orient.fas
    vsearch --orient ../SL/SN_R1.fas --db ../refs/SL.fas --fastaout ../SL/SN_R1_orient.fas

    For selected sample:
    vsearch --orient ../L1/SS_merged.fas --db ../refs/L1.fas --fastaout ../L1/SS_orient.fas
    vsearch --orient ../L2/SS_merged.fas --db ../refs/L2.fas --fastaout ../L1/SS_orient.fas

    SN = locus name
    L1 = locus name for amplicons based on paired-end reads
    L2 = locus name for amplicons with no mergeable R1/R2
    SL = selected locus
    """
    if rmenu in ["1", "1a"]:
        with open("scripts/orientloc1." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" +
                        main_stream_message(f"Orienting all merged reads in "
                        f"the same direction for each sample:\n"))
            for locus1b in loci1s:
                for sample in samples:
                    out.write(main_stream_message(f' {sample} VS {locus1b}...') +
                        f"vsearch --orient ../{locus1b}/{sample}_merged.fas --db ../refs/{locus1b}.fas --fastaout "
                        f"../{locus1b}/{sample}_orient.fas" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

        with open("scripts/orientloc2." + scriptExt, "w") as out1:
            out1.write(globalErrorOnStopCmd + "\n" +
                        main_stream_message(f"Orienting all single-end reads in "
                        f"the same direction for each sample:\n"))
            for locus2b in loci2s:
                for sample in samples:
                    out1.write(main_stream_message(f' {sample} VS {locus2b}...') +
                    f"vsearch --orient ../{locus2b}/{sample}_R1.fas --db ../refs/{locus2b}.fas --fastaout "
                    f"../{locus2b}/{sample}_R1_orient.fas" + localErrorOnStopCmd + "\n")
            out1.write(main_stream_message(f'\n\n'))

    elif rmenu == "1b":
        with open("scripts/orient_merged_1b." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" +
                main_stream_message(f'Orientation of all merged reads in the same '
                f'direction for selected locus {loc_sel1} '
                f'and for each sample:\n'))
            for sample in samples:
                out.write(main_stream_message(f' {sample}...') +
                f"vsearch --orient ../{loc_sel1}/{sample}_merged.fas --db ../refs/{loc_sel1}.fas --fastaout "
                f"../{loc_sel1}/{sample}_orient.fas" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

    elif rmenu == "1c":
        with open("scripts/orient_R1_1c." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" +
                main_stream_message(f'Orientation of all single-end reads in the '
                f'same direction for selected locus '
                f'{loc_sel2} and for each sample:\n'))
            for sample in samples:
                out.write(main_stream_message(f' {sample}...') +
                            f"vsearch --orient ../{loc_sel2}/{sample}_R1.fas --db ../refs/{loc_sel2}.fas --fastaout "
                            f"../{loc_sel2}/{sample}_R1_orient.fas" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

    if rmenu == "1d":
        with open("scripts/orient_merged_1d." + scriptExt, "w") as out:
            out.write(globalErrorOnStopCmd + "\n" +
                        main_stream_message(f'Orientation of all clusters of selected '
                        f'sample {sam_sel} for the loci:\n'))
            for locus1b in loci1s:
                out.write(main_stream_message(f' {locus1b}...') +
                f"vsearch --orient ../{locus1b}/{sam_sel}_merged.fas --db ../refs/{locus1b}.fas --fastaout "
                f"../{locus1b}/{sam_sel}_orient.fas" + localErrorOnStopCmd + "\n")
            out.write(main_stream_message(f'\n\n'))

        with open("scripts/orient_R1_1d." + scriptExt, "w") as out18:
            out18.write(main_stream_message(f'Orientation of all clusters of selected sample {sam_sel} for '
                                            f'the loci:\n'))
            for locus2b in loci2s:
                out18.write(main_stream_message(f' {locus2b}...') +
                            f"vsearch --orient ../{locus2b}/{sam_sel}_R1.fas --db ../refs/{locus2b}.fas --fastaout "
                            f"../{locus2b}/{sam_sel}_R1_orient.fas" + localErrorOnStopCmd + "\n")
            out18.write(main_stream_message(f'\n\n'))


# CREATION OF GLOBAL SCRIPTS ##########################################################################################
def runs_1x():
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
        print(f"Run of option {rmenu} was correctly achieved")

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
        print(f"Run of option {rmenu} was correctly achieved")

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
        print(f"Run of option {rmenu} was correctly achieved")

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
        print(f"Run of option {rmenu} was correctly achieved")

    if rmenu == "1d":
        with open(f"scripts/runall1d.{scriptExt}", "w") as out:
            if not winOS:
                out.write(start_log_redirect('../outputs/res1d.log') +
                "scriptArray=('"
                "./cluster_one_sample_1d." + scriptExt + "' '"
                "./chimera_one_sample_1d." + scriptExt + "' '"
                "./loci_merged_1d." + scriptExt + "' '"
                "./loci_R1_1d." + scriptExt + "' '"
                "./orient_merged_1." + scriptExt + "' '"
                ".orient_R1_1d." + scriptExt + "')\n" +
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
        print(f"Run of option {rmenu} was correctly achieved")

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
            print(f"\nMain analysis execution failed, please check {current_dir}/outputs/res1e.log")
            exit(1)
        print(f"Run of option {rmenu} was correctly achieved")
        sys.stdout.write("\nThe run option 1e is complete\n"
                         f"Statistical test files (*_quality.txt) are in\n"
                         f"Results in ---> {current_dir}/outputs\n\n")
        os.chdir(current_dir)


# GLOBAL STATISTICS ACCORDING THE OPTIONS ###########################################################################
def stats_1x():  # menu1 29 fin
    """Calculates the number of resulting sequences (reads, merged, dereplicates and clusters) according
    to the selected options 'minsize', minseqlength, alpha parameter' and 'identity'.
    """
    if rmenu == "1":
        os.chdir(current_dir)
        print("\nComputing statistics after using option 1, wait...\n"
              f"Results in ---> {current_dir}/outputs/Stats_option_n1.txt:")
        with open("outputs/Stats_option_1.txt", "w") as out:
            out.write("With option 1, parameters set to:\n\n"
                          f"Directory = {current_dir}\n"
                          f"Fastq R1 file name = {fastqr1_user}\n"
                          f"Fastq R2 file name = {fastqr2_user}\n"
                          f"Paired-end based loci = {loci1s}\n"
                          f"Single-end based (R1) loci = {loci2s}\n"
                          f"Samples = {samples}\n"
                          f"Minimum abundance for clusters = {minsize_user}\n"
                          f"Minimum length for sequences = {minseqlength}\n"
                          f"Alpha clustering parameter = {alpha}\n"
                          f"Identity for allocating clusters = {identity}\n\n")
            for sample in samples:
                sample = sample.rstrip()
            # Nb READS Calculated on R1.fa
                r1fa = open("tmp_files/" + sample + "_R1.fa", "rt")
                reads = r1fa.read()
                nb_reads = reads.count(">")
            # Nb merged .fa
                merged = open("tmp_files/" + sample + ".fa", "rt")
                mgd = merged.read()
                nb_merged = mgd.count(">")
                percent_merging = (nb_merged/nb_reads)*100
                percent = '{:.2f}%'.format(percent_merging)
            # Number of dereplicated sequences calculated on _derep.fas
                derepfas = open("tmp_files/" + sample + "_derep.fas", "rt")
                df = derepfas.read()
                nb_derpm = df.count("sample")
            # Nb clustermerged Calculated on _cluster.fas
                clusmerged = open("tmp_files/" + sample + "_cluster.fas", "rt")
                clusmer = clusmerged.read()
                nb_clusm = clusmer.count("sample")
            # Number of merged clusters without chimeras calculated on _cluster_OK.fas
                clusmergedok = open("tmp_files/" + sample + "_cluster_OK.fas", "rt")
                clusmerok = clusmergedok.read()
                nb_clusmok = clusmerok.count("sample")
            # Number of dereplicated sequences (R1) calculated on _derep_R1.fas
                derepfasr1 = open("tmp_files/" + sample + "_derep_R1.fas", "rt")
                dfr1 = derepfasr1.read()
                nb_derpr1 = dfr1.count("sample")
            # Nb clusterR1 Calculated on _cluster_R1.fas
                clusfasr1 = open("tmp_files/" + sample + "_cluster_R1.fas", "rt")
                clusr1 = clusfasr1.read()
                nb_clusr1 = clusr1.count("sample")
            # Nb clusterR1OK Calculated on _cluster_R1_OK.fas
                clusfasr1ok = open("tmp_files/" + sample + "_cluster_R1_OK.fas", "rt")
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
                    os.chdir(locus1)
                    refs = open(sample + "_merged.fas")
                    refsloc = refs.read()
                    nb_ref = refsloc.count("sample")
                    out.writelines(f"\t{nb_ref} clusters of merged sequences of {sample} affiliated "
                                    f"to locus {locus1}\n")
                    os.chdir("../")
                for locus2 in loci2s:
                    os.chdir(locus2)
                    refs2 = open(sample + "_R1.fas")
                    refsloc2 = refs2.read()
                    nb_ref2 = refsloc2.count("sample")
                    out.writelines(f"\t{nb_ref2} clusters of single-end (R1) sequences of {sample} affiliated "
                                    f"to locus {locus2}\n")
                    os.chdir(current_dir)
            print("Complete\n"
                  "\nThe MAIN MANDATORY ANALYSIS (option 1) is complete, you can continue by options 2x and 3x\n")

    elif rmenu == "1a":
        os.chdir(current_dir)
        print("\nComputing statistics using option 1a\n"
              f"Results in ---> {current_dir}/outputs/Stats_option_1a.txt:")
        with open("outputs/Stats_option_1a.txt", "w") as out:
            out.write("With option 1a, parameters set to:\n\n"
                       f"Directory = {current_dir}\n"
                       f"Fastq R1 file name = {fastqr1_user}\n"
                       f"Fastq R2 file name = {fastqr2_user}\n"
                       f"Paired-end based loci = {loci1s}\n"
                       f"Single-end based (R1) loci = {loci2s}\n"
                       f"Samples = {samples}\n"
                       f"Minimum abundance for clusters = {minsize_user}\n"
                       f"Minimum length for sequences = {minseqlength}\n"
                       f"Alpha clustering parameter = {alpha}\n"
                       f"Identity for allocating clusters = {identity}\n\n")

            for sample in samples:
                sample = sample.rstrip()
            # Nb READS Calculated on R1.fa
                r1fa = open("tmp_files/" + sample + "_R1.fa", "rt")
                reads = r1fa.read()
                nb_reads = reads.count(">")
            # Nb merged .fa
                merged = open("tmp_files/" + sample + ".fa", "rt")
                mgd = merged.read()
                nb_merged = mgd.count(">")
                percent_merging = (nb_merged/nb_reads)*100
                percent = '{:.2f}%'.format(percent_merging)
            # Number of dereplicated sequences calculated on _derep.fas
                derepfas = open("tmp_files/" + sample + "_derep.fas", "rt")
                df = derepfas.read()
                nb_derpm = df.count("sample")
            # Nb clustermerged Calculated on _cluster.fas
                clusmerged = open("tmp_files/" + sample + "_cluster.fas", "rt")
                clusmer = clusmerged.read()
                nb_clusm = clusmer.count("sample")
            # Number of merged clusters without chimeras calculated on _cluster_OK.fas
                clusmergedok = open("tmp_files/" + sample + "_cluster_OK.fas", "rt")
                clusmerok = clusmergedok.read()
                nb_clusmok = clusmerok.count("sample")
            # Number of dereplicated sequences (R1) calculated on _derep_R1.fas
                derepfasr1 = open("tmp_files/" + sample + "_derep_R1.fas", "rt")
                dfr1 = derepfasr1.read()
                nb_derpr1 = dfr1.count("sample")
            # Nb clusterR1 Calculated on _cluster_R1.fas
                clusfasr1 = open("tmp_files/" + sample + "_cluster_R1.fas", "rt")
                clusr1 = clusfasr1.read()
                nb_clusr1 = clusr1.count("sample")
            # Nb clusterR1OK Calculated on _cluster_R1_OK.fas
                clusfasr1ok = open("tmp_files/" + sample + "_cluster_R1_OK.fas", "rt")
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
                    os.chdir(locus1)
                    refs = open(sample + "_merged.fas")
                    refsloc = refs.read()
                    nb_ref = refsloc.count("sample")
                    out.writelines(f"\t{nb_ref} clusters of merged sequences of {sample} affiliated "
                                    f"to locus {locus1}\n")
                    os.chdir("../")
                for locus2 in loci2s:
                    os.chdir(locus2)
                    refs2 = open(sample + "_R1.fas")
                    refsloc2 = refs2.read()
                    nb_ref2 = refsloc2.count("sample")
                    out.writelines(f"\t{nb_ref2} clusters of single-end (R1) sequences of {sample} affiliated "
                                    f"to locus {locus2}\n")
                    os.chdir(current_dir)
            print("Complete\n")

    elif rmenu == "1b":
        os.chdir(current_dir)
        print("\nComputing statistics after using option 1b, wait...\n"
              f"Results in ---> {current_dir}/outputs/Stats_option_1b.txt:")
        with open("outputs/Stats_option_1b.txt", "w") as out:
            out.write("With option 1b, parameters set to:\n\n"
                       f"Directory = {current_dir}\n"
                       f"Fastq R1 file name = {fastqr1_user}\n"
                       f"Fastq R2 file name = {fastqr2_user}\n"
                       f"Paired-end based loci = {loci1s}\n"
                       f"Single-end based (R1) loci = {loci2s}\n"
                       f"Samples = {samples}\n"
                       f"Minimum abundance for clusters = {minsize_user}\n"
                       f"Minimum length for sequences = {minseqlength}\n"
                       f"Alpha clustering parameter = {alpha}\n"
                       f"Identity for allocating clusters = {identity}\n\n"
                       f"The selected locus is {loc_sel1}\n"
                       f"For {loc_sel1}:\n")
            os.chdir(loc_sel1)
            for sample in samples:
                a = open(sample + "_merged.fas")
                b = a.read()
                c = b.count("sample")
                out.write(f"\t{sample} has {c} clusters\n")
            os.chdir(current_dir)
        print("Complete\n"
              "\nThe run option 1b is complete, you can continue by options 2x and 3x\n")

    elif rmenu == "1c":
        os.chdir(current_dir)
        print("\nComputing statistics after using option 1c, wait...\n"
              f"Results in ---> {current_dir}/outputs/Stats_option_1c.txt")
        with open("outputs/Stats_option_1c.txt", "w") as out:
            out.write("With option 1c, parameters set to:\n\n"
                       f"Directory = {current_dir}\n"
                       f"Fastq R1 file name = {fastqr1_user}\n"
                       f"Fastq R2 file name = {fastqr2_user}\n"
                       f"Paired-end based loci = {loci1s}\n"
                       f"Single-end based (R1) loci = {loci2s}\n"
                       f"Samples = {samples}\n"
                       f"Minimum abundance for clusters = {minsize_user}\n"
                       f"Minimum length for sequences = {minseqlength}\n"
                       f"Alpha clustering parameter = {alpha}\n"
                       f"Identity for allocating clusters = {identity}\n\n"
                       f"The selected locus is {loc_sel2}\n"
                       f"For {loc_sel2}:\n")
            os.chdir(loc_sel2)
            for sample in samples:
                a = open(sample + "_R1.fas")
                b = a.read()
                c = b.count("sample")
                out.write(f"\t{sample} has {c} clusters\n")
            os.chdir(current_dir)
        print("Complete\n"
              "\nThe run option 1c is complete, you can continue by options 2x and 3x\n")

    elif rmenu == "1d":
        os.chdir(current_dir)
        print("\n\nComputing statistics after using option 1d, wait\n"
              f"Results in ---> {current_dir}/outputs/Stats_option_1d.txt:")
        with open("outputs/Stats_option_1d.txt", "w") as out:
            out.write("With option 1d, parameters set to:\n\n"
                       f"Directory = {current_dir}\n"
                       f"Fastq R1 file name = {fastqr1_user}\n"
                       f"Fastq R2 file name = {fastqr2_user}\n"
                       f"Paired-end based loci = {loci1s}\n"
                       f"Single-end based (R1) loci = {loci2s}\n"
                       f"Samples = {samples}\n"
                       f"Minimum abundance for clusters = {minsize_user}\n"
                       f"Minimum length for sequences = {minseqlength}\n"
                       f"Alpha clustering parameter = {alpha}\n"
                       f"Identity for allocating clusters = {identity}\n\n"
                       f"The selected sample is {sam_sel}\n"
                       f"For {sam_sel}:\n")

            for locus1 in loci1s:
                os.chdir(locus1)
                a = open(sam_sel + "_merged.fas")
                b = a.read()
                c = b.count("sample")
                out.write(f"\t{sam_sel} has {c} clusters for locus {locus1}\n")
                os.chdir(current_dir)

            for locus2 in loci2s:
                os.chdir(locus2)
                a = open(sam_sel + "_R1.fas")
                b = a.read()
                c = b.count("sample")
                out.write(f"\t{sam_sel} has {c} clusters for locus {locus2}\n")
                os.chdir(current_dir)
        print("Complete\n"
              "\nThe run option 1d is complete, you can continue by options 2x and 3x\n")


# TRIM PRIMERS AND SELECT SEQUENCES ACCORDING THRESHOLDS ##############################################################
def trim_2x():
    global loc2trim2a, loc2trim2b, loc2trim2c, loc2trim2d, ts, ts1, trim_left, trim_right, orient2c, orient2d
    if rmenu == "2a":
        while True:
            loc2trim2a = in_loc2trim_2x()
            os.chdir(loc2trim2a)
            trim_left = in_trim_left()
            trim_right = in_trim_right()
            ts = in_ts()
            stat_2a = open(f"../outputs/Stats_option_2a.txt", "w")
            stat_2a.write(f"Locus {loc2trim2a} trimmed {trim_left} bp (left) and {trim_right} bp (right) "
                          f"with threshold set at {ts1} percent\n")
            for sample in samples:
                with open(f"{sample}_orient.fas", 'r') as filin, open("trim-select." + scriptExt, "w") as out:
                    targets = [line for line in filin if "size" in line]
                    a = 0
                    for target in targets:
                        size = re.search('size=(.+?)$', target).group(1)
                        a = a + int(size)
                    b = int(a * float(ts1) + 1)
                    stat_2a.writelines(f"\tSum of sizes for {sample} = {a}\n"
                                       f"\tThe sizes > {b} for {sample} were conserved\n")
                    out.write(f"" + start_log_redirect('./' + sample + '.log') +
                                f"vsearch --fastx_filter {sample}_orient.fas --fastq_stripleft {trim_left} "
                                f"  --fastq_stripright {trim_right} --fastaout tmp --minsize {b}"
                                + localErrorOnStopCmd + "\n"
                                f"vsearch --derep_fulllength ./tmp --output {sample}_select.fas --sizein --sizeout"
                                + localErrorOnStopCmd + "\n"
                                f"" + end_log_redirect('./' + sample + '.log'))
                subprocess.run([shellCmd, "./trim-select." + scriptExt])
                selected = open(f'./{sample}_select.fas', 'r')
                nb_selected = selected.read().count('>')
                stat_2a.writelines(f"\tNumber of selected clusters for sample {sample}: {nb_selected}\n\n")
                sys.stdout.write(f"\n\nSum of sizes for {sample} at locus {loc2trim2a} = {a}\n"
                                 f"With threshold set at {ts}, sizes > {b} were conserved\n"
                                 f"Number of selected clusters for sample {sample}: {nb_selected}\n")
            os.remove("trim-select." + scriptExt)
            stat_2a.close()
            os.chdir(current_dir)

    if rmenu == "2b":
        while True:
            loc2trim2b = in_loc2trim_2x()
            os.chdir(loc2trim2b)
            trim_left = in_trim_left()
            trim_right = in_trim_right()
            ts = in_ts()
            stat_2b = open(f"../outputs/Stats_option_2b.txt", "w")
            stat_2b.write(f"Locus {loc2trim2b} trimmed {trim_left} bp (left) and {trim_right} bp (right) "
                          f"with threshold set at {ts} percent\n")
            for sample in samples:
                with open(sample + "_R1.fas", 'r') as filin, open("trim-select." + scriptExt, "w") as out:
                    targets = [line for line in filin if "size" in line]
                    a = 0
                    for target in targets:
                        size = re.search('size=(.+?)$', target).group(1)
                        a = a + int(size)
                    b = int(a * float(ts1) + 1)
                    stat_2b.writelines(f"\tSum of sizes for {sample} = {a}\n"
                                       f"\tThe sizes > {b} for {sample} were conserved\n")
                    out.write(f"" + start_log_redirect('./' + sample + '.log') +
                              f"vsearch --fastx_filter {sample}_R1_orient.fas --fastq_stripleft {trim_left} "
                              f"  --fastq_stripright {trim_right} --fastaout tmp --minsize {b}"
                              + localErrorOnStopCmd + "\n"
                              f"vsearch --derep_fulllength ./tmp --output {sample}_R1_select.fas --sizein --sizeout"
                              + localErrorOnStopCmd + "\n"
                                                      f"" + end_log_redirect('./' + sample + '.log'))
                subprocess.run([shellCmd, "./trim-select." + scriptExt])
                selected = open('./' + sample + '_R1_select.fas', 'r')
                nb_selected = selected.read().count('>')
                stat_2b.writelines(f"\tNumber of selected clusters for sample {sample}: {nb_selected}\n\n")
                sys.stdout.write(f"\nSum of sizes for {sample} at locus {loc2trim2b} = {a}\n"
                                 f"With threshold set at {ts}, sizes > {b} were conserved\n"
                                 f"Number of selected clusters for sample {sample}: {nb_selected}\n")
            os.remove("trim-select." + scriptExt)
            stat_2b.close()
            os.chdir(current_dir)

    if rmenu == "2c":
        while True:
            loc2trim2c = in_loc2trim_2x()
            os.chdir(f"{current_dir}/{loc2trim2c}")
            trim_left = in_trim_left()
            trim_right = in_trim_right()
            with open(f"../outputs/Stats_option_2d.txt", "w") as out:
                out.write(f"{loc2trim2c}: trimmed of {trim_left} bp (left) and {trim_right} bp (right)\n")
                while True:
                    orient2c = in_trim_sample2c()
                    ts = in_ts()
                    with open(orient2c + "_orient.fas", 'r') as filin, open("trim-select." + scriptExt,
                                                                                 "w") as out20:
                        targets = [line for line in filin if "size" in line]
                        a = 0
                        for target in targets:
                            size = re.search('size=(.+?)$', target).group(1)
                            a = a + int(size)
                        b = int(a * float(ts1) + 1)
                        out.writelines(f"\tSum of sizes for {loc2trim2c} = {a}\n"
                                       f"\tAt threshold {ts} sizes > {b} for {loc2trim2c} were conserved\n")
                        out20.writelines(f'' + start_log_redirect('./' + loc2trim2c + '.log') +
                            f' vsearch --fastx_filter {loc2trim2c}_orient.fas --fastq_stripleft {trim_left} '
                            f'  --fastq_stripright {trim_right} --fastaout tmp --minsize {b}\n'
                            + localErrorOnStopCmd + "\n"
                            f' vsearch --derep_fulllength ./tmp --output {loc2trim2c}_R1_select.fas --sizein '
                            f'--sizeout\n' + localErrorOnStopCmd + "\n"
                            f'' + end_log_redirect(
                            './' + orient2c + '.log'))
                    subprocess.run([shellCmd, "./trim-select." + scriptExt])
                    selected = open("./" + orient2c + "_select.fas", "r")
                    nb_selected = selected.read().count(">")
                    out.writelines(f"\t{loc2trim2c}: number of selected clusters is: {nb_selected}\n\n")
                    sys.stdout.write(f"\n{loc2trim2c}: sum of sizes = {a}\n"
                                     f"The sizes > {b} were conserved\n"
                                     f"{loc2trim2c}:number of selected clusters = {nb_selected}\n")
                    os.remove("trim-select." + scriptExt)
            out.close()
            os.chdir(current_dir)

    if rmenu == "2d":
        while True:
            loc2trim2d = in_loc2trim_2x()
            os.chdir(f"{current_dir}/{loc2trim2d}")
            trim_left = in_trim_left()
            trim_right = in_trim_right()
            with open(f"../outputs/Stats_option_2d.txt", "w") as out:
                out.write(f"{loc2trim2d}: trimmed of {trim_left} bp (left) and {trim_right} bp (right)\n")
                while True:
                    orient2d = in_trim_sample2d()
                    ts = in_ts()
                    with open(orient2d + "_R1_orient.fas", 'r') as filin, open("trim-select." + scriptExt,
                                                                                 "w") as out20:
                        targets = [line for line in filin if "size" in line]
                        a = 0
                        for target in targets:
                            size = re.search('size=(.+?)$', target).group(1)
                            a = a + int(size)
                        b = int(a * float(ts1) + 1)
                        out.writelines(f"\tSum of sizes for {loc2trim2d} = {a}\n"
                                       f"\tAt threshold {ts} sizes > {b} for {loc2trim2d} were conserved\n")
                        out20.writelines(f'' + start_log_redirect('./' + loc2trim2d + '.log') +
                            f' vsearch --fastx_filter {loc2trim2d}_R1_orient.fas --fastq_stripleft {trim_left} '
                            f'  --fastq_stripright {trim_right} --fastaout tmp --minsize {b}\n'
                            + localErrorOnStopCmd + "\n"
                            f' vsearch --derep_fulllength ./tmp --output {loc2trim2d}_R1_select.fas --sizein '
                            f'--sizeout\n' + localErrorOnStopCmd + "\n"
                            f'' + end_log_redirect(
                            './' + orient2d + '.log'))
                    subprocess.run([shellCmd, "./trim-select." + scriptExt])
                    selected = open("./" + orient2d + "_R1_select.fas", "r")
                    nb_selected = selected.read().count(">")
                    out.writelines(f"\t{loc2trim2d}: number of selected clusters is: {nb_selected}\n\n")
                    sys.stdout.write(f"\n{loc2trim2d}: sum of sizes = {a}\n"
                                     f"The sizes > {b} were conserved\n"
                                     f"{loc2trim2d}:number of selected clusters = {nb_selected}\n")
                    os.remove("trim-select." + scriptExt)
            out.close()
            os.chdir(current_dir)


# CLUSTERING OF ALL SAMPLES BY LOCUS FOR PHYLOGENETIC ANALYZES ########################################################
def concat_3():
    """CLUSTERING OF ALL SAMPLES BY LOCUS FOR PHYLOGENETIC ANALYZES PAIRED_END
    """
    os.system("cls" if winOS else "clear")
    global alloci, samples, loc2cat
    nb_samples = len(samples)
    alloci = loci1s + list(set(loci2s) - set(loci1s))
    if os.path.exists("outputs/Stats_option_3.txt"):
        os.remove("outputs/Stats_option_3.txt")
    while True:
        loc2cat = input("\n----- CONCATENATION OF ALL SAMPLES BY LOCUS FOR PHYLOGENETIC ANALYZES -----\n"
                        f"\nFor which (new) LOCUS do you want to cluster all sample sequences?\n"
                        f"among {alloci}?\n"
                        f"OR 'end' 'home' 'exit': ")
        while loc2cat not in alloci and loc2cat not in ["end", "home", "exit"]:
            loc2cat = input("\n-->ERROR: locus name is not valid, please enter a valid name \n"
                            f"among {alloci}: ")
        else:
            if loc2cat in ["end", "home"]:
                main()
            elif loc2cat == "exit":
                sys.stdout.write("\nSee the results of concatenation session\n"
                                 f"Results in ---> {current_dir}/outputs/Stats_option_3.txt\n"
                                 "\n\nCLUSTERING SESSION FOR PHYLOGENETIC PAIRED-END LOCI is COMPLETE\n")
                sys.stdout.write("\n\tBye!...\n\n")
                quit()
        stat_3 = open('./outputs/Stats_option_3.txt', 'a')
        os.chdir(loc2cat)
        files2cat = glob.glob('*_select.fas')
        with open(f"./{loc2cat}_allseq_select.fasta", "w") as out:
            for file in files2cat:
                with open(file, "r") as out2:
                    out.write(out2.read())
        tot = open("./" + loc2cat + "_allseq_select.fasta")
        nb_tot = tot.read().count(">")
        stat_3.writelines(f"The locus {loc2cat} has {nb_tot} sequences\n")
        sys.stdout.write(f"\nLocus {loc2cat}: {nb_tot} sequences from {nb_samples} samples have been concatenated\n"
                         f"Results in ---> {current_dir}/{loc2cat}/{loc2cat}_allseq_select.fasta\n")
        os.chdir('../')
        stat_3.close()


def prevent():
    global current_dir
    if os.path.isfile(f"{current_dir}/outputs/parameters_option_1.txt") is False:
        sys.stdout.write("\nYou have to previously run mandatory OPTION 1 before running this option \n")
        q = input("\nDo you want to run OPTION 1?, reply 'yes' 'no': ")
        if q == 'yes':
            main()
            quit()
        else:
            sys.stdout.write("\n\n")
            quit()


def main_menu1():
    os.system("cls" if winOS else "clear")
    global rmenu
    rmenu = input("\n----- BASIC ANALYSIS - only option 1 is strictly mandatory -----\n\n"
                  "1  -> NEW COMPLETE ANALYSIS (mandatory)\n"
                  "1a -> Re-analyze all loci, from the clustering step, modifying parameters\n"
                  "1b -> Re-analyze only one locus of paired-end amplicon (merged reads), modifying parameters\n"
                  "1c -> Re-analyze only one locus of single-end amplicon (R1 only), modifying parameters\n"
                  "1d -> Re-analyse only one sample, modifying parameters\n"
                  "1e -> Optional quality checking of fastq files (slow)\n\n"
                  "\nEnter '1' '1a' '1b' '1c' '1d' '1e' to run analysis\n"
                  "OR 'end' 'home' 'exit': ")
    while rmenu not in ['1', '1a', '1b', '1c', '1d', '1e', 'end', 'home', 'exit']:
        rmenu = input(f"\n-->ERROR: Please enter a CORRECT NAME of an option\n"
                      f" among '0' '1' '1a' '1b' '1c' '1d', 'home' 'exit': ")
    else:
        if rmenu in ["home", "end"]:
            main()
        elif rmenu == 'exit':
            sys.stdout.write("\n\tBye!...\n\n")
            quit()
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
    os.system("cls" if winOS else "clear")
    global rmenu
    rmenu = input("\n----- SELECTION OF MINIMUM SEQUENCE ABUNDANCES ACCORDING TO USER-DEFINED THRESHOLDS -----\n\n"
                  "2a -> Apply the SAME size threshold for ALL SAMPLES for the loci based on PAIRED-END reads "
                  "(R1/R2 merged)\n"
                  "\ti.e. you want to keep only sequences whose abundance (size)\n"
                  "\tis greater than x% of the total number of sequences for a given sample.\n"
                  "\tThis threshold of x% can be chosen for each locus.\n\n"
                  "2b -> Apply the SAME size threshold for ALL SAMPLES for the loci based on SINGLE-END reads "
                  "(R1 only)\n"
                  "\tsame as option 1 but only using the R1 reads instead of merged ones.\n\n"
                  "2c -> Apply a SPECIFIC size threshold for EACH SAMPLE, for the loci based on PAIRED-END reads "
                  "(R1/R2 merged)\n"
                  "\ti.e. you want to modulate the threshold of x% by locus but also by sample\n"
                  "\twithin a particular locus.\n\n"
                  "2d -> Apply a SPECIFIC size threshold for EACH SAMPLE, for the loci based on SINGLE-END reads "
                  "(R1 only)\n"
                  "\tsame as option 2c but only using the R1 sequences instead of merged ones.\n\n"
                  "\nEnter '2a' '2b' '2c' '2d' to run analysis\n"
                  "OR 'end' 'home' 'exit': ")
    while rmenu not in ['2a', '2b', '2c', '2d', 'end', 'home', 'exit']:
        rmenu = input(f"\n-->ERROR: Please enter a GOOD NAME of an option\n"
                      f" among '2a' '2b' '2c' '2d' 'end' 'home' 'exit': ")
    else:
        if rmenu in ["home", "end"]:
            main()
        elif rmenu == 'exit':
            sys.stdout.write("\n\tBye!...\n\n")
            quit()
        elif rmenu == "2a":
            menu2a()
        elif rmenu == "2b":
            menu2b()
        elif rmenu == "2c":
            menu2c()
        elif rmenu == "2d":
            menu2d()


def main_menu3():
    prevent()
    prev_param()
    concat_3()
    rerun()


def menu1():
    in_dir_fastq()
    in_fastqr1_user()
    in_fastqr2_user()
    in_loci1_user()
    in_loci2_user()
    in_sample_user()
    in_minsize_user()
    in_minseqlength()
    in_alpha()
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


def menu1a():
    prevent()
    prev_param()
    in_minsize_user()
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
    prevent()
    prev_param()
    in_loc_sel_merged()
    in_minsize_user()
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
    prevent()
    prev_param()
    in_loc_sel_r1()
    in_minsize_user()
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
    prevent()
    prev_param()
    in_sam_sel()
    in_minsize_user()
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
    prevent()
    prev_param()
    quality()
    runs_1x()
    rerun()


def menu2a():
    prevent()
    prev_param()
    trim_2x()
    rerun()


def menu2b():
    prevent()
    prev_param()
    trim_2x()
    rerun()


def menu2c():
    prevent()
    prev_param()
    trim_2x()
    rerun()


def menu2d():
    prevent()
    prev_param()
    trim_2x()
    rerun()


def menu3():
    prevent()
    prev_param()
    concat_3()
    rerun()


def rerun():
    global next_run
    next_run = input("Do you want continue with mbctools? enter 'yes' or 'no': ")
    if next_run == "yes":
        main()
    else:
        sys.stdout.write("\n\tBye!...\n\n")
        quit()


# Main menu
#######################################################################################################################
def main():
    os.system("cls" if winOS else "clear")

    global menu
    sys.stdout.write("-------------------- MAIN MENU --------------------\n"
                     "\nValidate without typing anything enters the default value, if any\n"
                     "Entering 'end' returns to the program upper level, if any\n"
                     "Entering 'home' returns to this main menu\n"
                     "Entering 'exit' leaves the program\n\n")
    menu = input("\n1 -> BASIC ANALYSES\n\n"
                 "2 -> SELECTION OF MINIMUM SEQUENCE ABUNDANCES ACCORDING TO USER-DEFINED THRESHOLDS\n\n"
                 "3 -> CONCATENATION OF ALL SAMPLES BY LOCUS FOR PHYLOGENETIC ANALYZES\n\n"
                 "\nEnter '1', '2', '3'\n"
                 "OR 'end' 'home' 'exit': ")
    while menu not in ['1', '2', '3', 'end', 'home', 'exit']:
        menu = input("\n-->ERROR: Please enter a CORRECT NAME of an option among '1', '2', '3'\n"
                     "OR 'end' 'home' 'exit': ")
    else:
        if menu in ['end', 'home', 'exit']:
            sys.stdout.write("\n\tBye!...\n\n")
            quit()
        elif menu == '1':
            main_menu1()
        elif menu == '2':
            main_menu2()
        elif menu == '3':
            main_menu3()
    return menu


if __name__ == "__main__":
    main()
