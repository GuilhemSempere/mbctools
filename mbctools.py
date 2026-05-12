#!/usr/bin/python3

"""mbctools a cross-platform toolkit to make use of VSEARCH easier and interactive, thus helping analyze
metabarcoding data in the best conditions. It assumes VSEARCH is pre-installed and consists in the following MAIN MENU:

1 -> INITIAL ANALYSIS (mandatory)
        merging, dereplication, clustering, chimera detection, affiliation of clusters to loci, sequence re-orientation
2 -> PRIMER REMOVAL AND SELECTION OF MINIMUM SEQUENCE ABUNDANCE LEVELS ACCORDING TO USER-DEFINED THRESHOLDS
3 -> GENERATION OF A UNIQUE SEQUENCE FILE FOR EACH LOCUS (comprising all samples' data)
4 -> CONVERTING ANALYSIS OUTPUTS FOR EXTERNAL TOOLS

Option 1 offers the following submenu:

1  -> INITIAL ANALYSIS (mandatory)
1a -> Re-analyze all loci, from the clustering step, modifying parameters
1b -> Re-analyze only one locus of paired-end amplicon (merged reads), modifying parameters
1c -> Re-analyze only one locus of single-end amplicon (R1 only), modifying parameters
1d -> Re-analyse a given sample, modifying parameters
1e -> Optional quality checking of fastq files (slow)

Option 2 offers the following submenu:

2a -> Apply the SAME size threshold for ALL SAMPLES for the loci based on PAIRED-END reads (R1/R2 merged)
2b -> Apply the SAME size threshold for ALL SAMPLES for the loci based on SINGLE-END reads (R1 only)
2c -> Apply a SPECIFIC size threshold for a given sample, for the loci based on PAIRED-END reads (R1/R2 merged)
2d -> Apply a SPECIFIC size threshold for a given sample, for the loci based on SINGLE-END reads (R1 only)

Option 4 offers the following submenu:

4a -> Generate sequence files
4b -> Generate assignment file
4c -> Build metaXplor-format sample metadata file from provided tabulated file
4d -> Build metaXplor archive
4e -> Generate MIAN taxonomy and build MIAN archive
4f -> Build Galaxy phylogeny pipeline archive


VSEARCH reference:
Rognes T, Flouri T, Nichols B, Quince C, Mahe F (2016). VSEARCH: a versatile open source tool for metagenomics
PeerJ 4:e2584 doi: 10.7717/peerj.2584 https://doi.org/10.7717/peerj.2584

Usage:
=====
        python3 mbctools.py
        (or simply "mbctools") if installed using pip
"""

__authors__ = "Christian Barnabé, Guilhem Sempéré"
__contact__ = "guilhem.sempere@cirad.fr"
__date__ = "2026-05-11"
__version__ = "2.0.0a5"
__copyright__ = "Copyright (c) 2024-2026 IRD, CIRAD"
__license__ = "This software is licensed under the MIT License. The full license text is available at https://github.com/GuilhemSempere/mbctools/blob/main/LICENSE"

import sys
import time
import datetime
from pathlib import Path
import os
import platform
import subprocess
import re
import glob
import configparser
import traceback
import io
import zipfile
import shutil
import math

from mbc_shell import start_log_redirect, end_log_redirect, main_stream_message, logFileMessage
from mbc_text import dos2unix, replaceInFile
from mbc_sequences import concat_sequences_by_locus, derep_based_on_ids, derep_several_fasta_files
import mbc_option1
import mbc_option2
import mbc_option3
import mbc_option4

winOS = "Windows" in platform.uname()
shellCmd = "powershell" if winOS else "bash"
scriptExt = "ps1" if winOS else "sh"
fileSep = "\\" if winOS else "/"
globalErrorOnStopCmd = "" if winOS else "set -e"
localErrorOnStopCmd = "; If ($LASTEXITCODE -gt 0) { exit $LASTEXITCODE }" if winOS else ""
date = datetime.datetime.now()
current_dir = os.getcwd()

errorStyle = "\033[91m"
warningStyle = "\033[93m"
normalStyle = "\033[0m"
titleStyle = "\033[94m\033[1m"
promptStyle = "\033[96m"
successStyle = "\033[92m"
citationStyle = "\033[1;33m\033[1m\033[3m"

global lociPEs, lociSEs, lociPE, lociSE, dir_fastq, fastq_R1, fastqr1s, fastq_R2, fastqr2s, Samples, \
        samples, alpha, identity, loc_sel1, loc_sel2, rmenu, sam_sel, minsize, minseqlength, loc2trim, trim_left, \
        trim_right, ts, ts1, sam2trim2c, sam2trim2d, all_loci, next_run, menu, ts2, loc2cat, loc2trim2a, loc2trim2b, \
        loc2trim2c, loc2trim2d

metaXplorFasta = "tmp_files/metaXplor_sequences.fasta"
metaXplorSequenceComposition = "tmp_files/metaXplor_sequences.tsv"
metaXplorAssignments = "tmp_files/metaXplor_assignments.tsv"
metaXplorSamples = "tmp_files/metaXplor_samples.tsv"


def promptUser(message, defaultResponse, validResponses, inputType, backFunction, exitMessage):
        """Prompts user input until a valid response is obtained. Input types are 1:text, 2:numeric, 3:file-path, 4:folder-path
        """
        firstAttempt = True
        response = None
        while response is None:
                if firstAttempt:
                        response = input("\n" + promptStyle + message + (("\n(default = " + normalStyle + defaultResponse + promptStyle + ")") if defaultResponse is not None else "") + ": " + normalStyle)
                        firstAttempt = False
                else:
                        response = input(errorStyle + "Wrong input, try again!" + promptStyle + " Accepted entries are " + normalStyle + ("a valid " + (("folder path" if inputType == 4 else ("file path" if inputType == 3 else "number")) + " or ") if inputType > 1 else "") + ", ".join(validResponses) + promptStyle + ": " + normalStyle)
                if response not in validResponses:
                        if response == "" and defaultResponse is not None:
                                response = defaultResponse
                        if inputType == 1:
                                response = None
                        elif inputType == 2 and response.isnumeric() is False:
                                response = None
                        elif inputType == 3 and Path(response).is_file() is False:
                                response = None
                        elif inputType == 4 and Path(response).is_dir() is False:
                                response = None

        if response in ["back", "home", "exit"]:
                if exitMessage is not None and exitMessage != "":
                        print(successStyle + "\n" + exitMessage + "\n" + normalStyle)
                        if response != "exit":
                                input("Press ENTER to continue ")
                if response == "back":
                        response = None
                        if backFunction is not None:
                                backFunction()
                elif response == "home":
                        response = None
                        main()
                elif response == "exit":
                        quit_mbctools()

        if response is not None:
                return response


def folders():
        """ Creates folders useful for the metabarcoding analyze
        """
        alreadyExisting = []

        global lociPEs, lociSEs, lociPE, lociSE
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
        folder = "results_by_locus"
        path = os.path.join(current_dir, folder)
        if Path(path).is_dir():
                alreadyExisting.append(folder)
        else:
                os.mkdir(path)

        sys.stdout.write("")
        for locus in list(set(lociPEs) | set(lociSEs)):
                path = os.path.join(current_dir + "/results_by_locus", locus)
                if Path(path).is_dir():
                        alreadyExisting.append("results_by_locus/" + locus)
                else:
                        os.chdir(f"{current_dir}{fileSep}results_by_locus")
                        os.mkdir(path)

        if len(alreadyExisting) > 0:
                sys.stdout.write(warningStyle + f"\nThe following folders already exist and will be used in the upcoming analysis: " + normalStyle + ", ".join(alreadyExisting))


def in_dir_fastq():
        """Input of the path containing the fastq files, option 1
        """
        global dir_fastq
        dir_fastq = promptUser("Enter the FULL PATH of the folder where fastq files are located", f"{current_dir}" + fileSep + "fastq", ["back", "home", "exit"], 4, main_menu1, "")


def in_fastq_R1():
        """Input of the file name containing the R1 fastq file names, option 1
        """
        global fastq_R1, fastqr1s
        fastq_R1 = promptUser("Enter the name of the file listing all R1 (or single-end) fastq file names", "fastqR1.txt", ["back", "home", "exit"], 3, main_menu1, "")
        with open(fastq_R1, "r") as out1:
                fastqr1s = out1.read().splitlines()


def in_fastq_R2():
        """Input of the file name containing the R2 fastq file names, option 1
        """
        global fastq_R2, fastqr2s
        fastq_R2 = promptUser("Enter the name of the file listing all R2 fastq file names (file must exist but may be empty if no R2 to process)", "fastqR2.txt", ["back", "home", "exit"], 3, main_menu1, "")
        with open(fastq_R2, "r") as out2:
                fastqr2s = out2.read().splitlines()


def in_lociPE():
        """Input of the file name containing the list of paired-end based loci, option 1
        """
        global lociPE, lociPEs
        lociPE = promptUser("Enter the name of the file listing names of loci based on paired-end reads (file may be empty if no R2 to process)", "lociPE.txt", ["back", "home", "exit"], 3, main_menu1, "")
        with open(lociPE, "r") as out:
                lociPEs = out.read().splitlines()
        for locusName in lociPEs:
                if not os.path.isfile(current_dir + "/refs/" + locusName + ".fas"):
                        print(errorStyle + "File does not exist: " + current_dir + "/refs/" + locusName + ".fas")
                        return in_lociPE()
        return lociPEs


def in_lociSE():
        """Input of the file name containing the list of single-end based loci, option 1
        """
        global lociSE, lociSEs
        lociSE = promptUser("Enter the name of the file containing loci based on R1/single-end reads (or unmerged R1)", "lociSE.txt", ["back", "home", "exit"], 3, main_menu1, "")
        with open(lociSE, "r") as out:
                lociSEs = out.read().splitlines()
        for locusName in lociSEs:
                if not os.path.isfile(current_dir + "/refs/" + locusName + ".fas"):
                        print(errorStyle + "File does not exist: " + current_dir + "/refs/" + locusName + ".fas")
                        return in_lociSE()
        return lociSEs


def in_Samples():
        """Input of the file name containing the list of samples, option 1
        """
        global Samples, samples
        Samples = promptUser("Enter the name of the file listing sample names", "samples.txt", ["back", "home", "exit"], 3, main_menu1, "")
        with open(Samples, "r") as out5:
                samples = out5.read().splitlines()


def in_minsize():
        """Input of the minimum abundance of sequences to retain for denoising/clustering, options 1, 1a, 1b, 1c and 1d
        """
        global minsize
        minsize = promptUser("Enter the minsize option value for clusters, i.e. the minimum sequence abundance of the retained clusters", "8", ["back", "home", "exit"], 2, main_menu1, "")


def in_minseqlength():
        """Input of the minimum length of sequences to keep for any locus, options 1, 1a, 1b, 1c and 1d
        """
        global minseqlength
        minseqlength = promptUser("Enter the minimum length of sequences to keep for any locus", "100", ["back", "home", "exit"], 2, main_menu1, "")


def in_alpha():
        """Input of the alpha parameter for denoising/clustering, options 1, 1a, 1b, 1c and 1d
        """
        global alpha
        alpha = promptUser("Enter alpha parameter (integer) for the clustering", "2", ["back", "home", "exit"], 2, main_menu1, "")


def in_identity():
        """Input of the identity parameter (0, 1.0) to match the clusters against reference sequences, for affiliating
        clusters to the different loci, options 1, 1a, 1b, 1c and 1d
        """
        global identity
        identity = promptUser("Enter identity parameter to match the clusters against references (as a percentage), enter an integer from 0 to 100", "70", ["back", "home", "exit"], 2, main_menu1, "")
        identity = int(identity) / 100


def in_loc_sel_merged():
        """Input of a selected locus based on paired-end reads to rerun for option 1b
        """
        global loc_sel1
        loc_sel1 = promptUser("Enter the name of the locus analysed by paired-end reads you want to rerun", None, lociPEs + ["back", "home", "exit"], 1, main_menu1, "")


def in_loc_sel_r1():
        """Input of a selected locus based on single-end read (R1) to rerun for option 1c
        """
        global loc_sel2, rmenu
        loc_sel2 = promptUser("Enter the name of the locus analysed by only single-end (R1) reads you want to rerun", None, lociSEs + ["back", "home", "exit"], 1, main_menu1, "")


def in_sam_sel():
        """Input of the sample name to rerun for option 1d
        """
        global sam_sel, samples
        sam_sel = promptUser("Enter the sample name you want to rerun", None, samples + ["back", "home", "exit"], 1, main_menu1, "")


def in_loc2trim_2x():
        """Input of the loci names for selection of minium sequence abundances according to user-defined thresholds,
        options 2a, 2b, 2c and 2d
        """
        global loc2trim2a, loc2trim2b, loc2trim2c, loc2trim2d, lociPEs, lociSEs
        if rmenu == "2a":
                loc2trim2a = promptUser("Enter a LOCUS name based on paired-end mergeable reads you want to analyze", None, lociPEs + ["back", "home", "exit"], 1, main_menu2, f"Selection of minimum sizes according to user-defined thresholds is complete.\nStatistical summary for paired-end based loci --> {current_dir}{fileSep}outputs{fileSep}Stats_option_2a.txt")
                return loc2trim2a

        if rmenu == "2b":
                loc2trim2b = promptUser("Enter a LOCUS name based on single-end R1 reads you want to analyze ", None, lociSEs + ["back", "home", "exit"], 1, main_menu2, f"Selection of minimum sizes according to user-defined thresholds is complete.\nStatistical summary for single-end based loci --> {current_dir}{fileSep}outputs{fileSep}Stats_option_2b.txt")
                return loc2trim2b

        if rmenu == "2c":
                loc2trim2c = promptUser("Enter a LOCUS name based on paired-end mergeable reads you want to analyze", None, lociPEs + ["back", "home", "exit"], 1, main_menu2, f"Selection of minimum sizes according to user-defined thresholds is complete.\nStatistical summary for paired-end based loci --> {current_dir}{fileSep}outputs{fileSep}Stats_option_2c.txt")
                return loc2trim2c

        if rmenu == "2d":
                loc2trim2d = promptUser("Enter a LOCUS name based on single-end R1 reads you want to analyze", None, lociSEs + ["back", "home", "exit"], 1, main_menu2, f"Selection of minimum sizes according to user-defined thresholds is complete.\nStatistical summary for single-end based loci --> {current_dir}{fileSep}outputs{fileSep}Stats_option_2d.txt")
                return loc2trim2d


def in_trim_sample2c(locus):
        """Input of sample names for option 2c
        """
        global samples, sam2trim2c
        sam2trim2c = promptUser("Enter the name of the sample you want to trim for locus " + locus, None, samples + ["back", "home", "exit"], 1, trim_2x, "")
        return sam2trim2c


def in_trim_sample2d(locus):
        """Input of sample names for option 2d
        """
        global samples, sam2trim2d
        sam2trim2d = promptUser("Enter the name of the sample you want to trim for locus " + locus, None, samples + ["back", "home", "exit"], 1, trim_2x, "")
        return sam2trim2d


def in_trim_left(orientFileSuffix):
        """Input of the number of bp corresponding to the forward primer to remove from the clusters,
        options 2a, 2b, 2c and 2d
        """
        global trim_left #, loc2trim2a, loc2trim2b, loc2trim2c
        if rmenu == "2a":
                trim_left = promptUser(f"Enter the number of bp of the forward primer for {loc2trim2a} (e.g. 20)", None, ["back", "home", "exit"], 2, main_menu2, "")

        if rmenu == "2b":
                if orientFileSuffix == "_plus":
                        trim_left = promptUser(f"Enter the number of bp of the forward primer for {loc2trim2b} (e.g. 20)", None, ["back", "home", "exit"], 2, main_menu2, "")
                elif orientFileSuffix == "_minus":
                        trim_left = promptUser(warningStyle + "You chose to keep antisense (-) clusters and therefore should safely be able to use 0 here\n" + promptStyle + f"Enter the number of bp of the forward primer for {loc2trim2b} (e.g. 0)", None, ["back", "home", "exit"], 2, main_menu2, "")
                else:
                        trim_left = promptUser(warningStyle + "Inspite of given advice, you chose to keep both sense (+) and antisense (-) clusters. Now you decide whether or no to trim forward primers ;-)\n" + promptStyle + f"Enter the number of bp of the forward primer for {loc2trim2b}", None, ["back", "home", "exit"], 2, main_menu2, "")

        if rmenu == "2c":
                trim_left = promptUser(f"Enter the number of bp of the forward primer for {loc2trim2c} (e.g. 20)", None, ["back", "home", "exit"], 2, main_menu2, "")

        if rmenu == "2d":
                if orientFileSuffix == "_plus":
                        trim_left = promptUser(f"Enter the number of bp of the forward primer for {loc2trim2d} (e.g. 20)", None, ["back", "home", "exit"], 2, main_menu2, "")
                elif orientFileSuffix == "_minus":
                        trim_left = promptUser(warningStyle + "You chose to keep antisense (-) clusters and therefore should safely be able to use 0 here\n" + promptStyle + f"Enter the number of bp of the forward primer for {loc2trim2d} (e.g. 0)", None, ["back", "home", "exit"], 2, main_menu2, "")
                else:
                        trim_left = promptUser(warningStyle + "Inspite of given advice, you chose to keep both sense (+) and antisense (-) clusters. Now you decide whether or no to trim forward primers ;-)\n" + promptStyle + f"Enter the number of bp of the forward primer for {loc2trim2d}", None, ["back", "home", "exit"], 2, main_menu2, "")
        return trim_left


def in_trim_right(orientFileSuffix):
        """Input of the number of bp corresponding to the reverse primer to remove from the clusters,
        options 2a, 2b, 2c and 2d
        """
        global trim_right #, loc2trim2a, loc2trim2b, loc2trim2c, loc2trim2d
        if rmenu == "2a":
                trim_right = promptUser(f"Enter the number of bp of the reverse primer for {loc2trim2a} (e.g. 22)", None, ["back", "home", "exit"], 2, main_menu2, "")

        if rmenu == "2b":
                if orientFileSuffix == "_minus":
                        trim_right = promptUser(f"Enter the number of bp of the reverse primer for {loc2trim2b} (e.g. 20)", None, ["back", "home", "exit"], 2, main_menu2, "")
                elif orientFileSuffix == "_plus":
                        trim_right = promptUser(warningStyle + "You chose to keep sense (+) sequences and therefore should safely be able to use 0 here\n" + promptStyle + f"Enter the number of bp of the reverse primer for {loc2trim2b} (e.g. 0)", None, ["back", "home", "exit"], 2, main_menu2, "")
                else:
                        trim_right = promptUser(warningStyle + "Inspite of given advice, you chose to keep both sense (+) and antisense (-) clusters. Now you decide whether or no to trim reverse primers ;-)\n" + promptStyle + f"Enter the number of bp of the reverse primer for {loc2trim2b}", None, ["back", "home", "exit"], 2, main_menu2, "")

        if rmenu == "2c":
                trim_right = promptUser(f"Enter the number of bp of the reverse primer for {loc2trim2c} (e.g. 22)", None, ["back", "home", "exit"], 2, main_menu2, "")

        if rmenu == "2d":
                if orientFileSuffix == "_minus":
                        trim_right = promptUser(f"Enter the number of bp of the reverse primer for {loc2trim2d} (e.g. 20)", None, ["back", "home", "exit"], 2, main_menu2, "")
                elif orientFileSuffix == "_plus":
                        trim_right = promptUser(warningStyle + "You chose to keep sense (+) sequences and therefore should safely be able to use 0 here\n" + promptStyle + f"Enter the number of bp of the reverse primer for {loc2trim2d} (e.g. 0)", None, ["back", "home", "exit"], 2, main_menu2, "")
                else:
                        trim_right = promptUser(warningStyle + "Inspite of given advice, you chose to keep both sense (+) and antisense (-) clusters. Now you decide whether or no to trim reverse primers ;-)\n" + promptStyle + f"Enter the number of bp of the reverse primer for {loc2trim2d}", None, ["back", "home", "exit"], 2, main_menu2, "")
        return trim_right


def in_ts():
        """Input of the user-defined threshold of minimum abundance of clusters to retain, options 2a, 2b, 2c and 2d
        """
        global ts, ts1, loc2trim2a, loc2trim2b, loc2trim2c, loc2trim2d
        if rmenu == "2a":
                ts = promptUser(f"Enter the THRESHOLD (integer between 0 and 100) you want to use for locus {loc2trim2a}.\nExample: " + normalStyle + f"if you want to keep only the clusters whose abundance is greater than 5% of the sum of cluster sizes for a given sample with {loc2trim2a}, enter 5" + promptStyle, None, ["back", "home", "exit"], 2, main_menu2, "")

        if rmenu == "2b":
                ts = promptUser(f"Enter the THRESHOLD (integer between 0 and 100) you want to use for locus {loc2trim2b}.\nExample: " + normalStyle + f"if you want to keep only the clusters whose abundance is greater than 5% of the sum of cluster sizes for a given sample with {loc2trim2b}, enter 5" + promptStyle, None, ["back", "home", "exit"], 2, main_menu2, "")

        if rmenu == "2c":
                ts = promptUser(f"Enter the THRESHOLD (integer between 0 and 100) you want to use for locus {loc2trim2c}.\nExample: " + normalStyle + f"if you want to keep only the clusters whose abundance is greater than 5% of the sum of cluster sizes for a given sample with {loc2trim2c}, enter 5" + promptStyle, None, ["back", "home", "exit"], 2, main_menu2, "")

        if rmenu == "2d":
                ts = promptUser(f"Enter the THRESHOLD (integer between 0 and 100) you want to use for locus {loc2trim2d}.\nExample: " + normalStyle + f"if you want to keep only the clusters whose abundance is greater than 5% of the sum of cluster sizes for a given sample with {loc2trim2d}, enter 5" + promptStyle, None, ["back", "home", "exit"], 2, main_menu2, "")

        ts1 = int(ts) / 100
        sys.stdout.write("\n")
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
                contents["lociPE"] = lociPE
                contents["lociSE"] = lociSE
                contents["Samples"] = Samples
                contents["alpha"] = alpha
                contents["identity"] = identity
        elif rmenu == '1b':
                contents["loc_sel1"] = lociPE
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
                global fastqr1s, fastqr2s, lociPEs, lociSEs, samples
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
                global lociPE
                lociPE = contents["lociPE"]
                global lociSE
                lociSE = contents["lociSE"]
                global Samples
                Samples = contents["Samples"]

                with open(fastq_R1, "r") as out1:
                        fastqr1s = out1.read().splitlines()

                with open(fastq_R2, "r") as out2:
                        fastqr2s = out2.read().splitlines()

                with open(lociPE, "r") as out3:
                        lociPEs = out3.read().splitlines()

                with open(lociSE, "r") as out4:
                        lociSEs = out4.read().splitlines()

                with open(Samples, "r") as out5:
                        samples = out5.read().splitlines()

                return fastqr1s, fastqr2s, lociPEs, lociSEs, samples
        except KeyError:
                print(errorStyle + "\nMissing parameter in configuration file " + fileToParse + " - " + normalStyle)
                traceback.print_exc(limit=0)
                customExit(1)


def quality():
        return mbc_option1.quality(sys.modules[__name__])


def merging():
        return mbc_option1.merging(sys.modules[__name__])



def fastq2fas():
        return mbc_option1.fastq2fas(sys.modules[__name__])


def derep_1():
        return mbc_option1.derep_1(sys.modules[__name__])


def cluster_1x():
        return mbc_option1.cluster_1x(sys.modules[__name__])


def chimera_remove():
        return mbc_option1.chimera_remove(sys.modules[__name__])


def runloc_merged():
        return mbc_option1.runloc_merged(sys.modules[__name__])


def runloc_r1():
        return mbc_option1.runloc_r1(sys.modules[__name__])


def runlocsel_merged():
        return mbc_option1.runlocsel_merged(sys.modules[__name__])


def runlocsel_r1():
        return mbc_option1.runlocsel_r1(sys.modules[__name__])


def runloc_one_sample_1d():
        return mbc_option1.runloc_one_sample_1d(sys.modules[__name__])


def orient_1x():
        return mbc_option1.orient_1x(sys.modules[__name__])


def getSingleSeqOrientFileSuffix(loci, samples, rmenu):
        return mbc_option2.get_single_seq_orient_file_suffix(sys.modules[__name__], loci, samples, rmenu)


def runs_1x():
        return mbc_option1.runs_1x(sys.modules[__name__])


def stats_1x():
        return mbc_option1.stats_1x(sys.modules[__name__])


def trim_2x():
        return mbc_option2.trim_2x(sys.modules[__name__])


def concat_3(loci):
        return mbc_option3.concat_3(sys.modules[__name__], loci)


def prevent():
        """Forces the user to run the mandatory option 1 before any other option
        """
        global current_dir, rmenu
        if os.path.isfile(f"{current_dir}{fileSep}outputs{fileSep}parameters_option_1.cfg") is False:
                sys.stdout.write("\nYou have to run mandatory OPTION 1 " + warningStyle + "before" + normalStyle + " running this option \n")
                q = promptUser("Do you want to run OPTION 1? " + normalStyle + "Reply yes or exit" + promptStyle, None, ["yes", "exit"], 1, None, "")
                if q == 'yes':
                        rmenu = '1'
                        menu1()
                quit_mbctools()


def quit_mbctools():
        print(successStyle + "\n\nThanks for using mbctools!" + normalStyle)
        printHowToCite()
        customExit(0)


def printHowToCite():
        print(warningStyle + "\nPlease cite this software as follows:" +
                        citationStyle + "\nmbctools:\tA User-Friendly Metabarcoding and Cross-Platform Pipeline for Analyzing\n\t\tMultiple Amplicon Sequencing Data across a Large Diversity of Organisms" + normalStyle
                        + "\nChristian Barnabé, Guilhem Sempéré, Vincent Manzanilla, Joel Moo Millan, Antoine Amblard-Rambert and Etienne Waleckx.\n" + citationStyle
                        + "https://github.com/GuilhemSempere/mbctools" + normalStyle + "  -  doi: 10.24072/pcjournal.501\n\n")


def main_menu1():
        return mbc_option1.main_menu1(sys.modules[__name__])


def main_menu2():
        return mbc_option2.main_menu2(sys.modules[__name__])


def main_menu3():
        return mbc_option3.main_menu3(sys.modules[__name__])


def main_menu4():
        return mbc_option4.main_menu4(sys.modules[__name__])


def menu1():
        return mbc_option1.menu1(sys.modules[__name__])


def menu1a():
        return mbc_option1.menu1a(sys.modules[__name__])


def menu1b():
        return mbc_option1.menu1b(sys.modules[__name__])


def menu1c():
        return mbc_option1.menu1c(sys.modules[__name__])


def menu1d():
        return mbc_option1.menu1d(sys.modules[__name__])


def menu1e():
        return mbc_option1.menu1e(sys.modules[__name__])


def menu2a():
        return mbc_option2.menu2a(sys.modules[__name__])


def menu2b():
        return mbc_option2.menu2b(sys.modules[__name__])


def menu2c():
        return mbc_option2.menu2c(sys.modules[__name__])


def menu2d():
        return mbc_option2.menu2d(sys.modules[__name__])


def menu3():
        return mbc_option3.menu3(sys.modules[__name__])


def derepBasedOnIDs(locusToFastaDict, outFastaName, outTsvName):
        return derep_based_on_ids(locusToFastaDict, outFastaName, outTsvName, samples, warningStyle, normalStyle, fileSep)


def derepSeveralFastaFiles(locusToFastaDict, outFastaName, outTsvName, logFile):
        return derep_several_fasta_files(locusToFastaDict, outFastaName, outTsvName, logFile, samples, warningStyle, normalStyle, fileSep)


def menu4a():
        return mbc_option4.menu4a(sys.modules[__name__])


def menu4b():
        return mbc_option4.menu4b(sys.modules[__name__])


def menu4c():
        return mbc_option4.menu4c(sys.modules[__name__])


def menu4d(invokedByUser):
        return mbc_option4.menu4d(sys.modules[__name__], invokedByUser)


def menu4e():
        return mbc_option4.menu4e(sys.modules[__name__])


def menu4f():
        return mbc_option4.menu4f(sys.modules[__name__])


def parse_date(date_str):
        return mbc_option4.parse_date(date_str)


def determineLatLon(cellArray, latLonIndex, latitudeIndex, longitudeIndex, sampleIndex):
        return mbc_option4.determine_lat_lon(sys.modules[__name__], cellArray, latLonIndex, latitudeIndex, longitudeIndex, sampleIndex)


def dmsToDecimal(dmsString):
        return mbc_option4.dms_to_decimal(dmsString)


def replaceInFile(filename, replacements):
        if os.path.exists(filename) and os.path.getsize(filename) > 0:
                with open(filename, 'r') as file:
                    content = file.read()

                for old_text, new_text in replacements:
                    content = content.replace(old_text, new_text)

                with open(filename, 'w') as file:
                    file.write(content)


def rerun(nextMenu):
        """Prompts the user to continue using mbctools or not
        """
        global next_run
        next_run = promptUser("Do you want to continue with mbctools? " + normalStyle + "Enter yes or no" + promptStyle, None, ["yes", "no"], 1, None, "")
        if next_run == "yes":
                if nextMenu is None:
                        main()
                else:
                        nextMenu()
        if next_run == "no":
                quit_mbctools()


def main():
        if winOS:
            permissionTestScript = "test_permissions." + scriptExt
            with open(permissionTestScript, "w") as out1:
                out1.write("dir\n")                    
            p = subprocess.run([shellCmd, "." + fileSep + permissionTestScript], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            os.remove(permissionTestScript)
           
            if p.returncode == 1:
                print(errorStyle + "Please adjust Execution-Policy to allow script execution" + normalStyle + ", then try again...")
                customExit(1)

        try:
                p = subprocess.run(["vsearch"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except FileNotFoundError:
                print(errorStyle + "Please install vsearch" + normalStyle + ", then try again...")
                customExit(1)

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

        sys.stdout.write(titleStyle + "\n------------------------------ mbctools v" + __version__ + " - MAIN MENU ------------------------------" + normalStyle + "\n")
        printHowToCite()
        sys.stdout.write(titleStyle + "NAVIGATION CONVENTIONS:\n" + normalStyle
                                         + "Entering '" + promptStyle + "back" + normalStyle + "' returns to the program upper level, if any\n"
                                         + "Entering '" + promptStyle + "home" + normalStyle + "' returns to this main menu\n"
                                         + "Entering '" + promptStyle + "exit" + normalStyle + "' leaves the program\n"
                                         + warningStyle + "Validating without typing anything applies the default value, if any\n" + normalStyle +"\n")
        print(titleStyle + "\nWe recommend executing procedures in the provided order:\n" + normalStyle)
        print(  "1 -> INITIAL ANALYSIS (" + warningStyle + "mandatory" + normalStyle + "): read merging, sample-level dereplication, sequence clustering,\n\tchimera detection, affiliation of sequences to loci, and sequence re-orientation\n\n"
                        "2 -> PRIMER REMOVAL AND SELECTION OF MINIMUM SEQUENCE ABUNDANCE LEVELS ACCORDING TO USER-DEFINED THRESHOLDS\n\n"
                        "3 -> GENERATION OF A UNIQUE SEQUENCE FILE FOR EACH LOCUS (comprising all samples' data)\n\n"
                        "4 -> CONVERTING ANALYSIS OUTPUTS FOR EXTERNAL TOOLS\n")

        global menu
        menu = promptUser("Please select an option among those listed above", None, ["1", "2", "3", "4", "exit"], 1, None, "")

        if menu == '1':
                main_menu1()
        elif menu == '2':
                main_menu2()
        elif menu == '3':
                main_menu3()
        elif menu == '4':
                main_menu4()


def customExit(code):
        input("\n(Press ENTER to exit) ")
        exit(code)


if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                print("\n")
