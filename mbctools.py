#!/usr/bin/python3

"""mbctools a cross-platform toolkit to make use of VSEARCH easier and interactive, thus helping analyze
metabarcoding data in the best conditions. It assumes VSEARCH is pre-installed and consists in the following MAIN MENU:

1 -> BASIC ANALYZES
2 -> REMOVAL OF PRIMERS AND SELECTION OF MINIMUM SEQUENCE ABUNDANCES ACCORDING TO USER-DEFINED THRESHOLDS
3 -> GENERATION OF A UNIQUE SEQUENCE FILE FOR EACH LOCUS (comprising all samples' data)
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
2c -> Apply a SPECIFIC size threshold for a given sample, for the loci based on PAIRED-END reads (R1/R2 merged)
2d -> Apply a SPECIFIC size threshold for a given sample, for the loci based on SINGLE-END reads (R1 only)

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

try:
	import dateutil.parser
except ModuleNotFoundError:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'python-dateutil'])
	import dateutil.parser

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


def start_log_redirect(filepath):
	"""Starts redirecting log messages
	"""
	if not winOS:
		return "exec 3>&1 4>&2 >" + filepath + " 2>&1\n"
	else:
		return "&{\n"


def end_log_redirect(filepath):
	"""Stops redirecting log messages
	"""
	if not winOS:
		return "exec 1>&3 2>&4\n"
	else:
		return "} 2> " + filepath + "\n"


def main_stream_message(message):
	"""Displays a main stream message on the console
	"""
	if not winOS:
		return "printf \"" + message + "\" >&3\n"
	else:
		return "Write-Host -NoNewline \"" + message + "\"\n"


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
	folder = "loci"
	path = os.path.join(current_dir, folder)
	if Path(path).is_dir():
		alreadyExisting.append(folder)
	else:
		os.mkdir(path)

	sys.stdout.write("")
	for locus in list(set(lociPEs) | set(lociSEs)):
		path = os.path.join(current_dir + "/loci", locus)
		if Path(path).is_dir():
			alreadyExisting.append("loci/" + locus)
		else:
			os.chdir(f"{current_dir}{fileSep}loci")
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
	lociSE = promptUser("Enter the name of the file containing loci based on single-end reads (or unmerged R1)", "lociSE.txt", ["back", "home", "exit"], 3, main_menu1, "")
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
	identity = promptUser("Enter identity parameter to match the clusters against references, i.e. the identity percentage, enter an integer from 0 to 100", "70", ["back", "home", "exit"], 2, main_menu1, "")
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


def in_trim_left():
	"""Input of the number of bp corresponding to the forward primer to remove from the clusters,
	options 2a, 2b, 2c and 2d
	"""
	global trim_left #, loc2trim2a, loc2trim2b, loc2trim2c
	if rmenu == "2a":
		trim_left = promptUser(f"Enter the number of bp of the forward primer for {loc2trim2a} (e.g. 20)", None, ["back", "home", "exit"], 2, main_menu2, "")

	if rmenu == "2b":
		trim_left = promptUser(f"Enter the number of bp of the forward primer for {loc2trim2b} (e.g. 20)", None, ["back", "home", "exit"], 2, main_menu2, "")

	if rmenu == "2c":
		trim_left = promptUser(f"Enter the number of bp of the forward primer for {loc2trim2c} (e.g. 20)", None, ["back", "home", "exit"], 2, main_menu2, "")

	if rmenu == "2d":
		trim_left = promptUser(f"Enter the number of bp of the forward primer for {loc2trim2d} (e.g. 20)", None, ["back", "home", "exit"], 2, main_menu2, "")
	return trim_left


def in_trim_right():
	"""Input of the number of bp corresponding to the reverse primer to remove from the clusters,
	options 2a, 2b, 2c and 2d
	"""
	global trim_right #, loc2trim2a, loc2trim2b, loc2trim2c, loc2trim2d
	if rmenu == "2a":
		trim_right = promptUser(f"Enter the number of bp of the reverse primer for {loc2trim2a} (e.g. 22)", None, ["back", "home", "exit"], 2, main_menu2, "")

	if rmenu == "2b":
		trim_right = promptUser(warningStyle + "You should safely be able to use 0 here as real single-end data has no reverse primer, and reverse reads will be filtered-out from unmerged R1 data\n" + promptStyle + f"Enter the number of bp of the reverse primer for {loc2trim2b} (e.g. 0)", None, ["back", "home", "exit"], 2, main_menu2, "")

	if rmenu == "2c":
		trim_right = promptUser(f"Enter the number of bp of the reverse primer for {loc2trim2c} (e.g. 22)", None, ["back", "home", "exit"], 2, main_menu2, "")

	if rmenu == "2d":
		trim_right = promptUser(warningStyle + "You should safely be able to use 0 here as real single-end data has no reverse primer, and reverse reads will be filtered-out from unmerged R1 data\n" + promptStyle + f"Enter the number of bp of the reverse primer for {loc2trim2d} (e.g. 0)", None, ["back", "home", "exit"], 2, main_menu2, "")
	return trim_right


def in_ts():
	"""Input of the user-defined threshold of minimum abundance of clusters to retain, options 2a, 2b, 2c and 2d
	"""
	global ts, ts1, loc2trim2a, loc2trim2b, loc2trim2c, loc2trim2d
	if rmenu == "2a":
		ts = promptUser(f"Enter the THRESHOLD (integer between 0 and 100) you want to use for this locus {loc2trim2a}.\nExample: " + normalStyle + f"if you want to keep only the clusters whose abundance is greater than 5% of the sum of sizes for a given sample with {loc2trim2a}, enter 5" + promptStyle, None, ["back", "home", "exit"], 2, main_menu2, "")

	if rmenu == "2b":
		ts = promptUser(f"Enter the THRESHOLD (integer between 0 and 100) you want to use for this locus {loc2trim2b}.\nExample: " + normalStyle + f"if you want to keep only the clusters whose abundance is greater than 5% of the sum of sizes for a given sample with {loc2trim2b}, enter 5" + promptStyle, None, ["back", "home", "exit"], 2, main_menu2, "")

	if rmenu == "2c":
		ts = promptUser(f"Enter the THRESHOLD (integer between 0 and 100) you want to use for this locus {loc2trim2c}.\nExample: " + normalStyle + f"if you want to keep only the clusters whose abundance is greater than 5% of the sum of sizes for a given sample with {loc2trim2c}, enter 5" + promptStyle, None, ["back", "home", "exit"], 2, main_menu2, "")

	if rmenu == "2d":
		ts = promptUser(f"Enter the THRESHOLD (integer between 0 and 100) you want to use for this locus {loc2trim2d}.\nExample: " + normalStyle + f"if you want to keep only the clusters whose abundance is greater than 5% of the sum of sizes for a given sample with {loc2trim2d}, enter 5" + promptStyle, None, ["back", "home", "exit"], 2, main_menu2, "")

	ts1 = int(ts) / 100
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
			out1.write(f"vsearch --fastq_qmax 93 --fastq_eestats2 \"{dir_fastq}{fileSep}{fastqr1}\" --output "
					   f"../outputs/{sample}_R1_quality.txt" + localErrorOnStopCmd + "\n")
		out1.write(main_stream_message(f'\n\n'))

	if len(lociPEs) > 0:
		with open("scripts/infor2." + scriptExt, "w") as out2:
			i = 0
			out2.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Quality statistical tests on R2 reads '
																		 f'for samples:\n'))
			while i < len(samples):
				sample = samples[i]
				fastqr2 = fastqr2s[i]
				i = i + 1
				out2.write(main_stream_message(f' {sample}...'))
				out2.write(f"vsearch --fastq_qmax 93 --fastq_eestats2 \"{dir_fastq}{fileSep}{fastqr2}\" --output "
						   f"../outputs/{sample}_R2_quality.txt" + localErrorOnStopCmd + "\n")
			out2.write(main_stream_message(f'\n\n'))


def merging():
	"""Merges paired-end reads into one sequence, when the length of the expected amplicon allows it (option 1)
	according the VSEARCH command:

	vsearch --fastq_mergepairs dir_fastq/fastqR1 --reverse dir_fastq/fastqR2 --fastaout
	../tmp_files/SN_pairedEnd.fa --fastq_allowmergestagger --relabel sample=SN_merged.

	dir_fastq =  directory containing the fastq files
	fastqR1 = complete file name of fastq R1 read
	fastqR2 = fastqR2 complete name
	option --fastq_allowmergestagger allows the merging of short fragments
	"""
	with open("scripts/merging." + scriptExt, "w") as out:
		i = 0
		out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Merging paired-end reads for samples:\n'))
		while i < len(samples):
			sample = samples[i]
			fastqr1 = fastqr1s[i]
			fastqr2 = fastqr2s[i]
			out.write(main_stream_message(
				f" {sample}...") +
					  f"vsearch --fastq_mergepairs \"{dir_fastq}{fileSep}{fastqr1}\" --reverse \"{dir_fastq}{fileSep}{fastqr2}\" "
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
		out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Converting FASTQ files into FASTA format for samples:\n'))
		while i < len(samples):
			sample = samples[i]
			fastqr1 = fastqr1s[i]
			out.write(main_stream_message(f' {sample}...') +
					  f"vsearch --fastq_qmax 93 --fastq_filter \"{dir_fastq}{fileSep}{fastqr1}\" --fastaout ../tmp_files/{sample}_singleEnd.fa --relabel sample={sample}_R1." + localErrorOnStopCmd + "\n")
			i = i + 1
		out.write(main_stream_message(f'\n\n'))


def derep_1():
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
	if len(lociPEs) > 0:
		with open("scripts/derep." + scriptExt, "w") as out:
			out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Dereplicating merged reads for samples:\n'))
			for sample in samples:
				out.write(main_stream_message(
					f' {sample}...') +
						  f"vsearch --fastx_uniques ../tmp_files/{sample}_pairedEnd.fa --fastaout "
						  f"../tmp_files/{sample}_pairedEnd_derep.fas --sizeout --strand both" + localErrorOnStopCmd + "\n")
			out.write(main_stream_message(f'\n\n'))

	if len(lociSEs) > 0:
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
	if len(lociPEs) > 0 and  rmenu in ["1", "1a", "1b"] and len(lociPEs) > 0:
		with open("scripts/cluster." + scriptExt, "w") as out:
			i = 0
			if rmenu != "1":
				print()
			out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Clustering merged reads for all samples:\n'))
			while i < len(samples):
				sample = samples[i]
				out.write(main_stream_message(
					f' {sample}...') +
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
			out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Clustering single-end or unmerged R1 reads for all samples:\n'))
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

	if rmenu == "1d":
		with open("scripts/cluster_one_sample_1d." + scriptExt, "w") as out:
			out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'\nClustering reads for selected sample {sam_sel}...'))
			if len(lociPEs) > 0:
				out.write(f"vsearch --cluster_unoise ../tmp_files/{sam_sel}_pairedEnd_derep.fas --sizein --centroids "
				  f"../tmp_files/{sam_sel}_pairedEnd_cluster.fas --strand both --minsize {minsize} --sizeout "
				  f"--unoise_alph {alpha} --minseqlength {minseqlength}" + localErrorOnStopCmd + "\n")

			if len(lociSEs) > 0:
				out.write(f"vsearch --cluster_unoise ../tmp_files/{sam_sel}_singleEnd_derep.fas --sizein --centroids "
				  f"../tmp_files/{sam_sel}_singleEnd_cluster.fas --strand both --minsize {minsize} "
				  f"--sizeout --unoise_alph {alpha} --minseqlength {minseqlength}" + localErrorOnStopCmd + "\n")
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
	
	if len(lociPEs) > 0 and rmenu in ["1", "1a", "1b"]:
		with open("scripts/chimera." + scriptExt, "w") as out:
			out.write(globalErrorOnStopCmd +
					  "\n" + main_stream_message(f'Detecting and removing chimeras within merged reads of samples:\n'))
			for sample in samples:
				out.write(main_stream_message(
					f' {sample}...') +
						  f"vsearch --uchime3_denovo ../tmp_files/{sample}_pairedEnd_cluster.fas --nonchimeras "
						  f"../tmp_files/{sample}_pairedEnd_cluster_OK.fas" + localErrorOnStopCmd + "\n")
			out.write(main_stream_message(f'\n\n'))

	if len(lociSEs) > 0 and rmenu in ["1", "1a", "1c"]:
		with open("scripts/chimera_r1." + scriptExt, "w") as out:
			out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Detecting and removing chimeras within single-end or unmerged R1 reads of samples:\n'))
			for sample in samples:
				out.write(main_stream_message(
					f' {sample}...') +
						  f"vsearch --uchime3_denovo ../tmp_files/{sample}_singleEnd_cluster.fas --nonchimeras"
						  f" ../tmp_files/{sample}_singleEnd_cluster_OK.fas" + localErrorOnStopCmd + "\n")
			out.write(main_stream_message(f'\n\n'))

	if rmenu == "1d":
		with open("scripts/chimera_one_sample_1d." + scriptExt, "w") as out:
			out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Detecting and removing chimeras for selected sample {sam_sel}...'))
			if len(lociPEs) > 0:
				out.write(f"vsearch --uchime3_denovo ../tmp_files/{sam_sel}_pairedEnd_cluster.fas --nonchimeras "
							f"../tmp_files/{sam_sel}_pairedEnd_cluster_OK.fas" + localErrorOnStopCmd + "\n")
			if len(lociSEs) > 0:
				out.write(f"vsearch --uchime3_denovo ../tmp_files/{sam_sel}_singleEnd_cluster.fas --nonchimeras "
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
		for lociPEb in lociPEs:
			for sample in samples:
				out.write(main_stream_message(
					f' {sample} vs {lociPEb}...') +
						  f"vsearch --usearch_global ../tmp_files/{sample}_pairedEnd_cluster_OK.fas --db ../refs/{lociPEb}.fas"
						  f" --matched ../loci/{lociPEb}/{sample}_pairedEnd.fas --id {identity} --strand both"
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
		for locusSEb in lociSEs:
			for sample in samples:
				out.write(main_stream_message(f' {sample} vs {locusSEb}...') +
						  f"vsearch --usearch_global ../tmp_files/{sample}_singleEnd_cluster_OK.fas --db "
						  f"../refs/{locusSEb}.fas --matched ../loci/{locusSEb}/{sample}_singleEnd.fas --id {identity} "
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
																	f'locus {loc_sel1} for samples:\n'))
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
																	f'locus {loc_sel2} for samples:\n'))
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
		for lociPEb in lociPEs:
			out.write(main_stream_message(f' {lociPEb}...') +
					  f"vsearch --usearch_global ../tmp_files/{sam_sel}_pairedEnd_cluster_OK.fas --db "
					  f"../refs/{lociPEb}.fas --matched ../loci/{lociPEb}/{sam_sel}_pairedEnd.fas --id {identity} "
					  f"--strand both" + localErrorOnStopCmd + "\n")
		out.write(main_stream_message(f'\n\n'))

	with open("scripts/loci_R1_1d." + scriptExt, "w") as out1:
		out1.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Affiliating clusters fo selected '
																	 f'sample {sam_sel} to the loci:\n'))
		for locusSEb in lociSEs:
			out1.write(main_stream_message(f' {locusSEb}...') +
					   f"vsearch --usearch_global ../tmp_files/{sam_sel}_singleEnd_cluster_OK.fas --db "
					   f"../refs/{locusSEb}.fas --matched ../loci/{locusSEb}/{sam_sel}_singleEnd.fas "
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
		if len(lociPEs) > 0:
			with open("scripts/orient_merged_1a." + scriptExt, "w") as out:
				out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f"Orienting all merged reads according to the reference:\n"))
				for locusPEb in lociPEs:
					for sample in samples:
						out.write(main_stream_message(f' {sample} vs {locusPEb}...') +
								  f"vsearch --orient ../loci/{locusPEb}/{sample}_pairedEnd.fas --db ../refs/{locusPEb}.fas "
								  f"--fastaout ../loci/{locusPEb}/{sample}_pairedEnd_orient.fas" + localErrorOnStopCmd + "\n")
				out.write(main_stream_message(f'\n\n'))

		if len(lociSEs) > 0:
			with open("scripts/orient_R1_1a." + scriptExt, "w") as out1:
				out1.write(globalErrorOnStopCmd + "\n" + main_stream_message(f"Orienting all R1/single-end reads according to the reference:\n"))
				for locusSEb in lociSEs:
					for sample in samples:
						out1.write(main_stream_message(f' {sample} vs {locusSEb}...') +
								   f"vsearch --orient ../loci/{locusSEb}/{sample}_singleEnd.fas --db ../refs/{locusSEb}.fas "
								   f"--fastaout ../loci/{locusSEb}/{sample}_singleEnd_orient.fas --tabbedout ../loci/{locusSEb}/{sample}_singleEnd_orient.tsv" + localErrorOnStopCmd + "\n")
				out1.write(main_stream_message(f'\n\n'))

	elif rmenu == "1b":
		with open("scripts/orient_merged_1b." + scriptExt, "w") as out:
			out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Orienting all merged reads according to the reference for selected locus {loc_sel1} and all samples:\n'))
			for sample in samples:
				out.write(main_stream_message(f' {sample}...') +
						  f"vsearch --orient ../loci/{loc_sel1}/{sample}_pairedEnd.fas --db "
						  f"../refs/{loc_sel1}.fas --fastaout ../loci/{loc_sel1}/{sample}_pairedEnd_orient.fas" + localErrorOnStopCmd + "\n")
			out.write(main_stream_message(f'\n\n'))

	elif rmenu == "1c":
		with open("scripts/orient_R1_1c." + scriptExt, "w") as out:
			out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Orienting all R1/single-end reads according to the reference for selected locus {loc_sel2} and all samples:\n'))
			for sample in samples:
				out.write(main_stream_message(f' {sample}...') +
						  f"vsearch --orient ../loci/{loc_sel2}/{sample}_singleEnd.fas --db ../refs/{loc_sel2}.fas "
						  f"--fastaout ../loci/{loc_sel2}/{sample}_singleEnd_orient.fas --tabbedout ../loci/{loc_sel2}/{sample}_singleEnd_orient.tsv" + localErrorOnStopCmd + "\n")
			out.write(main_stream_message(f'\n\n'))

	if rmenu == "1d":
		if len(lociPEs) > 0:
			with open("scripts/orient_merged_1d." + scriptExt, "w") as out:
				out.write(globalErrorOnStopCmd + "\n" + main_stream_message(f'Orienting all merged reads of selected sample {sam_sel} for all R1/R2 loci:\n'))
				for locusPEb in lociPEs:
					out.write(main_stream_message(f' {locusPEb}...') +
							  f"vsearch --orient ../loci/{locusPEb}/{sam_sel}_pairedEnd.fas --db ../refs/{locusPEb}.fas --fastaout "
							  f"../loci/{locusPEb}/{sam_sel}_pairedEnd_orient.fas" + localErrorOnStopCmd + "\n")
				out.write(main_stream_message(f'\n\n'))

		if len(lociSEs) > 0:
			with open("scripts/orient_R1_1d." + scriptExt, "w") as out18:
				out18.write(main_stream_message(f'Orienting all R1/single-end reads of selected sample {sam_sel} for all R1/single-end loci:\n'))
				for locusSEb in lociSEs:
					out18.write(main_stream_message(f' {locusSEb}...') +
								f"vsearch --orient ../loci/{locusSEb}/{sam_sel}_singleEnd.fas --db ../refs/{locusSEb}.fas "
								f"--fastaout ../loci/{locusSEb}/{sam_sel}_singleEnd_orient.fas --tabbedout ../loci/{locusSEb}/{sam_sel}_singleEnd_orient.tsv" + localErrorOnStopCmd + "\n")
				out18.write(main_stream_message(f'\n\n'))



def removeSingleEndReorientedSequences(loci, samples, rmenu):
	for locus in loci:
		for sample in samples:
			seqsToKeep = []
			with open(f"./loci/{locus}/{sample}_singleEnd_orient.tsv", "r") as orientInfo:
				for line in orientInfo.read().split("\n"):
					splitLine = line.split("\t")
					if len(splitLine) == 4 and splitLine[1] == '+':
						seqsToKeep.append(splitLine[0])
			print(warningStyle + f"Locus {locus} / sample {sample}: keeping only R1/single-end sequences that were initially aligned like the reference (" + str(len(seqsToKeep)) + ")\n" + normalStyle)
			with open(f"./loci/{locus}/{sample}_singleEnd_orient.tsv", "w") as orientInfo:
				for seq in seqsToKeep:
					orientInfo.write(seq + "\n")
			logFile = open(f'{current_dir}{fileSep}outputs{fileSep}res' + rmenu + '.log', 'w')
			subprocess.run(["vsearch", "--fastx_getseqs", f"./loci/{locus}/{sample}_singleEnd_orient.fas", "--labels", f"./loci/{locus}/{sample}_singleEnd_orient.tsv", "--fastaout", f"./loci/{locus}/{sample}_singleEnd_orient.fas"], stderr=logFile)
			logFile.close()


def runs_1x():
	"""Creation of global scripts for the options 1, 1a, 1b, 1c, Ad and 1e
	"""
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
							+ (("'./locimerged." + scriptExt + "' ") if len(lociPEs) > 0 else "")
							+ (("'./locir1." + scriptExt + "' ") if len(lociSEs) > 0 else "")
							+ (("'./orient_merged_1a." + scriptExt + "' ") if len(lociPEs) > 0 else "")
							+ (("'./orient_R1_1a." + scriptExt + "' ") if len(lociSEs) > 0 else "")
							+ ')\nfor script in "${scriptArray[@]}"\n' +
							  '	do\n' +
							  '	if ! ${script}; then\n' +
							  '		printf "Error executing ${script}\\n\\n" >&3\n' +
							  '		exit 1\n' +
							  '	fi\n' +
							  'done\n' +
						  end_log_redirect('../outputs/res1.log'))
			else:
				out.write(start_log_redirect('../outputs/res1.log') +
							'   $scriptArray = @(' + (('"./merging.' + scriptExt + '", ') if len(lociPEs) > 0 else "")
							+ '"./fqtofas.' + scriptExt + '", '
							+ (('"./derep.' + scriptExt + '", ') if len(lociPEs) > 0 else "")
							+ (('"./derep_r1.' + scriptExt + '", ') if len(lociSEs) > 0 else "")
							+ (('"./cluster.' + scriptExt + '", ') if len(lociPEs) > 0 else "")
							+ (('"./cluster_r1.' + scriptExt + '", ') if len(lociSEs) > 0 else "")
							+ (('"./chimera.' + scriptExt + '", ') if len(lociPEs) > 0 else "")
							+ (('"./chimera_r1.' + scriptExt + '", ') if len(lociSEs) > 0 else "")
							+ (('"./locimerged.' + scriptExt + '", ') if len(lociPEs) > 0 else "")
							+ (('"./locir1.' + scriptExt + '", ') if len(lociSEs) > 0 else "")
							+ (('"./orient_merged_1a.' + scriptExt + '"' + (', ' if len(lociSEs) > 0 else "")) if len(lociPEs) > 0 else "")
							+ (('"./orient_R1_1a.' + scriptExt + '"') if len(lociSEs) > 0 else "")
							+ ')\n   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
							'		$script = $scriptArray[$i]\n' +
							'		& "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; '
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
		removeSingleEndReorientedSequences(lociSEs, samples, rmenu)
		print(successStyle + f"Step {rmenu} ended successfully" + normalStyle)

	if rmenu == "1a":
		with open(f"scripts/runall1a.{scriptExt}", "w") as out:
			if not winOS:
				out.write(start_log_redirect('../outputs/res1a.log') + "scriptArray=("
													+ (("'./cluster." + scriptExt + "' ") if len(lociPEs) > 0 else "")
													+ (("'./cluster_r1." + scriptExt + "' ") if len(lociSEs) > 0 else "")
													+ (("'./chimera." + scriptExt + "' ") if len(lociPEs) > 0 else "")
													+ (("'./chimera_r1." + scriptExt + "' ") if len(lociSEs) > 0 else "")
													+ (("'./locimerged." + scriptExt + "' ") if len(lociPEs) > 0 else "")
													+ (("'./locir1." + scriptExt + "' ") if len(lociSEs) > 0 else "")
													+ (("'./orient_merged_1a." + scriptExt + "' ") if len(lociPEs) > 0 else "")
													+ (("'./orient_R1_1a." + scriptExt + "' ") if len(lociSEs) > 0 else "")
													+ ')\nfor script in "${scriptArray[@]}"\n' +
													  '	do\n' +
													  '	if ! ${script}; then\n' +
													  '		printf "Error executing ${script}\\n\\n" >&3\n' +
													  '		exit 1\n' +
													  '	fi\n' +
													  'done\n' +
						  end_log_redirect('../outputs/res1a.log'))
			else:
				out.write(start_log_redirect('../outputs/res1a.log') +
						  '   $scriptArray = @('
							+ (('"./cluster.' + scriptExt + '", ') if len(lociPEs) > 0 else "")
							+ (('"./cluster_r1.' + scriptExt + '", ') if len(lociSEs) > 0 else "")
							+ (('"./chimera.' + scriptExt + '", ') if len(lociPEs) > 0 else "")
							+ (('"./chimera_r1.' + scriptExt + '", ') if len(lociSEs) > 0 else "")
							+ (('"./locimerged.' + scriptExt + '", ') if len(lociPEs) > 0 else "")
							+ (('"./locir1.' + scriptExt + '", ') if len(lociSEs) > 0 else "")
							+ (('"./orient_merged_1a.' + scriptExt + '"' + (', ' if len(lociSEs) > 0 else "")) if len(lociPEs) > 0 else "")
							+ (('"./orient_R1_1a.' + scriptExt + '"') if len(lociSEs) > 0 else "")
							+ ')\n   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
							  '		$script = $scriptArray[$i]\n' +
							  '		& "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; '
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
		removeSingleEndReorientedSequences(lociSEs, samples, rmenu)
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
													  '	do\n' +
													  '	if ! ${script}; then\n' +
													  '		printf "Error executing ${script}\\n\\n" >&3\n' +
													  '		exit 1\n' +
													  '	fi\n' +
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
													  '		$script = $scriptArray[$i]\n' +
													  '		& "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; '
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
																	   "./locir1_sel." + scriptExt + "' '"
																	   "./orient_R1_1c." + scriptExt + "')\n" +
																		  'for script in "${scriptArray[@]}"\n' +
																		  '	do\n' +
																		  '	if ! ${script}; then\n' +
																		  '		printf "Error executing ${script}\\n\\n" >&3\n' +
																		  '		exit 1\n' +
																		  '	fi\n' +
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
														  '		$script = $scriptArray[$i]\n' +
														  '		& "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; '
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
		removeSingleEndReorientedSequences(lociSEs, samples, rmenu)
		print(successStyle + f"Step {rmenu} ended successfully" + normalStyle)

	if rmenu == "1d":
		with open(f"scripts/runall1d.{scriptExt}", "w") as out:
			if not winOS:
				out.write(start_log_redirect('../outputs/res1d.log') +
						  "scriptArray=('./cluster_one_sample_1d." + scriptExt + "' './chimera_one_sample_1d." + scriptExt + "' "
																+ (("'./loci_merged_1d." + scriptExt + "' ") if len(lociPEs) > 0 else "")
																+ (("'./loci_R1_1d." + scriptExt + "' ") if len(lociSEs) > 0 else "")
																+ (("'./orient_merged_1d." + scriptExt + "' ") if len(lociPEs) > 0 else "")
																+ (("'./orient_R1_1d." + scriptExt + "' ") if len(lociSEs) > 0 else "")
																+	 ')\nfor script in "${scriptArray[@]}"\n' +
																	  '	do\n' +
																	  '	if ! ${script}; then\n' +
																	  '		printf "Error executing ${script}\\n\\n" >&3\n' +
																	  '		exit 1\n' +
																	  '	fi\n' +
																	  'done\n' +
						  end_log_redirect('../outputs/res1d.log'))
			else:
				out.write(start_log_redirect('../outputs/res1d.log') +
						  '   $scriptArray = @('
							+ '"./cluster_one_sample_1d.' + scriptExt + '", "./chimera_one_sample_1d.' + scriptExt + '", '
							+ (('"./loci_merged_1d.' + scriptExt + '", ') if len(lociPEs) > 0 else "")
							+ (('"./loci_R1_1d.' + scriptExt + '", ') if len(lociSEs) > 0 else "")
							+ (('"./orient_merged_1d.' + scriptExt + '"' + (', ' if len(lociSEs) > 0 else "")) if len(lociPEs) > 0 else "")
							+ (('"./orient_R1_1d.' + scriptExt + '"') if len(lociSEs) > 0 else "")
							+ ')\n   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
							  '		$script = $scriptArray[$i]\n' +
							  '		& "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; '
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
		removeSingleEndReorientedSequences(lociSEs, [sam_sel], rmenu)
		print(successStyle + f"Step {rmenu} ended successfully" + normalStyle)

	if rmenu == "1e":
		with open("scripts/runall1e." + scriptExt, "w") as out:
			if not winOS:
				out.write(start_log_redirect('../outputs/res1e.log') +
						  "scriptArray=('"
						  "./infor1." + scriptExt + "' '"
													"./infor2." + scriptExt + "')\n" +
						  'for script in "${scriptArray[@]}"\n' +
						  '	do\n' +
						  '	if ! ${script}; then\n' +
						  '		printf "Error executing ${script}\\n\\n" >&3\n' +
						  '		exit 1\n' +
						  '	fi\n' +
						  'done\n' + end_log_redirect('../outputs/res1e.log'))
			else:
				out.write(start_log_redirect('../outputs/res1e.log') +
						  '   $scriptArray = @("'
						  './infor1.' + scriptExt + '", "'
													'./infor2.' + scriptExt + '")\n' +
						  '   For ($i=0; $i -lt $scriptArray.Length; $i++) {\n' +
						  '		$script = $scriptArray[$i]\n' +
						  '		& "$script" ; If ($LASTEXITCODE -gt 0) { "Error executing $script"; exit $LASTEXITCODE }\n' +
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


def stats_1x():
	"""Calculates the number of resulting sequences (reads, merged, dereplicates and clusters) according
	to the selected options 'minsize', 'minseqlength', 'alpha parameter' and 'identity', options 1, 1a, 1b, 1c and 1d
	"""
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

				# Nb READS Calculated on singleEnd.fa
				r1fa = open(f"tmp_files/{sample}_singleEnd.fa", "rt")
				reads = r1fa.read()
				nb_reads = reads.count(">")

				try:
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
				except FileNotFoundError:
					nb_derpr1 = 0
					nb_clusr1 = 0
					nb_clusr1ok = 0

				try:
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
					os.chdir(f"loci/{locusPE}")
					refs = open(f"{sample}_pairedEnd.fas")
					refsloc = refs.read()
					nb_ref = refsloc.count("sample")
					out.writelines(f"\t{nb_ref} clusters of merged sequences of {sample} affiliated "
								   f"to locus {locusPE}\n")
					os.chdir(current_dir)
				for locusSE in lociSEs:
					os.chdir(f"loci/{locusSE}")
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
			os.chdir(f"./loci/{loc_sel1}")
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
			os.chdir(f"./loci/{loc_sel2}")
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
				os.chdir(f"./loci/{locusPE}")
				a = open(sam_sel + "_pairedEnd.fas")
				b = a.read()
				c = b.count("sample")
				out.write(f"\t{sam_sel} has {c} clusters for locus {locusPE} based on paired-end reads\n")
				os.chdir(current_dir)

			for locusSE in lociSEs:
				os.chdir(f"./loci/{locusSE}")
				a = open(sam_sel + "_singleEnd.fas")
				b = a.read()
				c = b.count("sample")
				out.write(f"\t{sam_sel} has {c} clusters for locus {locusSE} based on single-end reads\n")
				os.chdir(current_dir)
		print(successStyle + "\nExecution of option 1d is complete\n" + normalStyle)
		out.close()


def trim_2x():
	"""Removes primers and selects clusters according to minimum abundance thresholds fixed by the user
	"""
	global loc2trim2a, loc2trim2b, loc2trim2c, loc2trim2d, ts, ts1, trim_left, trim_right, sam2trim2c, sam2trim2d
	if rmenu == "2a":
		os.chdir(current_dir)
		stat_2a = open(f"{current_dir}{fileSep}outputs{fileSep}Stats_option_2a.txt", "w")
		while True:
			loc2trim2a = in_loc2trim_2x()
			os.chdir(f"./loci/{loc2trim2a}")
			trim_left = in_trim_left()
			trim_right = in_trim_right()
			ts = in_ts()
			stat_2a.write(f"Locus {loc2trim2a} trimmed {trim_left} bp (left) and {trim_right} bp (right) with threshold set at {ts1}\n")
			for sample in samples:
				with open(f"{sample}_pairedEnd_orient.fas", 'r') as filin, open("trim-select." + scriptExt, "w") as out:
					targets = [line for line in filin if "size" in line]
					a = 0
					for target in targets:
						size = re.search('size=(.+?)$', target).group(1)
						a = a + int(size)
					b = int(a * float(ts1) + 1)
					stat_2a.writelines(f"\tSum of sizes for {sample} = {a}\n"
									   f"\tThe sizes > {b} for {sample} were retained\n")
					out.write(start_log_redirect('./' + sample + '.log') +
							  f"vsearch --fastx_filter {sample}_pairedEnd_orient.fas --fastq_stripleft {trim_left} "
							  f" --fastq_stripright {trim_right} --fastaout tmp --minsize {b}"
							  + localErrorOnStopCmd + "\n"
							  f"vsearch --derep_fulllength ./tmp --output {sample}_pairedEnd_select.fas "
							  f'--sizein --sizeout\n' + localErrorOnStopCmd + "\n"
							  + end_log_redirect('./' + sample + '.log'))
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
			stat_2a.flush()
		stat_2a.close()	

	if rmenu == "2b":
		os.chdir(current_dir)
		stat_2b = open(f"{current_dir}{fileSep}outputs{fileSep}Stats_option_2b.txt", "w")
		while True:
			loc2trim2b = in_loc2trim_2x()
			os.chdir(f"./loci/{loc2trim2b}")
			trim_left = in_trim_left()
			trim_right = in_trim_right()
			ts = in_ts()
			stat_2b.write(f"Locus {loc2trim2b} trimmed {trim_left} bp (left) and {trim_right} bp (right) with threshold set at {ts1}\n")
			for sample in samples:
				with open(sample + "_singleEnd_orient.fas", 'r') as filin, open("trim-select." + scriptExt, "w") as out:
					targets = [line for line in filin if "size" in line]
					a = 0
					for target in targets:
						size = re.search('size=(.+?)$', target).group(1)
						a = a + int(size)
					b = int(a * float(ts1) + 1)
					stat_2b.writelines(f"\tSum of sizes for {sample} = {a}\n"
									   f"\tThe sizes > {b} for {sample} were retained\n")
					out.write(start_log_redirect('./' + sample + '.log') +
							  f"vsearch --fastx_filter {sample}_singleEnd_orient.fas --fastq_stripleft {trim_left} "
							  f" --fastq_stripright {trim_right} --fastaout tmp --minsize {b}"
							  + localErrorOnStopCmd + "\n"
							  f"vsearch --derep_fulllength ./tmp --output {sample}_singleEnd_select.fas "
							  f'--sizein --sizeout\n' + localErrorOnStopCmd + "\n"
							  + end_log_redirect('./' + sample + '.log'))
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
			stat_2b.flush()
		stat_2b.close()			

	if rmenu == "2c":
		os.chdir(current_dir)
		stat_2c = open(f"{current_dir}{fileSep}outputs{fileSep}Stats_option_2c.txt", "a")
		while True:
			loc2trim2c = in_loc2trim_2x()
			os.chdir(f"{current_dir}{fileSep}loci{fileSep}{loc2trim2c}")
			trim_left = in_trim_left()
			trim_right = in_trim_right()
			stat_2c.write(f"Locus {loc2trim2c} trimmed at {trim_left} bp (left) and {trim_right} bp (right)\n")
			while True:
				sam2trim2c = in_trim_sample2c(loc2trim2c)
				ts = in_ts()
				with open(sam2trim2c + "_pairedEnd_orient.fas", "r") as filin, open("trim-select." + scriptExt, "w") as filout:
					targets = [line for line in filin if "size" in line]
					a = 0
					for target in targets:
						size = re.search('size=(.+?)$', target).group(1)
						a = a + int(size)
					b = int(a * float(ts1) + 1)
					stat_2c.write(f"\tSum of sizes for {sam2trim2c} at locus {loc2trim2c} = {a}\n"
								  f"\tAt threshold {ts} sizes > {b} for {sam2trim2c} were retained\n")
					filout.writelines(start_log_redirect('./' + loc2trim2c + '.log') +
									  f' vsearch --fastx_filter {sam2trim2c}_pairedEnd_orient.fas --fastq_stripleft {trim_left} '
									  f' --fastq_stripright {trim_right} --fastaout tmp --minsize {b}\n'
									  + localErrorOnStopCmd + "\n"
									  f' vsearch --derep_fulllength ./tmp --output {sam2trim2c}_pairedEnd_select.fas '
									  f'--sizein --sizeout\n' + localErrorOnStopCmd + "\n"
									  + end_log_redirect('./' + sam2trim2c + '.log'))
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
				stat_2c.flush()
		stat_2c.close()

	if rmenu == "2d":
		os.chdir(current_dir)
		stat_2d = open(f"{current_dir}{fileSep}outputs{fileSep}Stats_option_2d.txt", "a")
		while True:
			loc2trim2d = in_loc2trim_2x()
			os.chdir(f"{current_dir}{fileSep}loci{fileSep}{loc2trim2d}")
			trim_left = in_trim_left()
			trim_right = in_trim_right()
			stat_2d.write(f"Locus {loc2trim2d} trimmed at {trim_left} bp (left) and {trim_right} bp (right)\n")
			while True:
				sam2trim2d = in_trim_sample2d(loc2trim2d)
				ts = in_ts()
				with open(sam2trim2d + "_singleEnd_orient.fas", "r") as filin, open("trim-select." + scriptExt, "w") as filout:
					targets = [line for line in filin if "size" in line]
					a = 0
					for target in targets:
						size = re.search('size=(.+?)$', target).group(1)
						a = a + int(size)
					b = int(a * float(ts1) + 1)
					stat_2d.write(f"\tSum of sizes for {sam2trim2d} at locus {loc2trim2d} = {a}\n"
								  f"\tAt threshold {ts} sizes > {b} for {sam2trim2d} were retained\n")
					filout.writelines(start_log_redirect('./' + loc2trim2d + '.log') +
									  f' vsearch --fastx_filter {sam2trim2d}_singleEnd_orient.fas --fastq_stripleft {trim_left} '
									  f' --fastq_stripright {trim_right} --fastaout tmp --minsize {b}\n'
									  + localErrorOnStopCmd + "\n"
									  f' vsearch --derep_fulllength ./tmp --output {sam2trim2d}_singleEnd_select.fas '
									  f'--sizein --sizeout\n' + localErrorOnStopCmd + "\n"
									  + end_log_redirect('./' + sam2trim2d + '.log'))
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
				stat_2d.flush()
		stat_2d.close()


def concat_3(loci):
	"""Compiling of all sample sequences by locus in a unique file (option 3)
	"""
	global all_loci, samples, loc2cat
	nb_samples = len(samples)
	all_loci = lociPEs + list(set(lociSEs) - set(lociPEs))
	if os.path.exists("outputs/Stats_option_3.txt"):
		os.remove("outputs/Stats_option_3.txt")

	skippedLoci = []
	lociSelectedByUser = []
	locusIndex = 0
	anySeqsFound = False
	while loci is None or locusIndex < len(loci):
		loc2cat = promptUser("For which LOCUS do you want to generate a unique sequence file?", None, all_loci + ["back", "home", "exit"], 1, main, f"Results of concatenation session --> {current_dir}{fileSep}outputs{fileSep}Stats_option_3.txt") if loci is None else loci[locusIndex]

		stat_3 = open(f'{current_dir}{fileSep}outputs{fileSep}Stats_option_3.txt', 'a')
		os.chdir(f"{current_dir}{fileSep}loci{fileSep}{loc2cat}")
		files2cat = glob.glob('*_select.fas')
		if len(files2cat) == 0:
			print(errorStyle + f"\nSequences for locus {loc2cat} have not been filtered using an abundance threshold." + warningStyle + " Please run step 2 on all loci for which you want to run step 3" + normalStyle)
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
			sys.stdout.write(successStyle + f"\nLocus {loc2cat}: {nb_tot} sequences from {nb_samples} samples were added to a unique fasta file\n" + normalStyle)
			if nb_tot > 0:
				sys.stdout.write(f"Results in --> {current_dir}{fileSep}loci{fileSep}{loc2cat}/{loc2cat}_allseq_select.fasta\n")

			if loci is not None:
				locusIndex = locusIndex + 1
			elif nb_tot > 0 and "yes" == promptUser("Do you want to generate dereplicated versions of the above mentioned fasta file? " + normalStyle + "Enter yes or no" + promptStyle, None, ["yes", "no"], 1, None, ""):
				print()
				logFile = open(f'{current_dir}{fileSep}outputs{fileSep}res3.log', 'w')
				outFasta = f"{current_dir}{fileSep}loci{fileSep}{loc2cat}/{loc2cat}" + '_allseq_select_derep.fasta'
				outTsv = f"{current_dir}{fileSep}loci{fileSep}{loc2cat}/{loc2cat}" + '_allseq_select_derep.tsv'
				derepResults = derepSeveralFastaFiles({loc2cat : loc2cat + '_allseq_select.fasta'}, loc2cat + '_allseq_select_derep.fasta', loc2cat + '_allseq_select_derep.tsv', logFile)
				print(successStyle + str(derepResults[0]) + " distinct sequences were dereplicated into fasta and tsv files: " + normalStyle + outFasta + ", " + outTsv + "\n")
				logFile.close()
		stat_3.close()

	if loci is not None and anySeqsFound is True and "yes" == promptUser("Do you want to generate dereplicated versions of the above mentioned fasta files? " + normalStyle + "Enter yes or no" + promptStyle, None, ["yes", "no"], 1, None, ""):
		print()
		logFile = open(f'{current_dir}{fileSep}outputs{fileSep}res3.log', 'w')
		for locus in loci:
			if locus not in skippedLoci:
				outFasta = f"{current_dir}{fileSep}loci{fileSep}{locus}/{locus}" + '_allseq_select_derep.fasta'
				outTsv = f"{current_dir}{fileSep}loci{fileSep}{locus}/{locus}" + '_allseq_select_derep.tsv'
				derepResults = derepSeveralFastaFiles({locus : f"{current_dir}{fileSep}loci{fileSep}{locus}/{locus}" + '_allseq_select.fasta'}, outFasta, outTsv, logFile)
				if derepResults[0] > 0:
					print(successStyle + str(derepResults[0]) + " distinct sequences were dereplicated into fasta and tsv files: " + normalStyle + outFasta + ", " + outTsv + normalStyle + "\n")
				else:
					print(warningStyle + "No sequences to process" + normalStyle + "\n")
		logFile.close()

	input("\nPress ENTER to continue ")


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
	print("\nPlease cite this software as follows:" +
		  citationStyle + "\nmbctools: A User-Friendly Metabarcoding and Cross-Platform Pipeline for Analyzing Amplicon Sequencing Data across a Large Diversity of Organisms"
						  "\nChristian Barnabé, Guilhem Sempéré, Vincent Manzanilla and Etienne Waleckx. https://github.com/GuilhemSempere/mbctools" + normalStyle + "\n\n")


def main_menu1():
	"""Displays submenu 1
	"""
	os.system("cls" if winOS else "clear")
	print(titleStyle + "\n---- MENU 1 - BASIC ANALYSIS - only option 1 is strictly mandatory ----" + normalStyle + "\n\n"
															   "1  -> NEW COMPLETE ANALYSIS (" + warningStyle + "mandatory" + normalStyle + ")\n"
															   "1a -> Re-analyze all loci, from the clustering step, modifying parameters\n"
															   "1b -> Re-analyze only one locus of paired-end amplicon (merged reads), modifying parameters\n"
															   "1c -> Re-analyze only one locus of single-end amplicon (R1 only), modifying parameters\n"
															   "1d -> Re-analyse a given sample, modifying parameters\n"
															   "1e -> " + warningStyle + "Optional" + normalStyle + " quality checking of fastq files (slow)\n" + normalStyle)
	global rmenu
	rmenu = promptUser("Please select an option among those listed above", None, ["1", "1a", "1b", "1c", "1d", "1e", "back", "home", "exit"], 1, main, "")

	if rmenu == "1e":
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
	print(titleStyle + "\n---- MENU 2 - REMOVAL OF PRIMERS, SELECTION OF MINIMUM SEQUENCE ABUNDANCES ACCORDING TO USER-DEFINED THRESHOLDS ----" + normalStyle + "\n\n"
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
																"\tsame as option 2c but only using the R1 sequences instead of merged ones.\n" + normalStyle)

	rmenu = promptUser("Please select an option among those listed above", None, ["2a", "2b", "2c", "2d", "back", "home", "exit"], 1, main, "")

	if rmenu == "2a":
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
	os.system("cls" if winOS else "clear")
	print(titleStyle + "\n---- MENU 3 - GENERATION OF A UNIQUE SEQUENCE FILE FOR EACH LOCUS (comprising all samples' data) ----" + normalStyle +
																"\n\n3a -> Process all loci at once\n"
																"3b -> Process loci one by one\n")

	rmenu = promptUser("Please select an option among those listed above", None, ["3a", "3b", "back", "home", "exit"], 1, main, "")
	concat_3(lociPEs + list(set(lociSEs) - set(lociPEs)) if rmenu == '3a' else None)
	main_menu3()


def main_menu4():
	"""Displays submenu 4
	"""
	os.system("cls" if winOS else "clear")
	global rmenu
	print(titleStyle + "\n---- MENU 4 - CONVERSION OF ANALYSIS RESULTS INTO metaXplor IMPORT FORMAT ----" + normalStyle + "\n\n"
																  "4a -> Generate sequence files\n"
																  "\tCompiles all sequences selected for all loci into a single fasta\n"
																  "\tOutputs a .tsv file indicating samples weights for each sequence\n\n"
																  "4b -> Generate assignment file\n"
																  "\tConverts blastn results (obtained from blasting above-mentioned fasta file) from 'Hit table (text)'"
																  "\n\t(format #7) into metaXplor format\n\n"
																  "4c -> Builds metaXplor-format sample metadata file from provided tabulated file\n\n"
																  "4d -> Compresses all metaXplor files into a final, ready to import, zip archive\n" + normalStyle)

	rmenu = promptUser("Please select an option among those listed above", None, ["4a", "4b", "4c", "4d", "back", "home", "exit"], 1, main, "")

	if rmenu == "4a":
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
		lociPE
	except NameError:
		in_lociPE()
	try:
		lociSE
	except NameError:
		in_lociSE()
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
	if len(lociPEs) > 0:
		merging()
	fastq2fas()
	derep_1()
	cluster_1x()
	chimera_remove()
	if len(lociPEs) > 0:
		runloc_merged()
	if len(lociSEs) > 0:
		runloc_r1()
	orient_1x()
	sys.stdout.write("\n\n")
	runs_1x()
	stats_1x()
	rerun(None)
	if len(sys.argv) == 0:
		rerun(None)


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
	print()
	cluster_1x()
	chimera_remove()
	runloc_merged()
	runloc_r1()
	orient_1x()
	runs_1x()
	stats_1x()
	rerun(None)


def menu1b():
	"""Runs option 1b
	"""
	prevent()
	prev_param(None)
	if len(lociPEs) == 0:
		print(errorStyle + "No paired-end reads in current selection!" + normalStyle)
		input("\nPress ENTER to continue ")
		main_menu1()
	else:
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
		rerun(None)


def menu1c():
	"""Runs option 1c
	"""
	prevent()
	prev_param(None)
	if len(lociSEs) == 0:
		print(errorStyle + "No single-end or unmerged R1 reads in current selection!" + normalStyle)
		input("\nPress ENTER to continue ")
		main_menu1()
	else:
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
		rerun(None)


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
	rerun(None)


def menu1e():
	"""Runs option 1e
	"""
	prevent()
	prev_param(None)
	quality()
	runs_1x()
	rerun(None)


def menu2a():
	"""Runs option 2a
	"""
	prevent()
	prev_param(None)
	if len(lociPEs) == 0:
		print(errorStyle + "No paired-end reads in current selection!" + normalStyle)
		input("\nPress ENTER to continue ")
		main_menu2()
	else:
		trim_2x()
		rerun(None)


def menu2b():
	"""Runs option 2b
	"""
	prevent()
	prev_param(None)
	if len(lociSEs) == 0:
		print(errorStyle + "No single-end or unmerged R1 reads in current selection!" + normalStyle)
		input("\nPress ENTER to continue ")
		main_menu2()
	else:
		trim_2x()
		rerun(None)


def menu2c():
	"""Runs option 2c
	"""
	prevent()
	prev_param(None)
	if len(lociPEs) == 0:
		print(errorStyle + "No paired-end reads in current selection!" + normalStyle)
		input("\nPress ENTER to continue ")
		main_menu2()
	else:
		if os.path.exists("outputs/Stats_option_2c.txt"):
			os.remove("outputs/Stats_option_2c.txt")
		trim_2x()
		rerun(None)


def menu2d():
	"""Runs option 2d
	"""
	prevent()
	prev_param(None)
	if len(lociSEs) == 0:
		print(errorStyle + "No single-end or unmerged R1 reads in current selection!" + normalStyle)
		input("\nPress ENTER to continue ")
		main_menu2()
	else:
		if os.path.exists("outputs/Stats_option_2d.txt"):
			os.remove("outputs/Stats_option_2d.txt")
		trim_2x()
		rerun(None)


def menu3():
	"""Runs option 3
	"""
	prevent()
	prev_param(None)
	concat_3()
	rerun(None)


def derepBasedOnIDs(locusToFastaDict, outFastaName, outTsvName):
	fastaContents = ""
	with open(outFastaName, "w") as fastaFile, open(outTsvName, "w") as seqCompositionFile:
		i = 0
		seqCompositionFile.write("qseqid")
		while i < len(samples):
			seqCompositionFile.write("\t" + samples[i])
			i += 1

		seqSampleAbundances = {}
		seqHashToMeaningfulNamesDict = {}
		skippedLoci = []
		totalDistinctSeqCount = 0
		for locus in locusToFastaDict:
			print("Dereplicating sample sequences by locus for " + locus + "...")
			try:
				file1 = open(locusToFastaDict[locus], 'r')
				lines = file1.readlines()

				j = 0
				seqCount = 0
				skipActualSequence = False # will be true when the sequence being read has already been written out (dereplication)
				while j < len(lines):
					line = lines[j].replace("\r\n", "").replace("\n", "")
					if line.startswith(">"):
						seqCount += 1
						splitIdLine = line[1:].split(" ")
						sampleAndCount = re.sub(r'_(merged|R1).*;', ';', splitIdLine[1].replace("sample=", "")).split(";size=")
						if splitIdLine[0] not in seqSampleAbundances:
							fastaContents += ">" + splitIdLine[0] + "\n"
							seqSampleAbundances[splitIdLine[0]] = ['0'] * len(samples)
							totalDistinctSeqCount += 1
							skipActualSequence = False
						else:
							skipActualSequence = True
						seqSampleAbundances[splitIdLine[0]][samples.index(sampleAndCount[0])] = sampleAndCount[1]
						if splitIdLine[0] not in seqHashToMeaningfulNamesDict:
							seqHashToMeaningfulNamesDict[splitIdLine[0]] = []
						seqHashToMeaningfulNamesDict[splitIdLine[0]].append(sampleAndCount[0])
					elif skipActualSequence is False:
						fastaContents += lines[j]
					j += 1
			except FileNotFoundError:
				print(warningStyle + "File not found: loci/" + locus + fileSep + locus + '_allseq_select.fasta: skipping locus ' + locus + normalStyle)
				skippedLoci.append(locus)
			i += 1
		
		i = 1
		padLevel = len(str(len(seqSampleAbundances)))
		for seqId in seqSampleAbundances:
			seqName = "seq" + str(i).zfill(padLevel) + "." + "_".join(seqHashToMeaningfulNamesDict[seqId][0:5]) + ("" if len(seqHashToMeaningfulNamesDict[seqId]) <= 5 else "...")
			fastaContents = fastaContents.replace(">" + seqId + "\n", ">" + seqName + "\n")
			seqCompositionFile.write("\n" + seqName + "\t" + "\t".join(seqSampleAbundances[seqId]))
			i += 1

		fastaFile.write(fastaContents)

	return [totalDistinctSeqCount, skippedLoci]


def derepSeveralFastaFiles(locusToFastaDict, outFastaName, outTsvName, logFile):
	with open(outFastaName + ".tmp1", 'wb') as wfd:
		for sp in locusToFastaDict:
			with open(locusToFastaDict[sp],'rb') as fd:
				shutil.copyfileobj(fd, wfd)

	subprocess.run(["vsearch", "--fastx_filter", outFastaName + ".tmp1", "--relabel_sha1", "--relabel_keep", "--fastaout", outFastaName + ".tmp2"], stderr=logFile)
	os.remove(outFastaName + ".tmp1")
	retVal = derepBasedOnIDs({", ".join(locusToFastaDict.keys()) : outFastaName + ".tmp2"}, outFastaName, outTsvName)
	os.remove(outFastaName + ".tmp2")

	return retVal


def menu4a():
	"""Runs option 4a
	"""
	prevent()
	prev_param(None)
	locusToFastaDict = {}
	for locus in list(set(lociPEs + lociSEs)):
		files2cat = glob.glob("loci/" + locus + fileSep + '*_allseq_select.fasta')
		if len(files2cat) > 0:
			locusToFastaDict[locus] = "loci/" + locus + fileSep + locus + '_allseq_select.fasta'

	print()
	logFile = open(f'{current_dir}{fileSep}outputs{fileSep}res4a.log', 'w')
	derepResults = derepSeveralFastaFiles(locusToFastaDict, metaXplorFasta, metaXplorSequenceComposition, logFile)
	logFile.close()

	if derepResults[0] == 0:
		print(errorStyle + "\nNo concatenated sequences found for any loci. Please run step 3 before retrying" + normalStyle)
		os.remove(metaXplorFasta)
		os.remove(metaXplorSequenceComposition)
		rerun(main_menu3)
	else:
		print(successStyle + "\n" + str(derepResults[0]) + " distinct sequences were compiled into .fasta and .tsv files")
		if len(derepResults[1]) > 0:
			print(
				warningStyle + "Warning: not all sequences could be included because concatenation step (#3) was not run on some loci: " + ", ".join(derepResults[1]) + successStyle)
		print("You may now run blastn on " + metaXplorFasta + ", download all results as 'Hit table (text)' (format #7), then come back and launch step 4b" + normalStyle)
		print(warningStyle + "\nNB: The appropriate format required by step 4b may be obtained by either:" + normalStyle
			+ "\n - if using NCBI online BLAST interface, selecting 'Hit table (text)' from the 'Download All' dropdown list"
			+ "\n - if executing command line BLAST, specifying the following argument: -outfmt 7\n"
			+ "(Its header must contain a line starting with: " + citationStyle + "'# Fields: query acc.ver, subject acc.ver, '" + normalStyle + ")")
		rerun(main_menu4)


def menu4b():
	"""Runs option 4b
	"""
	prevent()
	prev_param(None)

	print(warningStyle + "\nNB: The appropriate format required by step 4b may be obtained by either:" + normalStyle
		+ "\n - if using NCBI online BLAST interface, selecting 'Hit table (text)' from the 'Download All' dropdown list"
		+ "\n - if executing command line BLAST, specifying the following argument: -outfmt 7\n"
		+ "(Its header must contain a line starting with: " + citationStyle + "'# Fields: query acc.ver, subject acc.ver, '" + normalStyle + ")")

	blastTextHitTable = promptUser("Enter path to blastn hit-table (text format #7)", None, ["back", "home", "exit"], 3, main_menu4, "")
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


def menu4c():
	"""Runs option 4c
	"""
	prevent()
	prev_param(None)

	sampleMetadataFile = None
	print("\n\nYou must now provide a tabulated metadata file for your samples. A header field named 'Sample' is expected for the column featuring sample names")
	print("Any field with 'date' in its header name will be considered to be the sample collection date")
	print("Collection location may be specified:")
	print("\t- either as commma-separated decimal-format values in a single field named 'LatLon' (e.g. -17.7127, -67.9905)")
	print("\t- or in separate columns named 'Latitude' and 'Longitude', in decimal format (e.g. -17.7127) or DMS format (e.g. 16°42'45.6\"S)")
	print("Any additional columns will remain named as provided")

	sampleMetadataFile = promptUser("Enter path to tabulated sample metadata file", None, ["back", "home", "exit"], 3, main_menu4, "")
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


def menu4d(invokedByUser):
	"""Runs option 4d
	"""
	if not os.path.isfile(metaXplorFasta) or os.path.getsize(metaXplorFasta) == 0:
		if invokedByUser:
			print(errorStyle + "\nFile " + metaXplorFasta + " is missing or empty. Please run step 4a" + normalStyle)
		rerun(main_menu4)
	if not os.path.isfile(metaXplorSequenceComposition) or os.path.getsize(metaXplorSequenceComposition) == 0:
		if invokedByUser:
			print(errorStyle + "\nFile " + metaXplorSequenceComposition + " is missing or empty. Please run step 4a" + normalStyle)
		rerun(main_menu4)
	if not os.path.isfile(metaXplorAssignments) or os.path.getsize(metaXplorAssignments) == 0:
		if invokedByUser:
			print(errorStyle + "\nFile " + metaXplorAssignments + " is missing or empty. Please run step 4b" + normalStyle)
		rerun(main_menu4)
	if not os.path.isfile(metaXplorSamples) or os.path.getsize(metaXplorSamples) == 0:
		if invokedByUser:
			print(errorStyle + "\nFile " + metaXplorSamples + " is missing or empty. Please run step 4c" + normalStyle)
		rerun(main_menu4)

	global zipNow
	print("\n\nAll metaXplor files seem to be ready for zipping.")
	zipNow = promptUser("Zip them now to create the final import file? " + normalStyle + "Enter yes or no" + promptStyle, None, ["yes", "no"], 1, None, "")
	if zipNow == "yes":
		b = io.BytesIO()
		zf = zipfile.ZipFile(b, mode='w')
		zf.write(metaXplorSamples, os.path.basename(metaXplorSamples))
		zf.write(metaXplorAssignments, os.path.basename(metaXplorAssignments))
		zf.write(metaXplorFasta, os.path.basename(metaXplorFasta))
		zf.write(metaXplorSequenceComposition, os.path.basename(metaXplorSequenceComposition))
		zf.close()
		zipFileName = 'mbctools_metaXplor_export_' + date.strftime('%Y%m%d') + '.zip'
		open(zipFileName, 'wb').write(b.getbuffer())
		print(successStyle + "\n\nmetaXplor import archive was successfully created as " + current_dir + fileSep + zipFileName + normalStyle)
		rerun(None)
	else:
		main_menu4()


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
		command = 'powershell.exe Get-ExecutionPolicy'
		process = subprocess.Popen(command, stdout=subprocess.PIPE)
		result = process.communicate()[0].strip()
		execution_policy = result.decode()
		if execution_policy != 'Unrestricted':
			print(errorStyle + "\nPlease run the following command in Windows prompt to be able to use mbctools:\n" + normalStyle + "\nSet-ExecutionPolicy Unrestricted -Scope CurrentUser -Force\n")
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
			customExit(0)
		else:
			print(errorStyle + "\nUnexisting configuration file: " + sys.argv[1] + "\n" + normalStyle)
			customExit(1)

	os.system("cls" if winOS else "clear")

	sys.stdout.write(titleStyle + "\n------------------------------ mbctools - MAIN MENU ------------------------------" + normalStyle + "\n")
	printHowToCite()
	sys.stdout.write(titleStyle + "NAVIGATION CONVENTIONS:\n" + warningStyle + "Validating without typing anything applies the default value, if any\n" + normalStyle +
					 "Entering '" + promptStyle + "back" + normalStyle + "' returns to the program upper level, if any\n"
					 "Entering '" + promptStyle + "home" + normalStyle + "' returns to this main menu\n"
					 "Entering '" + promptStyle + "exit" + normalStyle + "' leaves the program\n\n")
	print(titleStyle + "\nWe recommend executing procedures in the provided order:\n" + normalStyle)
	print(	"1 -> BASIC ANALYZES\n\n"
			"2 -> REMOVAL OF PRIMERS AND SELECTION OF MINIMUM SEQUENCE ABUNDANCES ACCORDING TO USER-DEFINED THRESHOLDS\n\n"
			"3 -> GENERATION OF A UNIQUE SEQUENCE FILE FOR EACH LOCUS (comprising all samples' data)\n\n"
			"4 -> CONVERSION OF ANALYSIS RESULTS INTO metaXplor IMPORT FORMAT\n")

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
	input("\n(Press ENTER to exit)")
	exit(code)


if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		print("\n")
