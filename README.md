# mbctools

This crossplatform (Linux, OSX, Windows) Python program is a toolkit to make the use of VSEARCH easier and interactive, helping to analyze your
Metabarcoding NGS data in best conditions. It features the following MAIN MENU:
<pre>
1 -> INITIAL ANALYSIS (mandatory): read merging, sample-level dereplication, sequence clustering,
	chimera detection, affiliation of sequences to loci, and sequence re-orientation

2 -> PRIMER REMOVAL AND SELECTION OF MINIMUM SEQUENCE ABUNDANCE LEVELS ACCORDING TO USER-DEFINED THRESHOLDS

3 -> GENERATION OF A UNIQUE SEQUENCE FILE FOR EACH LOCUS (comprising all samples' data)

4 -> EXPORTING ANALYSIS RESULTS INTO metaXplor FORMAT
</pre>

#### mbctools reference:
Christian Barnabé†, Guilhem Sempéré†, Vincent Manzanilla, Joel Moo Millan, Antoine Amblard-Rambert and Etienne Waleckx (†Co–first authors).
<strong><em>mbctools: A User-Friendly Metabarcoding and Cross-Platform Pipeline for Analyzing Multiple Amplicon Sequencing Data across a Large Diversity of Organisms.</em></strong>
bioRxiv doi: https://doi.org/10.1101/2024.02.08.579441

  
#### VSEARCH reference:
Rognes T, Flouri T, Nichols B, Quince C, Mahe F (2016).
<strong><em>VSEARCH: a versatile open source tool for metagenomics.</em></strong>
PeerJ 4:e2584 doi: 10.7717/peerj.2584 https://doi.org/10.7717/peerj.2584

#### metaXplor reference:
Guilhem Sempéré, Adrien Pétel, Magsen Abbé, Pierre Lefeuvre, Philippe Roumagnac, Frédéric Mahé, Gaël Baurens, Denis Filloux.
<strong><em>metaXplor: an interactive viral and microbial metagenomic data manager.</em></strong>
GigaScience, Volume 10, Issue 2, February 2021, giab001, https://doi.org/10.1093/gigascience/giab001

---

## Requirements

- Python v3.7 or higher
- VSEARCH v2.19.0 or higher must be installed and accessible in the PATH (binaries for all platforms available at https://github.com/torognes/vsearch, along with installation procedure descriptions)
- Windows users will need to enable Powershell script execution using the Set-ExecutionPolicy command 
