# mbctools

This program is a toolkit to make the use of VSEARCH easier and interactive, to analyze your
Metabarcoding NGS data in best conditions. It proposes the following MAIN MENU:
<pre>
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
</pre>
  
VSEARCH reference:
Rognes T, Flouri T, Nichols B, Quince C, Mahe F (2016)
VSEARCH: a versatile open source tool for metagenomics
PeerJ 4:e2584 doi: 10.7717/peerj.2584 https://doi.org/10.7717/peerj.2584
