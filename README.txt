The scripts contained in this directory are ment to provide a high resolution SNP analysis of closely related isolates.

This script is the second script in a two script workflow.  Script 2 genotypes Mycobacterium tuberculosis complex and Brucella species from SNP data contained in VCFs.  It operates on VCFs generated with the same reference output from script 1.  VCFs are collected into a single working directory.  Comparisons are output as SNP tables and alignment files to view as trees in your program of choice.

Script 2 will run and output tables and alignments from just the data generated from script 1, however the data will be more informative if additional files are provide.  Those files are:
1) A file that contains positions to cluster individual isolates/VCFs into groups, subgroups and clades.
2) Files that contain positions to remove from the analysis.

Paradigm
1) Once a SNP occurs and establishes in a population it does not revert back
2) Observations of homoplasy are rare
3) Group, subgroup and clade clusters only show parsimony informative SNPs for the isolates within that cluster
4) SNPs observed in a single isolate are less informative than SNPs seen in multiple isolates and therefore established in a population

