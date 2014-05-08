The scripts contained in this directory are meant to provide a low and high resolution SNP analysis of closely related isolates.  Individuals should customize scripts to cater to the throughput of their lab.  The two main scripts are: processZips.sh and vcftofasta.sh.  processZips outputs BAM and VCFs from Illumina paired-end data.  processZips.sh is called on a working directory containing paired files with “R1” and “R2” designations.  VCFs output from processZips.sh are collected into a single working directory and vcftofasta.sh is called.  Output are alignment files to be viewed in a tree viewer of chose, and SNP tables which provide a high resolution view of SNP data.

The first script ran on 20 isolates from a typical MiSeq run using 2 x 250 chemistry outputting 3.3 mb genomes with approximately 150X average depth of coverage on an Apple Mac Pro with 2 x 6-core Xeon processors and 48 GB of memory will complete in 1.5 hours.  Following the first script a collection of a 1000 VCFs will analyze through vcftofasta.sh in 2 hours.  Although more processing and memory is favorable, scripts have been tested successfully using a Macbook laptop with an Intel Core 2 Duo processor and 4 GB of memory.

Currently the main use of this workflow is for genotyping Mycobacterium tuberculosis complex and Brucella species. 

Paradigm
1) Once a SNP occurs and establishes in a population it does not revert back
2) Observations of homoplasy are rare
3) Group, subgroup and clade clusters only show informative SNPs for the isolates within that cluster
4) SNPs observed in a single isolate are less informative than SNPs seen in multiple isolates and therefore established in a population

