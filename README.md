# TARGET-seq

"Unravelling intratumoral heterogeneity through high-sensitivity single-cell mutational analysis and parallel RNA-sequencing"

Rodriguez-Meira et al. (citation pending)

# TARGET-seq single cell genotyping pipeline

This pipeline provides an integrated framework for the analysis of single-cell targeted sequencing data obtained using TARGET-seq method. The pipeline processes raw fastq files or aligned and sorted bam files, and outputs counts tables for every mutation assessed. 

# How to run 

perl SCgenotype.pl -c SCpipe.conf

The pipeline requires python version 2.7

# Input files
Directories and input files are specified through SCpipe.conf file:

/pathto/fastqdirectory/
