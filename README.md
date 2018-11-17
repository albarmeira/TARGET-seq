# TARGET-seq

"Unravelling intratumoral heterogeneity through high-sensitivity single-cell mutational analysis and parallel RNA-sequencing"

Rodriguez-Meira et al. (citation pending)

# TARGET-seq single cell genotyping pipeline

This pipeline provides an integrated framework for the analysis of single-cell targeted sequencing data obtained using TARGET-seq method. The pipeline processes raw fastq files (which can be gzip compressed) or aligned and sorted bam files, and outputs counts tables for every mutation assessed. 

# How to run 

To run the pipeline, open the terminal and type:
```
perl SCgenotype.pl -c SCpipe.conf
```

The pipeline requires python version 2.7, STAR version 2.4.2a and samtools version 1.1.

# Input files and parameters

Directories and input files are specified through SCpipe.conf file:

FASTQ_DIR /pathto/fastqdirectory/ - path to directory in which all fastq files to be analyzed are stored

FASTQ_FWD_SUFFIX/ *_LXX_RX_.fastq.gz - suffix of read 1 files, in the case of MiSeq platform this is _L001_R1_001.fastq.gz

FASTQ_REV_SUFFIX  *_LXX_RX_.fastq.gz -  suffix of read 2 files, in the case of MiSeq platform this is _L001_R2_001.fastq.gz

gPRIMERS_BED /pathto/gPRIMERS.bed - path to bed file containing chromosomal coordinates for each primer used to target genomic DNA amplicons. An example file used for analysis of 3'-TARGETseq dataset can be downloaded from this page ("gPRIMERS.bed").

mPRIMERS_BED /pathto/mPRIMERS.bed - path to bed file containing chromosomal coordinates for each primer used to target genomic DNA amplicons. An example file used for analysis of 3'-TARGETseq dataset can be downloaded from this page ("mPRIMERS.bed").

VARIANTS /pathto/variants.tsv - tsv file containing coordinates of each point mutation/indel to be assessed. Start and end coordinates should be the same for all mutations, including indels. Position of insertions and deletions should be - 1bp the annotated position (identified through a bulk sequencing mutational calling pipeline). For example, ASXL1 c.2728_2729del, annotated in chr20:31023242, should be introduced in the variants.tsv file as chr20:31023241 (annotated position minus 1 base).

# Output directory 

# Executable modules

STAR_EXECUTABLE /pathto/package/rna-star/2.4.2a/bin/STAR - path to STAR executable module

STAR_INDEX_DIR /pathto/Homo_sapiens/UCSC/hg19/Sequence/STAR - path to STAR index directory

STAR_NUM_THREADS 4 - number of threads to use when running STAR

SAMTOOLS_EXECUTABLE /pathto/package/samtools/1.1/bin/samtools - path to samtools executable module

SAMTOOLS_NUM_THREADS 4 - number of threads to use when running samtools

SEPARATE_gDNA_SCRIPT ./scripts/genome_reads_gDNA.pl - number of threads to use when running samtools

SEPARATE_mRNA_SCRIPT ./scripts/genome_reads_mRNA.pl

SUMMARISE_SCRIPT ./scripts/getcounts_mpileup_TARGET.py
