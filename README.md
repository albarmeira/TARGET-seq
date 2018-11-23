# TARGET-seq <img align="left" width="100" height="100" src="TARGET-seq/target.pdf">

"Unravelling intratumoral heterogeneity through high-sensitivity single-cell mutational analysis and parallel RNA-sequencing"

Rodriguez-Meira, A., Buck, G., Clark, S.-A., Povinelli, B.J., Alcolea, V., Louka, E., McGowan, S., Hamblin, A., Sousos, N., Barkas, N., et al. (2018). Unraveling intratumoral heterogeneity through high-sensitivity single-cell mutational analysis and parallel RNA-sequencing. bioRxiv.

https://doi.org/10.1101/474734

# TARGET-seq single cell genotyping pipeline (SCpipeline)

This pipeline provides an integrated framework for the analysis of single-cell targeted sequencing data obtained using TARGET-seq method. The pipeline processes raw fastq files (which can be gzip compressed) or aligned and sorted bam files, and outputs counts tables for every mutation assessed. 

Authors: Alba Rodriguez-Meira, Simon McGowan.

## Availability

SCpipeline is freely available under a GPL3 license.

## How to run 

To run the pipeline, open the terminal and type:
```
perl SCgenotype.pl -c SCpipe.conf
```

The pipeline requires python version 2.7, STAR version 2.4.2a and samtools version 1.1.

## Input files and parameters

Directories and input files are specified through SCpipe.conf file:

```FASTQ_DIR </pathto/fastqdirectory/>```
Replace </pathto/fastqdirectory/> with the path to directory in which all fastq files to be analyzed are stored

```FASTQ_FWD_SUFFIX <_LXX_RX_.fastq.gz>``` 
Replace <_LXX_R1.fastq.gz> with the suffix of read 1 files; in the case of MiSeq platform this is _L001_R1_001.fastq.gz

```FASTQ_REV_SUFFIX  <_LXX_RX_.fastq.gz>``` 
Replace <_LXX_R2.fastq.gz> with the suffix of read 2 files; in the case of MiSeq platform this is _L001_R2_001.fastq.gz

```gPRIMERS_BED </pathto/gPRIMERS.bed>``` 
Replace </pathto/gPRIMERS.bed> with the path to the bed file containing chromosomal coordinates for each primer used to target genomic DNA amplicons. An example file used for analysis of 3'-TARGETseq dataset can be downloaded from this page (See "gPRIMERS.bed").

```mPRIMERS_BED </pathto/mPRIMERS.bed>``` 
Replace </pathto/mPRIMERS.bed> with the path to bed file containing chromosomal coordinates for each primer used to target genomic DNA amplicons. An example file used for analysis of 3'-TARGETseq dataset can be downloaded from this page (See "mPRIMERS.bed").

```VARIANTS </pathto/variants.tsv>``` 
Replace </pathto/variants.tsv> with a tsv file containing coordinates of each point mutation/indel to be assessed. Start and end coordinates should be the same for all mutations, including indels. Position of insertions and deletions should be - 1bp the annotated position (identified through a bulk sequencing mutational calling pipeline). For example, ASXL1 c.2728_2729del, annotated in chr20:31023242, should be introduced in the variants.tsv file as chr20:31023241 (annotated position minus 1 base). An example file used for analysis of 3'-TARGETseq dataset can be downloaded from this page (See "variants.tsv").

## Output directory 
```ANALYSIS_DIR </pathto/outdir/>``` 
Replace </pathto/outdir/> with the path to the directory where output count files should be stored. Intermediate directories used for analysis will also be created in here.

## Executable modules

```STAR_EXECUTABLE </pathto/package/rna-star/2.4.2a/bin/STAR>``` 
Replace </pathto/package/rna-star/2.4.2a/bin/STAR> with path to STAR executable module.

```STAR_INDEX_DIR </pathto/Homo_sapiens/UCSC/hg19/Sequence/STAR>``` Replace </pathto/Homo_sapiens/UCSC/hg19/Sequence/STAR> with path to STAR index directory. In the analysis of TARGET-seq paper, we used hg19.

```STAR_NUM_THREADS <number>``` 
Replace <number> with the number of threads to use when running STAR

```SAMTOOLS_EXECUTABLE </pathto/package/samtools/1.1/bin/samtools>``` 
Replace </pathto/package/samtools/1.1/bin/samtools> with the path to samtools executable module. 

```SAMTOOLS_NUM_THREADS <number>``` 
Replace <number> with the number of threads to use when running samtools

```SEPARATE_gDNA_SCRIPT </pathto/scripts/genome_reads_gDNA.pl>``` 
Replace </pathto/scripts/genome_reads_gDNA.pl> with the path to genome_reads_gDNA.pl script.

```SEPARATE_mRNA_SCRIPT </scripts/genome_reads_mRNA.pl>```
Replace </pathto/scripts/genome_reads_mRNA.pl> with the path to genome_reads_mRNA.pl script.

```SUMMARISE_SCRIPT </pathto/scripts/getcounts_mpileup_TARGET.py>```
Replace </pathto/scripts/getcounts_mpileup_TARGET.py> with the path to getcounts_mpileup_TARGET.py script.
