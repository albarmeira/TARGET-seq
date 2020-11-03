#High throuhgput single cell genotyping of 3'-TARGETseq dataset
#Author:Alba Rodriguez Meira
#TARGET-seq : Unravelling intratumoral heterogeneity through high-sensitivity single-cell mutational analysis and parallel RNA-sequencing


library(ggplot2)
library(stringr)

metadata<-read.table('metadata_HT_TARGETseq.txt',header = T,sep = '\t')
rownames(metadata)<-metadata$cell.id

############# ANALYSIS OF gJAK2 #############

gJAK2<-read.table('counts/gDNA_JAK2_V617F.txt',header = T,sep = '\t')
#Remove _S pattern from the sequencer
gJAK2$id<-sapply(strsplit(as.character(gJAK2$id),"_",fixed = TRUE),"[[",1)
gJAK2$id<-gsub("-","_",gJAK2$id)
rownames(gJAK2)<-gJAK2$id
#Calculate coverage
gJAK2$coverage<-(gJAK2$C+gJAK2$A+gJAK2$G+gJAK2$T+gJAK2$del+gJAK2$ins+str_count(gJAK2$ambiguous, ">"))
#Select cells sequenced for the analysis, merge with metadata
metadata<-metadata[gJAK2$id,]
gJAK2<-merge(gJAK2,metadata,by='row.names',all=TRUE)

# Calculate VAF from mutant and WT alleles


gJAK2$MUTANT<-(gJAK2$T)/(gJAK2$coverage)
gJAK2$WT<-(gJAK2$G)/(gJAK2$coverage)

#Plot coverage per donor; set cut-off at 3% and 97% VAF. Regions [0.03-0.1] and [0.9-0.97]
#are usually considered indetermined.

ggplot(gJAK2, aes(x=gJAK2$MUTANT, y=gJAK2$coverage,colour=donor))+geom_point()+
  geom_vline(xintercept=0.03, colour="red")+
  geom_vline(xintercept=0.1, colour="black")+
  geom_vline(xintercept=0.9, colour="black")+
  geom_vline(xintercept=0.97,colour="red")+
  theme_bw()+ggsave('gJAK2_coverage_vs_VAF.pdf',height=7,width=6)

#Fix the limit of detection based on blank samples

gJAK2$Detection<-gJAK2$coverage>175

############# ANALYSIS OF mJAK2

mJAK2<-read.table('counts/mRNA_JAK2_V617F.txt',header = T,sep = '\t')
mJAK2$id<-sapply(strsplit(as.character(mJAK2$id),"_",fixed = TRUE),"[[",1)
mJAK2$id<-gsub("-","_",mJAK2$id)
rownames(mJAK2)<-mJAK2$id
mJAK2$coverage<-(mJAK2$C+mJAK2$A+mJAK2$G+mJAK2$T+mJAK2$del+mJAK2$ins+str_count(mJAK2$ambiguous, ">"))
metadata<-metadata[mJAK2$id,]
mJAK2<-merge(mJAK2,metadata,by='row.names',all=TRUE)

#Calculate mutant and WT vaf

mJAK2$MUTANT<-(mJAK2$T)/(mJAK2$coverage)
mJAK2$WT<-(mJAK2$G)/(mJAK2$coverage)

ggplot(mJAK2, aes(x=mJAK2$MUTANT, y=mJAK2$coverage,colour=donor))+geom_point()+
  geom_vline(xintercept=0.03, colour="red")+
  geom_vline(xintercept=0.1, colour="black")+
  geom_vline(xintercept=0.9, colour="black")+
  geom_vline(xintercept=0.97,colour="red")+
  theme_bw()+ggsave('mJAK2_coverage_vs_VAF.pdf',height=7,width=6)

#Fix the limit of detection based on blank samples

mJAK2$Detection<-mJAK2$coverage>150

##############################################
#Merge mRNA and gDNA information for each gene
##############################################

#Assign genotypes for gDNA amplicon
gJAK2$gJAK2genotype[gJAK2$MUTANT<0.03]<-'WT' 
gJAK2$gJAK2genotype[gJAK2$MUTANT>0.03 & gJAK2$MUTANT<=0.1]<-'NA' 
gJAK2$gJAK2genotype[gJAK2$MUTANT>0.1 & gJAK2$MUTANT<=0.90]<-'HET'
gJAK2$gJAK2genotype[gJAK2$MUTANT>0.90 & gJAK2$MUTANT<=0.97]<-'NA' 
gJAK2$gJAK2genotype[gJAK2$MUTANT>0.97]<-'HOM'
gJAK2$gJAK2genotype[gJAK2$Detection == "FALSE"] <- NA

#Assign genotypes for mRNA amplicon
mJAK2$mJAK2genotype[mJAK2$MUTANT<0.03]<-'WT' 
mJAK2$mJAK2genotype[mJAK2$MUTANT>0.03 & mJAK2$MUTANT<=0.1]<-'NA' 
mJAK2$mJAK2genotype[mJAK2$MUTANT>0.1 & mJAK2$MUTANT<=0.90]<-'HET'
mJAK2$mJAK2genotype[mJAK2$MUTANT>0.90 & mJAK2$MUTANT<=0.97]<-'NA' 
mJAK2$mJAK2genotype[mJAK2$MUTANT>0.97]<-'HOM'
mJAK2$mJAK2genotype[mJAK2$Detection == "FALSE"] <- NA

#Merge genotype from gDNA and mRNA amplicons
gJAK2$JAK2all<-paste(gJAK2$gJAK2genotype,mJAK2$mJAK2genotype,sep="_")

#Now assing a consensus genotype based on gDNA and mRNA amplicons

gJAK2$JAK2all[grepl("HET",gJAK2$JAK2all)]<-'JAK2_HET'
gJAK2$JAK2all[grepl("HOM_WT|WT_HOM",gJAK2$JAK2all)]<-'JAK2_HET'
gJAK2$JAK2all[grepl("HOM_HOM|HOM_NA",gJAK2$JAK2all)]<-'JAK2_HOM'
gJAK2$JAK2all[grepl("WT_WT|WT_NA",gJAK2$JAK2all)]<-'JAK2_WT'
gJAK2$JAK2all[grepl("NA_HOM|NA_WT",gJAK2$JAK2all)]<-NA #If gDNA is not detected but mRNA is WT or HOM, 
#we can't trust the genotype so it'll be unassigned (due to the high allelic dropout rate of mRNA amplicons)
gJAK2$JAK2all[grepl("NA_NA",gJAK2$JAK2all)]<- NA

############# ANALYSIS OF U2AF1 ##############

###### ANALYSIS OF gU2AF1 #######

gU2AF1<-read.table('counts/gDNA_U2AF1_pQ157.txt',header = T,sep = '\t')
gU2AF1$id<-sapply(strsplit(as.character(gU2AF1$id),"_",fixed = TRUE),"[[",1)
gU2AF1$id<-gsub("-","_",gU2AF1$id)
rownames(gU2AF1)<-gU2AF1$id
gU2AF1$coverage<-(gU2AF1$C+gU2AF1$A+gU2AF1$G+gU2AF1$T+gU2AF1$del+gU2AF1$ins+str_count(gU2AF1$ambiguous, ">"))
metadata<-metadata[gU2AF1$id,]
gU2AF1<-merge(gU2AF1,metadata,by='row.names',all=TRUE)

# 

gU2AF1$MUTANT<-(gU2AF1$C)/(gU2AF1$coverage)
gU2AF1$WT<-(gU2AF1$T)/(gU2AF1$coverage)

ggplot(gU2AF1, aes(x=gU2AF1$MUTANT, y=gU2AF1$coverage,colour=donor))+geom_point()+
  geom_vline(xintercept=0.035, colour="red")+
  geom_vline(xintercept=0.1, colour="black")+
  geom_vline(xintercept=0.9, colour="black")+
  geom_vline(xintercept=0.965,colour="red")+
  theme_bw()+ggsave('gU2AF1_coverage_vs_VAF.pdf',height=7,width=6)

#Fix the limit of detection based on blank samples

gU2AF1$Detection<-gU2AF1$coverage>55

############# ANALYSIS OF mU2AF1

mU2AF1<-read.table('counts/mRNA_U2AF1_pQ157.txt',header = T,sep = '\t')
#Remove _S pattern from the sequencer
mU2AF1$id<-sapply(strsplit(as.character(mU2AF1$id),"_",fixed = TRUE),"[[",1)
mU2AF1$id<-gsub("-","_",mU2AF1$id)
rownames(mU2AF1)<-mU2AF1$id
#Calculate coverage
mU2AF1$coverage<-(mU2AF1$C+mU2AF1$A+mU2AF1$G+mU2AF1$T+mU2AF1$del+mU2AF1$ins+str_count(mU2AF1$ambiguous, ">"))
#Select cells sequenced for the analysis, merge with metadata
metadata<-metadata[mU2AF1$id,]
mU2AF1<-merge(mU2AF1,metadata,by='row.names',all=TRUE)

# Calculate VAF from mutant and WT alleles
mU2AF1$MUTANT<-(mU2AF1$C)/(mU2AF1$coverage)
mU2AF1$WT<-(mU2AF1$T)/(mU2AF1$coverage)

#Plot coverage versus donor
ggplot(mU2AF1, aes(x=mU2AF1$MUTANT, y=mU2AF1$coverage,colour=donor))+geom_point()+
  geom_hline(yintercept=100)+ #limit of detection
  geom_vline(xintercept=0.03, colour="red")+
  geom_vline(xintercept=0.1, colour="black")+
  geom_vline(xintercept=0.9, colour="black")+
  geom_vline(xintercept=0.97,colour="red")+
  theme_bw()+ggsave('output/mU2AF1_coverage_vs_VAF.pdf',height=7,width=6)

#Fix the limit of detection based on blanks samples

mU2AF1$Detection<-mU2AF1$coverage>100


##############################################
#Merge mRNA and gDNA information for each gene
##############################################

#Assign genotypes for gDNA amplicon
gU2AF1$gU2AF1genotype[gU2AF1$MUTANT<0.03]<-'WT' 
gU2AF1$gU2AF1genotype[gU2AF1$MUTANT>0.03 & gU2AF1$MUTANT<=0.1]<-'NA' 
gU2AF1$gU2AF1genotype[gU2AF1$MUTANT>0.1 & gU2AF1$MUTANT<=0.90]<-'HET'
gU2AF1$gU2AF1genotype[gU2AF1$MUTANT>0.90 & gU2AF1$MUTANT<=0.97]<-'NA' 
gU2AF1$gU2AF1genotype[gU2AF1$MUTANT>0.97]<-'HOM'
gU2AF1$gU2AF1genotype[gU2AF1$Detection == "FALSE"] <- NA

#Assign genotypes for mRNA amplicon
mU2AF1$mU2AF1genotype[mU2AF1$MUTANT<0.03]<-'WT' 
mU2AF1$mU2AF1genotype[mU2AF1$MUTANT>0.03 & mU2AF1$MUTANT<=0.1]<-'NA' 
mU2AF1$mU2AF1genotype[mU2AF1$MUTANT>0.1 & mU2AF1$MUTANT<=0.90]<-'HET'
mU2AF1$mU2AF1genotype[mU2AF1$MUTANT>0.90 & mU2AF1$MUTANT<=0.97]<-'NA' 
mU2AF1$mU2AF1genotype[mU2AF1$MUTANT>0.97]<-'HOM'
mU2AF1$mU2AF1genotype[mU2AF1$Detection == "FALSE"] <- NA

#Merge genotype from gDNA and mRNA amplicons
gU2AF1$U2AF1all<-paste(gU2AF1$gU2AF1genotype,mU2AF1$mU2AF1genotype,sep="_")

#Now assing a consensus genotype based on gDNA and mRNA amplicons

gU2AF1$U2AF1all[grepl("HET",gU2AF1$U2AF1all)]<-'U2AF1_HET'
gU2AF1$U2AF1all[grepl("HOM_WT|WT_HOM",gU2AF1$U2AF1all)]<-'U2AF1_HET'
gU2AF1$U2AF1all[grepl("HOM_HOM|HOM_NA",gU2AF1$U2AF1all)]<-'U2AF1_HOM'
gU2AF1$U2AF1all[grepl("WT_WT|WT_NA",gU2AF1$U2AF1all)]<-'U2AF1_WT'
gU2AF1$U2AF1all[grepl("NA_HOM|NA_WT",gU2AF1$U2AF1all)]<-NA #If gDNA is not detected but mRNA is WT or HOM, 
#we can't trust the genotype so it'll be unassigned (due to the high allelic dropout rate of mRNA amplicons)
gU2AF1$U2AF1all[grepl("NA_NA",gU2AF1$U2AF1all)]<- NA


############################################################
## MERGING PROCESS
############################################################
metadata$genotypesall <- paste(gJAK2$JAK2all,
                   gU2AF1$U2AF1all,
                   sep="_")

#Remove cells with missing genotypes
metadata$genotypesall<-gsub("NA_","",metadata$genotypesall)
metadata$genotypesall<-gsub("_NA","",metadata$genotypesall)

metadata$qcgenotype<-"FALSE"
metadata$qcgenotype[grepl("U2AF1",metadata$genotypesall) & grepl("JAK2",metadata$genotypesall)]<-"TRUE"

passqc<-subset(metadata,qcgenotype=="TRUE")
table(passqc$genotypesall)

############################################################
IF0137<-metadata[metadata$donor=="IF0137" & metadata$qcgenotype=="TRUE",] #778 cells pass QC
table(IF0137$genotypesall)
IF0137.summary<-as.data.frame(table(IF0137$genotypesall))
write.table(IF0137.summary,file="output/IF0137.summary.txt",sep="\t",row.names = F)