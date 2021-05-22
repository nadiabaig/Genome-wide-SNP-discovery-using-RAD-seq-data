## Single nucleotide polymorphism (SNP) detection using RAD-seq data

This bash script allows to do  a basic hand's own on your NGS RAD-seq data.


# Background
Restriction site-associated DNA sequencing (RADseq) is a fractional genome strategy designed to allow analysis of some subsets of genomes of the individuals. Working principle of RAD-seq is based upon the fragmentation of target genome using restriction enzyme digestion followed by sequencing of fragments using Next generation sequencing (NGS) platforms such as Illumina. Depending upon the project/goal sequencing can be done in single or paired end mode. The sequenced data can be utilized to identify genetic variants known as single nucleotide polymorphisms (SNP's) in a genome using various publically available NGS tools. 

# Objectives

    Designing a bash script based variant analysis work flow using raw RAD_seq data
    RAD-seq data manipulation and pre-treatment
    SNP calling from RAD-sequencing data
    Hard filtering of SNPs
    
 ## Variant analysis workflow
 
 # Demultiplexing reads

The sequences are raw sequences from the sequencing machine, without any pretreatments. They need to be demultiplexed to use for the further analysis.

# Quality control (Adapter trimming or read trimming)
It is better to get a quick impression of your quality controlled data before proceeding to the next steps of variant analysis. Re-cheking the data quality   reduces the likelihood of false positive SNP calls. If you feel your data is not cleaned properly, try to tweak the parameters of the quality control step.
 
# Genome indexing
Indexing a reference genome is pre-requisite of read mapping.


# Mapping, checking bam file statistics, variant calling and filtering
To make sense of the reads, they must be aligned to the reference genome of interest. The repeat masked reference genome can be downloaded from the following link: http://ensembl.gramene.org/Zea_mays/Info/Index. Note: this is for maize genome. Choose the reference for your reads accordingly. If your genome is not available in concatenated form (i.e sequence of individual chromosomes is available) then concatenate it using linux command 'cat'.

It is important to check the stats (alignment %) of your reads before processing your files further. Make sure to check the statistics, then do variant calling. Sometime the SNPs in your VCF file seems to be false positives. It's important to clean your data before publishing and applying any statistics or doing any analysis.





