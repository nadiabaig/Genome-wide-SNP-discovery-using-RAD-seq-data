## Single nucleotide polymorphism (SNP) detection using RAD-seq data

This bash script allows to do  a basic hand's on on your NGS RAD-seq data.


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


## Installing NGS data processing tools

# Samtools
    Download the current version of samtools from: http://www.htslib.org/download/
    cd samtools-1.x    # and similarly for bcftools and htslib
    ./configure --prefix=/where/to/install
    make
    make install
    
    Modify your .bashrc file so that when you type "samtools" it calls the program:
    export PATH=$PATH:/directory/samtools-version
  
 # BWA
      Download the latest version of BWA and unpack the file using
      tar xvzf bwa-0.7.12.tar.bz2
      followed by
      cd bwa-0.7.12
      make
      
 # Picard
     Download the latest version from their website using wget ----
      tar xjf picard-tools-2.4.1.zip 
     Picard tools are distributed as a pre-compiled Java executable (jar file) so there is no need to compile them. You can use them by directly 
     specifying the   path of .jar file of picard or can set an enviroment variable by using following one liner
     PICARD = "~/src/jars/picard.jar"
     
 # GATK
       tar xjf GenomeAnalysisTK-3.3-0.tar.bz2 
       No need of compilation as GATK is distributed as a Java executable .jar file
       

# STACKS
    Download from here: https://catchenlab.life.illinois.edu/stacks/
    
    tar xfvz stacks-2.xx.tar.gz
    cd stacks-2.xx
    ./configure
    make -j 8
    sudo make install

## References

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., ... & Durbin, R. (2009). The sequence alignment/map format and SAMtools. Bioinformatics, 25(16), 2078-2079.

Catchen J, Hohenlohe PA, Bassham S, Amores A, Cresko WA: Stacks: an analysis tool set for population genomics. Mol Ecol. 2013, 22: 3124-3140. 10.1111/mec.12354.

Jo, H., & Koh, G. (2015). Faster single-end alignment generation utilizing multi-thread for BWA. Bio-medical materials and engineering, 26(s1), S1791-S1796.

https://github.com/broadinstitute/gatk.









