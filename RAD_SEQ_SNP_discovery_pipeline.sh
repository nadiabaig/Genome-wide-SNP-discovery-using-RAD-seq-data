#!/bin/bash

#PBS -l select=1:ncpus=6:mem=50G
#PBS -l walltime=20:00:00
#PBS -A  "yourprojectname_on_hpc_cluster"


module load gcc/8.1.0    #Stacks pre-requisite
module load Java/1.8.0  #Pre-requisite of fastqc and GATK
module load BWA/0.7.15 # needed to do indexing and alignment
module load SamTools/1.6 # needed to convert sam into bam and other statistics
module load Stacks/2.5


pw_dir="/gpfs/project/projects/-----/maze_pipeline"   #path to files where processed data will be saved
lib="/gpfs/project/projects/----/data_lib" # raw data
cd $pw_dir
mkdir data_1   
cd "$pw_dir/data_1"
mkdir p_tags  #for storing processed rad-tags

#####################################################
#Step-01: Data de-multiplexing-RAD-tag processing   #
#####################################################
st_tool="/gpfs/project/projects/..../src/stacks-2.5/process_radtags"  #if tool isn't available on your hpc cluster install it in your src folder
path_rads="$pw_dir/data_1/p_tags"
cd $path_rads

barcode_file="$pw_dir/barcodes_new.csv" # change the path for barcode files to process other datasets
for j in 1 ; #repeat the same for other folders,just change the last index of folder name, e.g data_lib2 then do for j in 2 and so on.
do
lib_new=$lib$j
  $st_tool -p ${lib_new} -o $path_rads -b $barcode_file -e kpnI -r -c -D -i gzfastq --retain_header --barcode_dist_1 2 -q --filter_illumina
done

wait
#################################################
#step-02:Check the quality of reads using fastqc#
################################################
cd "$pw_dir/data1"
mkdir quality
cd "$pw_dir/data1/quality"
Fastqc_tool="/gpfs/project/src/FastQC/fastqc" #install it, if not present on HPC
for k in $( ls $path_rads/*.gz)
do
$Fastqc_tool ${k} -t 10 --extract
done

wait

####################################
#Step-03:Trimming low quality reads#
####################################
Trimmomatic="/gpfs/project/projects/----/src/Trimmomatic-0.39/trimmomatic-0.39.jar"
adapter="/gpfs/project/projects/---/src/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"  #choose type of adapter according to your settings of reads, in my case I used paired-end reads to test the pipeline
trimmer="ILLUMINACLIP:$adapter:2:30:10 SLIDINGWINDOW:4:15 MINLEN:75"
folder= "$pw_dir/data1"
cd $folder
mkdir trimmed_reads
cd "$folder/trimmed_reads"
for  files in $(ls $path_rads/*.fq.gz); #we have processed rad_seq data here..step1 results
do
java -jar $Trimmomatic SE -threads 4 -phred33 ${files} "`basename $files .fq.gz`.trimmomatic_out.fq.gz"  $trimmer &>> output.file
done

wait
#################################################
#step-04: Reference genome indexing using bwa   # 
#################################################
Ref="/gpfs/project/projects/----/maze_pipeline/reference" # download your own reference. I have mentioned link in .README file
cd $Ref

bwa index $Ref/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz

wait

###############################################
#step-05: Alignment                           #
###############################################
trimmed_path=$folder/trimmed_reads
Index=$Ref/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz
cd $folder  #path where datalib1 folder is present
mkdir alignments
aln="$folder/alignments" #output file to save sam/bam files
cd $aln
for line in $(ls $trimmed_path/*.fq.gz | rev | cut -c 23-| rev | uniq)
do
bwa mem -M -t 2 -R "@RG\tID:${line}\tPU:${line}$\tSM:${line}\tPL:ILLUMINA\tLB:${line}" $Index ${line}.trimmomatic_out.fq.gz > $line.sam
done

wait
#######################################################
#step:06- conversion to bam ,sorting and indexing      #
########################################################
cd $aln
for files in $(ls $aln/*.sam | rev | cut -c 5-| rev | uniq)
do
samtools view -S -b ${files}.sam > ${files}.bam
wait
samtools sort  ${files}.bam -o ${files}_sorted.bam
wait
samtools index ${files}_sorted.bam
done

wait
######or alternatively you can also use this script for your step6...
#or you can also add this code
#!/bin/bash

#PBS -l select=1:ncpus=1:mem=1G
#PBS -l walltime=00:30:00
#PBS -A  "xyz"

path="/gpfs/project/----/alignments"
for x in {1,2,3,4}; #your files where alignment is stored for each sample e.g your file names are aln1,aln2,aln3
do

cd /gpfs/project/----/alignments/$x/

echo "#!/bin/bash
        #PBS -l select=1:ncpus=1:mem=2gb
        #PBS -l walltime=3:59:00
        #PBS -A  ''

        module load SamTools/1.6 
        
        cd /gpfs/project/----/alignments/$x/ 
         
        for i in \$(ls *.bam);
         do
  
          samtools view -h -q 20 -b \${i} > filtered_\${i}
          wait
          samtools index filtered_\${i}
          wait
          samtools stats filtered_\${i} > filtered_stats_\${i}
   
     
  
        done" > /gpfs/project/---/mq.$x.sh
cd /gpfs/project/..../alignments/
qsub mq.$x.sh
done

######################################
#creating sequence dictionary of reference #
#####################################
cd $Ref  #path where reference sequence is stored

java -jar /software/picard/1.130/picard.jar CreateSequenceDictionary \
REFERENCE=$Ref/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz \
OUTPUT=Zea_mays.B73_RefGen_v4.dna.toplevel.dict

wait
#####################################
#step:07-variant calling using GATK #
#####################################
cd $aln #alignment folder to store vcf
java -jar /software/gatk/3.7/GenomeAnalysisTK.jar \
 -T HaplotypeCaller -R $Ref/Zea_mays.B73_RefGen_v4.dna.toplevel.fa \
 -I list1.list \ #create a list of proceesed bam files to do variant calling collectively for all samples
 --genotyping_mode 'DISCOVERY'  -ploidy 2  -mbq 30 -mmq 40 -stand_call_conf 30 \
 -o raw_variants.vcf

#####################################
#step:07-Variant filtering          #
#####################################
cd $aln
bcf="/gpfs/project/projects/----/src/bcftools-1.3.1/bcftools"
$bcf view -m2 -M2 -v snps  $aln/raw_variants.vcf > bi_allelic.vcf #using bcf tools as vcftools --minalleles 2 parameter is not working. You can also write a python script to to filter

#########################################
#other options for biallelic filetering #
#########################################
#!/bin/bash

#PBS -l select=1:ncpus=1:mem=1G
#PBS -l walltime=1:00:00
#PBS -A  "xyz"


module load GATK/3.7
module load SamTools/1.6
module load Java/1.8.0_151

cd /gpfs/project/projects/-/raw_vcf  #change paths according to your folders settings
reference=".../.../.fasta"
for x in {1,2,3,4};
do
cd /gpfs/project/---/raw_vcf/$x
echo "#!/bin/bash
#PBS -l select=1:ncpus=2:mem=10G
#PBS -l walltime=36:59:00
#PBS -A "xyz"

module load GATK/3.7
module load SamTools/1.6
module load Java/1.8.0_151

cd /gpfs/project/..../raw_vcf/$x

for files in \$(ls *.vcf.gz)
do
java -jar /software/gatk/3.7/GenomeAnalysisTK.jar -T SelectVariants \
  -R $reference \
  -o \$files._biallelic.vcf \
  --variant \$files \
  -restrictAllelesTo BIALLELIC
done " > /gpfs/project/..../raw_vcf/bi-allelic.$x.sh
cd /gpfs/project/..../raw_vcf
qsub bi-allelic.$x.sh

done
wait

#set filters
MAF=0.05
Depth=4.0
quality=40.0
cd $aln
v_tool="/gpfs/project/projects/src/vcftools_0.1.13/bin/vcftools"

# perform the filtering with vcftools 
$v_tool --vcf $aln/bi_allelic.vcf --out Filtered.vcf \
 --maf $MAF \
 --minDP $Depth \
 --minQ $quality  --recode
path="/gpfs/project/-----/RAD-Seq-Maize/SNPs_final/maxDP50_maf0.05_MaxMis80_max95missSamples_QD_2_MQ_40_hard_filterd_SNPs.recode.vcf"
$v_tool stats $path



               
