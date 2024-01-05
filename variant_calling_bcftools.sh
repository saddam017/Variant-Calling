#!/bin/bash

echo "Run Prep files..."


#conda activate variants_calling_bcf
#download data
wget -p /home/saddam/Desktop/variant_calling_bcftools/reads http://ftp.sra.ebi.ac.uk/vol1/run/ERR333/ERR3335404/P7741_R1.fastq.gz
wget -p /home/saddam/Desktop/variant_calling_bcftools/reads http://ftp.sra.ebi.ac.uk/vol1/run/ERR333/ERR3335404/P7741_R2.fastq.gz


###################################################### VARIANT CALLING STEPS ##########################################


# directories


reads="/home/saddam/Desktop/variant_calling_bcftools/reads"
ref="/home/saddam/Desktop/variant_calling_bcftools/reference_genome/Agy99.fasta" #downloaded fro NCBI nucliotide Mycobacterium ucerans Agy99, complete genome
ref2="/home/saddam/Desktop/variant_calling_bcftools/reference_genome/Agy99.fasta.fai"
aligne_reads="/home/saddam/Desktop/variant_calling_bcftools/aligne_reads"
results="/home/saddam/Desktop/variant_calling_bcftools/results"

#-----------------
# STEP 1: QC - run fastqc
#-----------------------

echo "STEP 1: QC - run fastqc"

fastqc ${reads}/P7741_R1.fastq.gz -o ${reads}/
fastqc ${reads}/P7741_R2.fastq.gz -o ${reads}/

#trimmed bad quality read
echo "STEP 1: QC - run sickle for trimming"
sickle pe -f ${reads}/P7741_R1.fastq.gz -r ${reads}/P7741_R2.fastq.gz -t sanger -q 20 -l 20 -g -o ${reads}/trimmed_R1.fastq.gz -p ${reads}/trimmed_R2.fastq.gz -s ${reads}/trimmed_S.fastq.gz

echo "STEP 1: QC - run fastqc to check quality of trimmed reads"
fastqc ${reads}/trimmed_R1.fastq.gz -o ${reads}/
fastqc ${reads}/trimmed_R2.fastq.gz -o ${reads}/

#-----------------
# STEP 2: Map to reference using BWA-MEM
#-----------------------

echo "STEP 2: Map to reference using BWA-MEM"

# BWA index reference

bwa index ${ref}

# BWA allignment 
bwa mem -t 4 ${ref} ${reads}/trimmed_R1.fastq.gz ${reads}/trimmed_R2.fastq.gz > ${aligne_reads}/output.sam


#-----------------
# STEP 3: Create BAM file 
#-----------------------

echo "create BAM file from SAM file"
samtools view -S -b ${aligne_reads}/output.sam > ${aligne_reads}/output.bam


#-----------------
# STEP 4: shorted the BAM file 
#-----------------------

echo "STEP 4: shorted the BAM file " 
samtools sort -o ${aligne_reads}/output.sorted.bam ${aligne_reads}/output.bam

samtools flagstat ${aligne_reads}/output.sorted.bam > ${aligne_reads}/mappingstats.txt #to see the summary of bam file


#-----------------
# STEP 4: shorted the BAM file 
#-----------------------

echo "creating bcf fie"

samtools faidx ~/Desktop/variant_calling_bcftools/reference_genome/Agy99.fasta
bcftools mpileup -O b -o ${results}/raw.bcf -f ${ref} --threads 4 -q 20 -Q 30 ${aligne_reads}/output.sorted.bam



#-----------------
# STEP 5: Final steps Variants Calling 
#-----------------------
# variant calling code

bcftools call --ploidy 1 -m -v -o ${results}/variants.raw.bcf ${results}/raw.bcf

grep -v -c '^#' ${results}/variants.raw.bcf #to check the number of variants in the bcf file

bcftools view -v snps ${results}/variants.raw.bcf | grep -v -c '^#' # to check the number of snp in the bcf file
bcftools view -v indels ${results}/variants.raw.bcf | grep -v -c '^#' # to check the number of indels in the bcf file

bcftools query -f '%POS\n' ${results}/variants.raw.bcf >${results}/position.txt #positions of the variants


echo "variants calling is done by BCFTooLs"







