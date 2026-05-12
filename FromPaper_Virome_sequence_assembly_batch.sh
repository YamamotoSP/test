#!/bin/bash

sra=$1
genome_fasta=$2
input="input_examples/"
output="output_examples/"

mkdir ${input}
mkdir ${output}

# 1. prepare the genome index 
mkdir ${input}${sra}

#hisat2-build \
#	${input}${genome_fasta} \
#	${output}${sra}/genome_index

# 2. download the RNA-seq data
prefetch ${sra}

# convert sra to fastq
fasterq-dump --split-files ${sra}/${sra}.sra \
	--force -O ${output}${sra}

# 3. quality control of the RNA-seq data
fastp \
	-i  ${output}${sra}/${sra}.sra_1.fastq \
	-I ${output}${sra}/${sra}.sra_2.fastq \
	-o ${output}${sra}/${sra}.sra_1_fastp.fastq \
	-O ${output}${sra}/${sra}.sra_2_fastp.fastq \
	-l 35 -y -3 -W 3 -M 15 -x

# 4. mapping analysis to the host genome
mkdir ${output}${sra}/mapped

hisat2 \
	-x ${output}${sra}/genome_index \
	-1 ${output}${sra}/${sra}.sra_1_fastp.fastq \
	-2 ${output}${sra}/${sra}.sra_2_fastp.fastq \
	-S ${output}${sra}/mapped/Hisat2.out.sam

# 5. extract unmapped reads  
samtools sort -O BAM -o ${output}${sra}/mapped/Hisat2.out.sorted.bam \
	${output}${sra}/mapped/Hisat2.out.sam 
	
samtools view -u -f 12 -F 256 \
	${output}${sra}/mapped/Hisat2.out.sorted.bam \
	> ${output}${sra}/mapped/Hisat2.out.sorted.unmapped.bam

bedtools bamtofastq \
	-i ${output}${sra}/mapped/Hisat2.out.sorted.unmapped.bam \
	-fq ${output}${sra}/mapped/${sra}_R1_fastp.hisat.unmapped.fastq \
	-fq2 ${output}${sra}/mapped/${sra}_R2_fastp.hisat.unmapped.fastq

# 6. sequence assembly
mkdir ${output}${sra}/assembly
spades.py \
	-1 ${output}${sra}/mapped/${sra}_R1_fastp.hisat.unmapped.fastq  \
	-2 ${output}${sra}/mapped/${sra}_R2_fastp.hisat.unmapped.fastq \
	-k 21,33,55,77,99 \
	-o ${output}${sra}/assembly/spades/

metaspades.py \
	-1 ${output}${sra}/mapped/${sra}_R1_fastp.hisat.unmapped.fastq  \
	-2 ${output}${sra}/mapped/${sra}_R2_fastp.hisat.unmapped.fastq \
	-k 21,33,55,77,99 \
	-o ${output}${sra}/assembly/metaspades

# 7. clustering sequences
cat \
	${output}${sra}/assembly/spades/scaffolds.fasta \
	${output}${sra}/assembly/metaspades/scaffolds.fasta \
	> ${output}${sra}/assembly/scaffolds.fasta 

seqkit seq -m 500 \
	${output}${sra}/assembly/scaffolds.fasta \
	> ${output}${sra}/assembly/scaffolds.500.fasta 

cd-hit-est \
	-i  ${output}${sra}/assembly/scaffolds.500.fasta \
	-o  ${output}${sra}_scaffolds.fasta \
	-c 0.95 -n 10
