#!/bin/bash

#metagenomedata: /home/virus_i/NGS_data/DL/
 #e.g.           /home/virus_i/NGS_data/DL/SRR7755285/SRR7755285.sra
 
#HISAT2DB: /home/virus_i/NGS_DB/HISAT2/Zothers/
 #e.g.     /home/virus_i/NGS_DB/HISAT2/Zothers/GCF_001704415.1_ARS1_genomic.fna/genome_index.1.ht2
 
#OutputData: /home/virus_i/NGS_data/output/
 #e.g.       /home/virus_i/NGS_data/output/SRR7755285/

sra=$1
genome_fasta=$2
input="/home/virus_i/NGS_data/DL/"
output="/home/virus_i/NGS_data/output/"
hisat2DB="/home/virus_i/NGS_DB/HISAT2/Zothers/"${genome_fasta}


mkdir ${output}${sra}

#The steps of this group have been performed ahead of time.
	# 1. prepare the genome index 
	#hisat2-build \
	#	${input}${genome_fasta} \
	#	${output}${sra}/genome_index
	# 2. download the RNA-seq data
	#prefetch ${sra}

# convert sra to fastq
fasterq-dump --split-files ${input}/${sra}/${sra}.sra \
	--force -O ${output}${sra}

# 3. quality control of the RNA-seq data
fastp \
	-i  ${output}${sra}/${sra}_1.fastq \
	-I ${output}${sra}/${sra}_2.fastq \
	-o ${output}${sra}/${sra}_1_fastp.fastq \
	-O ${output}${sra}/${sra}_2_fastp.fastq \
	-l 35 -y -3 -W 3 -M 15 -x

# 4. mapping analysis to the host genome
mkdir ${output}${sra}/mapped

#hisat2 \
#	-x ${hisat2DB}/genome_index \
#	-1 ${output}${sra}/${sra}_1_fastp.fastq \
#	-2 ${output}${sra}/${sra}_2_fastp.fastq \
#	-S ${output}${sra}/mapped/Hisat2.out.sam

# 5. extract unmapped reads  
#samtools sort -O BAM -o ${output}${sra}/mapped/Hisat2.out.sorted.bam \
#	${output}${sra}/mapped/Hisat2.out.sam 
	
#samtools view -u -f 12 -F 256 \
#	${output}${sra}/mapped/Hisat2.out.sorted.bam \
#	> ${output}${sra}/mapped/Hisat2.out.sorted.unmapped.bam
    
#bedtools bamtofastq \
#	-i ${output}${sra}/mapped/Hisat2.out.sorted.unmapped.bam \
#	-fq ${output}${sra}/mapped/${sra}_R1_fastp.hisat.unmapped.fastq \
#	-fq2 ${output}${sra}/mapped/${sra}_R2_fastp.hisat.unmapped.fastq

# 4 & 5 were combined
hisat2 -p 16 \
    -x "${hisat2DB}/genome_index" \
    -1 "${output}${sra}/${sra}_1_fastp.fastq" \
    -2 "${output}${sra}/${sra}_2_fastp.fastq" \
| samtools view -u -f 12 -F 256 \
| samtools sort -n -O BAM \
| bedtools bamtofastq \
    -i - \
    -fq  "${output}${sra}/mapped/${sra}_R1_fastp.hisat.unmapped.fastq" \
    -fq2 "${output}${sra}/mapped/${sra}_R2_fastp.hisat.unmapped.fastq"
    
    
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

