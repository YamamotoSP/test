#!/bin/bash

#metagenomedata: /home/virus_i/NGS_data/DL
 #e.g.           /home/virus_i/NGS_data/DL/SRR7755285/SRR7755285.sra
 
#HISAT2DB: /home/virus_i/NGS_DB/HISAT2/Zothers
 #e.g.     /home/virus_i/NGS_DB/HISAT2/Zothers/GCF_001704415.1_ARS1_genomic.fna/genome_index.1.ht2
 
#OutputData: /home/virus_i/NGS_data/output
 #e.g.       /home/virus_i/NGS_data/output/SRR7755285

cd /home/virus_i/Github/test

sra="SRR7755285" 
genome_fasta="GCF_001704415.1_ARS1_genomic.fna"

bash A1_Virome_seq_assembly.sh \
	${sra}  \
	${genome_fasta}  

