#!/bin/bash
#Script used to remove ribosomal RNA still left in the samples after library preparation
#Made by Jonas N. Søndergaard
#Made on 181215

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A uppmax_proj_number
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 6:00:00
#SBATCH -J 181215_HiSAT2_rRNA
#SBATCH --output=181215_HiSAT2_rRNA.out
#SBATCH --error=181215_HiSAT2_rRNA.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load HISAT2/2.1.0 

#file paths
FQ_PATH=/proj/trimmedFQfiles
OUTPUT_PATH=/proj
REF_PATH=/proj/ref_genomes/rRNA.dir

#loop to run HiSAT2 alignment to rRNA, and keeping unaligned reads.
for i in {1..34}; do \
	FILE_NAME=`sed "${i}q;d" Name.list`
	
	hisat2 \
		-p 8 \
		--rna-strandness RF \
		--un-conc-gz ${OUTPUT_PATH}/FQfiles_without_rRNA/${FILE_NAME}_tc_rmrRNA.fastq.gz \
        	--summary-file ${OUTPUT_PATH}/results_rRNA_alignment/${FILE_NAME}_tc_rRNAcontam.txt \
		-x ${REF_PATH}/hrRNA_combined \
		-1 ${FQ_PATH}/${FILE_NAME}_tc_R1.fastq.gz \
		-2 ${FQ_PATH}/${FILE_NAME}_tc_R2.fastq.gz \
		-S ${OUTPUT_PATH}/results_rRNA_alignment/${FILE_NAME}_tc_rRNA.sam \
		>> ${OUTPUT_PATH}/results_rRNA_alignment/${FILE_NAME}_rRNAalign_stdout.stderr.txt 2>&1
done


#Readme:
#-p specifies the number of computational cores/threads that will be used by the program
#--rna-strandness: strand-specific information. Needs to be RF if using Illumina Truseq library preparation.
#--un-conc-gz: Write paired-end reads that fail to align concordantly to file(s) at <path>. Useful for rRNA removal
#--summary-file: Print alignment summary to this file
#-x path to the pre-built genome index. Note that the index consists of multiple files ending in .ht2 , and only the shared part of the filename should be indicated (e.g. genome if the files are called genome.1.ht2 , genome.2.ht2 , etc).
#hrRNA_combined is made by combining all RNA fasta files from the databases HGNC, ENA, and SILVA on 171009.
#-1 the first-read mate FASTQ file
#-2 the second-read mate FASTQ file
#-S name of the result file that will be created
#>> send all messages from HISAT2 (including errors and warnings) into the specified file
