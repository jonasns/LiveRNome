#!/bin/bash
#Script used to align RNA-seq reads to a reference genome with HiSAT2
#Made by Jonas N. Søndergaard
#Made on 181215

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A uppmax_proj_number
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH -J 181215_HiSAT2_nsa_ref_genome
#SBATCH --output=181215_HiSAT2_nsa_ref_genome.out
#SBATCH --error=181215_HiSAT2_nsa_ref_genome.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load HISAT2/2.1.0 

#file paths
FQ_PATH=/proj/FQfiles_without_rRNA
OUTPUT_PATH=/proj/SAMfiles_nsa
REF_PATH=/proj/ref_genomes

#loop to run HiSAT2 alignment to a reference genome.
for i in {1..34}; do \
	FILE_NAME=`sed "${i}q;d" Name.list`

	hisat2 \
		-p 8 \
		--no-spliced-alignment \
		--rna-strandness RF \
		--dta \
		-k 5 \
       		--summary-file ${OUTPUT_PATH}/align_stats/${FILE_NAME}_tc_rmRNA_alignStats.txt \
		-x ${REF_PATH}/GRCh38.p10.genome \
		-1 ${FQ_PATH}/${FILE_NAME}_tc_rmrRNA.fastq.1.gz \
		-2 ${FQ_PATH}/${FILE_NAME}_tc_rmrRNA.fastq.2.gz \
		-S ${OUTPUT_PATH}/${FILE_NAME}_tc_rmrRNA.sam \
		>> ${OUTPUT_PATH}/align_stats/${FILE_NAME}_hg38align_stdout.stderr.txt 2>&1
done


#Readme:
#-p: specifies the number of computational cores/threads that will be used by the program
#--no-spliced-alignment: Disable spliced alignments (the algorithm does not attempt to look for reads that are split by overlapping an intron).
#--rna-strandness: strand-specific information. Needs to be RF if using Illumina Truseq library preparation.
#--dta: Report alignments tailored for transcript assemblers including StringTie
-k: number of accepted multi mappings.
#--summary-file: Print alignment summary to this file
#-x: path to the pre-built genome index. Note that the index consists of multiple files ending in .ht2 , and only the shared part of the filename should be indicated (e.g. genome if the files are called genome.1.ht2 , genome.2.ht2 , etc).
#-1: the first-read mate FASTQ file
#-2: the second-read mate FASTQ file
#-S: name of the result file that will be created
#--un-conc-gz: Write paired-end reads that fail to align concordantly to file(s) at <path>. Useful for rRNA removal
#>> send all messages from HISAT2 (including errors and warnings) into the specified file

