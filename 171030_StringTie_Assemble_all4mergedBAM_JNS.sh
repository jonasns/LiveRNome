#!/bin/bash
#script for assembling novel GTF files from BAM files
#Made by Jonas N. Søndergaard
#Made on 171030

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A uppmax_proj_number
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 4:00:00
#SBATCH -J 171030_StringTie_Assemble
#SBATCH --output=171030_StringTie_Assemble.out
#SBATCH --error=171030_StringTie_Assemble.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load StringTie/1.3.3

#file paths
BAM_PATH=/proj/BAMfiles_sorted
OUTPUT_PATH=/proj/StringTie/
REF_PATH=/proj/ref_genomes

#loop to assemble novel GTFs for 11 files
for i in {1..11}; do \
	FILE_NAME=`sed "${i}q;d" siRNA.list`
	
	stringtie \
		-p 8 \
		--rf \
		-G ${REF_PATH}/gencode.v27.annotation.gtf \
		-o ${OUTPUT_PATH}/${FILE_NAME}.gtf \
		${BAM_PATH}/${FILE_NAME}*.bam 
done

#Readme:
#-p: number of computational cores used to run the script
#--rf: Assumes a stranded library fr-firststrand as found with Illumina Truseq library prep protocol.
#-G: known annotations
#-o: output file name
