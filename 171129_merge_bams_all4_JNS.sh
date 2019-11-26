#!/bin/bash
#script to merge bam files
#Made by Jonas N. Søndergaard
#Made on 171129

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A uppmax_proj_number
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 8:00:00
#SBATCH -J 171129_merge_bams
#SBATCH --output=171129_merge_bams.out
#SBATCH --error=171129_merge_bams.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages 
module load bioinfo-tools
module load samtools/1.5

#file paths
BAM_PATH=/proj/BAMfiles_sorted_nsa
OUTPUT_PATH=/proj/Merged_bam_all4

#loop to merge 4 bam files into 1 for 11 conditions. And subsequent sorting and indexing. 
for i in {1..11}; do \
	FILE_NAME=`sed "${i}q;d" siRNA.list`

	samtools merge \
		${OUTPUT_PATH}/${FILE_NAME}_nsa_all4.bam \
		${BAM_PATH}/*${FILE_NAME}*.bam

	samtools sort \
		${OUTPUT_PATH}/${FILE_NAME}_nsa_all4.bam \
		-o ${OUTPUT_PATH}/${FILE_NAME}_nsa_all4.sorted.bam

	samtools index \
		${OUTPUT_PATH}/${FILE_NAME}_nsa_all4.sorted.bam

done

