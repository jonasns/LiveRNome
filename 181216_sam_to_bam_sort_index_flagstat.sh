#!/bin/bash
#script used for converting sam files to bam files and subsequently sort and index the bam files.
#Made by Jonas N. Søndergaard
#Made on 181216 

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A uppmax_proj_number
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH -J 181216_sam_to_bam_sort_index_flagstat
#SBATCH --output=181216_sam_to_bam_sort_index_flagstat.out
#SBATCH --error=181216_sam_to_bam_sort_index_flagstat.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load samtools/1.5

#file paths
SAM_DIR=/proj/SAMfiles_nsa
BAM_DIR=/proj/BAMfiles_nsa
BAM_DIR_SORTED=/proj/BAMfiles_sorted_nsa
FLAGSTAT_DIR=/proj/Flagstat_nsa

#loop to make BAM files from SAM files, and sorting, indexing and generating flagstats of the BAM files.
for i in {1..34}; do \
	FILE_NAME=`sed "${i}q;d" Name.list`

	samtools view \
		-bS ${SAM_DIR}/${FILE_NAME}_tc_rmrRNA.sam \
		> ${BAM_DIR}/${FILE_NAME}_tc_rmrRNA.bam \

	samtools sort \
		${BAM_DIR}/${FILE_NAME}_tc_rmrRNA.bam \
		-o ${BAM_DIR_SORTED}/${FILE_NAME}_tc_rmrRNA.sorted.bam

	samtools index \
		${BAM_DIR_SORTED}/${FILE_NAME}_tc_rmrRNA.sorted.bam

	samtools flagstat \
		${BAM_DIR_SORTED}/${FILE_NAME}_tc_rmrRNA.sorted.bam \
		> ${FLAGSTAT_DIR}/${FILE_NAME}.flagstat

	date
done
