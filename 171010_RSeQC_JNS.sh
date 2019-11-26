#!/bin/bash
#script to run RseQC (RNA-seq Quality Control Package)
#Made by Jonas N. Søndergaard
#Made on 171010

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A uppmax_proj_number
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 20:00:00
#SBATCH -J 171010_RSeQC
#SBATCH --output=171010_RSeQC.out
#SBATCH --error=171010_RSeQC.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load rseqc/2.6.4

#file paths
BAM_PATH=/proj/BAMfiles_sorted
BED_PATH=/proj/ref_genomes
OUTPUT_PATH=/proj/RseQC

#Loop to run RseQC (RNA-seq Quality Control Package)
for i in {56..99}; do \
	FILE_NAME=`sed "${i}q;d" Name.list.rev`

	read_distribution.py \
		-r ${BED_PATH}/gencode.v27.annotation.noMYX.bed \
		-i ${BAM_PATH}/CK00${i}*.sorted.bam \
		> ${OUTPUT_PATH}/read_distribution/${FILE_NAME}_read_distribution.out

	geneBody_coverage.py \
		-r ${BED_PATH}/hg38.HouseKeepingGenes.bed \
		-i ${BAM_PATH}/CK00${i}*.sorted.bam \
		-o ${OUTPUT_PATH}/geneBody_coverage/${FILE_NAME}_geneBody_coverage.out

	inner_distance.py \
		-r ${BED_PATH}/gencode.v27.annotation.noMYX.bed \
		-i ${BAM_PATH}/CK00${i}*.sorted.bam \
		-o ${OUTPUT_PATH}/inner_distance/${FILE_NAME}_inner_distance.out

	read_duplication.py \
		-i ${BAM_PATH}/CK00${i}*.sorted.bam \
		-o ${OUTPUT_PATH}/read_duplication/${FILE_NAME}_read_duplication.out

	read_GC.py \
		-i ${BAM_PATH}/CK00${i}*.sorted.bam \
		-o ${OUTPUT_PATH}/read_GC/${FILE_NAME}_read_GC.out

	read_NVC.py \
		-i ${BAM_PATH}/CK00${i}*.sorted.bam \
		-o ${OUTPUT_PATH}/read_NVC/${FILE_NAME}_read_NVC.out

	read_quality.py \
                -i ${BAM_PATH}/CK00${i}*.sorted.bam \
                -o ${OUTPUT_PATH}/read_quality/${FILE_NAME}_read_quality.out

	date
done
