#!/bin/bash
#script used to make bedgraph files using bedtools 
#Made by Jonas N. Søndergaard
#Made on 171130

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A uppmax_proj_number
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -J 171130_bedgraphGen_bedtools_fwdrev_merged_JNS
#SBATCH --output=171130_bedgraphGen_bedtools_fwdrev_merged_JNS.out
#SBATCH --error=171130_bedgraphGen_bedtools_fwdrev_merged_JNS.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load BEDTools/2.26.0

#file paths
OUTPUT_PATH=/proj/bedGraph_nsa_merged
BAM_PATH=/proj/splitBams

#loop to run make bedgraph files for 11 BAM files
for i in {1..11}; do \
	FILE_NAME=`sed "${i}q;d" siRNA.list`
	
	bedtools genomecov \
		-ibam ${BAM_PATH}/${FILE_NAME}*.fwd.sorted.bam \
		-bg \
		> ${OUTPUT_PATH}/${FILE_NAME}.fwd.bg	

        bedtools genomecov \
                -ibam ${BAM_PATH}/${FILE_NAME}*.rev.sorted.bam \
                -bg \
                > ${OUTPUT_PATH}/${FILE_NAME}.rev.bg

	#keep only reads aligning to chromosome 1-22
	awk '$1~/chr[1-9]/ || $1~/#/' \
		${OUTPUT_PATH}/${FILE_NAME}.fwd.bg \
		> ${OUTPUT_PATH}/${FILE_NAME}.fwd.noMYX.bg

        awk '$1~/chr[1-9]/ || $1~/#/' \
                ${OUTPUT_PATH}/${FILE_NAME}.rev.bg \
                > ${OUTPUT_PATH}/${FILE_NAME}.rev.noMYX.bg

	#add header to the bedgraph file 
	echo 'track type=bedGraph visibility=full color=0,0,255 name="'${FILE_NAME}'.fwd"' \
                | cat - ${OUTPUT_PATH}/${FILE_NAME}.fwd.noMYX.bg \
                > ${FILE_NAME}.temp \
                && mv ${FILE_NAME}.temp ${OUTPUT_PATH}/${FILE_NAME}.fwd.noMYX.bg \

	echo 'track type=bedGraph visibility=full color=255,0,0 name="'${FILE_NAME}'.rev"' \
                | cat - ${OUTPUT_PATH}/${FILE_NAME}.rev.noMYX.bg \
                > ${FILE_NAME}.temp \
                && mv ${FILE_NAME}.temp ${OUTPUT_PATH}/${FILE_NAME}.rev.noMYX.bg \
	
	#gzip bedgraph files	
	gzip ${OUTPUT_PATH}/${FILE_NAME}.*.noMYX.bg

done
