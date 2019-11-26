#!/bin/bash
#script for generating bedgraph files using the homer package
#Made by Jonas N. Søndergaard
#Made on 181218

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A uppmax_proj_number 
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 10:00:00
#SBATCH -J 181218_bedgraphGenIndiv_homer_fwdrev
#SBATCH --output=181218_bedgraphGenIndiv_homer_fwdrev.out
#SBATCH --error=181218_bedgraphGenIndiv_homer_fwdrev.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load samtools/1.5

#file paths
OUTPUT_PATH=/proj/sllstore2017022/EXP_18_CA6218/bedGraph_nsa_individual/TagDirectories
BAM_PATH=/proj/sllstore2017022/EXP_18_CA6218/splitBams

#loop to generate bedgraph files for 34 files
for i in {1..34}; do \
	FILE_NAME=`sed "${i}q;d" Name.list`
	
#make TagDirectory for each bam file:
	/home/homer/bin/makeTagDirectory \
	${OUTPUT_PATH}/${FILE_NAME}.fwd \
	${BAM_PATH}/${FILE_NAME}*.fwd.sorted.bam

        /home/homer/bin/makeTagDirectory \
        ${OUTPUT_PATH}/${FILE_NAME}.rev \
        ${BAM_PATH}/${FILE_NAME}*.rev.sorted.bam

#remove reads outside Chr1-22:
	/home/homer/bin/removeOutOfBoundsReads.pl \
	${OUTPUT_PATH}/${FILE_NAME}.fwd \
	hg38

        /home/homer/bin/removeOutOfBoundsReads.pl \
        ${OUTPUT_PATH}/${FILE_NAME}.rev \
        hg38

#make bedgraph file:
	/home/homer/bin/makeUCSCfile \
	${OUTPUT_PATH}/${FILE_NAME}.fwd \
	-o auto \
        -norm 1e7 \
	-normLength 100 \
	-style rnaseq \
	-strand both \
	-color 0,0,255

        /home/homer/bin/makeUCSCfile \
        ${OUTPUT_PATH}/${FILE_NAME}.rev \
        -o auto \
	-norm 1e7 \
        -normLength 100 \
        -style rnaseq \
	-strand both \
	-color 255,0,0

done
