#!/bin/bash
#script made to count the number of reads aligned to genes and non-coding features in a reference genome file. 
#Made by Jonas N. Søndergaard
#Made on 190425

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A uppmax_proj_number
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 2:00:00
#SBATCH -J 190425_FeatureCounts
#SBATCH --output=190425_FeatureCounts.out
#SBATCH --error=190425_FeatureCounts.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load subread/1.5.2

#file paths
BAM_PATH=/proj/BAMfiles_sorted_nsa
OUTPUT_PATH=/proj/FeatureCounts_nsa
GTF_PATH=/proj/ref_genomes

#run featureCounts
featureCounts \
	-t exon \
	-g gene_id \
	-s 2 \
	-T 2 \
	-p \
	-M \
	-O \
	--fraction \
	-C \
	-a ${GTF_PATH}/190425_gencode.v27.annotation.2novel.noMYX.gtf \
	-o ${OUTPUT_PATH}/190425_nsa_2novel.readCount \
	${BAM_PATH}/*.sorted.bam \
	&> ${OUTPUT_PATH}/190425_nsa_2novel.readCount.log


#Readme
#-t: Specify feature type in GTF annotation. Here I chose exon, because on transcript level we might not get the lncRNAs
#-g: Specify attribute type in GTF annotation. Here we could chose e.g. transcript ID or gene ID. I chose gene ID, because I want to do DE analysis on gene level.
#-s: use '-s 2' if reversely stranded (as is the case for the Illumina Truseq library prep protocol)
#-T: Number of computational cores/threads used for the analysis
#-p: The experiment is paired end
#-M: Multi-mapping reads will also be counted. Each alignment will have 1 count or a fractional count if --fraction is specified
#-O: Allow reads that overlaps multiple features to be counted
#-C: If specified, the chimeric fragments (those fragments that have their two ends aligned to different chromosomes) will NOT be counted.
#-a: Name of the annotation file. Here it's gencode v27 chr1-22 (e.g. no M, X, or Y chromosomes). Additionally the two novel genes MSTRG.28468 and MSTRG.12891 has been added manually.
#-o: output file name