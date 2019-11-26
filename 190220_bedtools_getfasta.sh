#!/bin/bash
#script used to getting a fasta file using a bed file
#Made by Jonas N. Søndergaard
#Made on 190220

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A uppmax_proj_number
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J 190220_bedtools_getfasta
#SBATCH --output=190220_bedtools_getfasta.out
#SBATCH --error=190220_bedtools_getfasta.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load BEDTools/2.27.1

#file paths
BED_PATH=/proj/promoter_BED
OUTPUT_PATH=/proj/promoter_fasta
GTF_PATH=/proj/ref_genomes

#run bedtools
bedtools getfasta \
	-fi ${GTF_PATH}/GRCh38.p10.genome.fa \
	-bed ${BED_PATH}/190220_lncFARED1_prom_all.bed \
	-fo ${OUTPUT_PATH}/190220_lncFARED1_prom_all.fa