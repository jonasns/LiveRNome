#!/bin/bash
#script to build ribosomal RNA reference genome for HiSAT2 
#Made by Jonas N. Søndergaard
#Made on 171003

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A uppmax_proj_number 
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 3:00:00

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load HISAT2/2.1.0 

#use HiSAT2 to build reference genome
hisat2-build \
	-f /proj/ref_genomes/rRNA.dir/hrRNA_combined.fa \
	/proj/ref_genomes/rRNA.dir/hrRNA_combined

#Readme
#hrRNA_combined.fa contains rRNA sequences from the HGNC, ENA, and SILVA databases, as well as some manually curated ones from Pubmed
