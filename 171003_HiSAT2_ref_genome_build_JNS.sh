#!/bin/bash
#script to build reference genome for HiSAT2 
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
	-f /proj/ref_genomes/GRCh38.p10.genome.fa \
	/proj/ref_genomes/GRCh38.p10.genome
