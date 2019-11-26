#!/bin/bash
#script to run FastQC on fastq files
#Made by Jonas N. Søndergaard
#Made on 181214

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A uppmax_proj_number
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 6:00:00
#SBATCH -J 181214_FastQC
#SBATCH --output=181214_FastQC.out
#SBATCH --error=181214_FastQC.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load FastQC/0.11.5

#run FastQC
fastqc \
	/proj/FQfiles/CK* \
	--outdir /proj/FastQC
