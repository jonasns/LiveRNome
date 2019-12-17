#!/bin/bash
#script to run MultiQC
#Made by Jonas N. Søndergaard
#Made on 170928

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A uppmax_proj_number 
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 1:00:00

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load MultiQC/1.2

#run MultiQC
multiqc \
	/proj/FQfiles
