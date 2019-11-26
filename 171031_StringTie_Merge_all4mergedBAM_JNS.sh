#!/bin/bash
#script for merging gtf files and cleaning up afterwards 
#Made by Jonas N. Søndergaard
#Made on 171031

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A snic2017-7-154 
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 3:00:00
#SBATCH -J 171031_StringTie_Merge_all4mergedBAM_JNS
#SBATCH --output=171031_StringTie_Merge_all4mergedBAM_JNS.out
#SBATCH --error=171031_StringTie_Merge_all4mergedBAM_JNS.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load StringTie/1.3.3

#file paths
OUTPUT_PATH=/proj/StringTie
REF_PATH=/proj/ref_genomes

stringtie \
	--merge \
	-p 8 \
	-G ${REF_PATH}/gencode.v27.annotation.gtf \
	-o ${OUTPUT_PATH}/stringtie_merged.gtf \
	${OUTPUT_PATH}/mergelist.txt

#Remove header and select chr1-22: 
tail -n+3 stringtie_merged.gtf | awk '$1~/chr[1-9]/' > tmp.gtf 

#Clean up with Python script made by Hassan Foroughi Asl:
module load python/3.6.0  

python3.6 /proj/9.PrepareMergedGTF.denovoStringTie.py tmp.gtf 171031_gencode.v27.annotation.novel.noMYX.gtf   

#Readme:
#Here mergelist.txt is a text file that has the names of the gene transfer format (GTF) files created in the previous step, with each file name on a single line. If you do not run StringTie in the same directory in which all the GTF files are, then you also need to include the full path in each GTF file name in mergelist.txt
