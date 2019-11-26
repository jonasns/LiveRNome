#!/bin/bash
#script for splitting bam files strand-wise
#Made by Jonas N. Søndergaard
#Made on 181218
          
#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A uppmax_proj_number
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 18:00:00
#SBATCH -J 181218_splitBam_fwdrev
#SBATCH --output=181218_splitBam_fwdrev.out
#SBATCH --error=181218_splitBam_fwdrev.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load samtools/1.5

#file paths
OUTPUT_PATH=/proj/splitBams
BAM_PATH=/proj/BAMfiles_sorted_nsa

#loop to split BAM files into forward and reverse strand for 34 files
for i in {1..34}; do \
        FILE_NAME=`sed "${i}q;d" Name.list`

	# include reads that are 2nd in a pair (128);
	# exclude reads that are mapped to the reverse strand (16)
	samtools view -b -f 128 -F 16 ${BAM_PATH}/$FILE_NAME*.bam > ${OUTPUT_PATH}/a.fwd1.bam

	# exclude reads that are mapped to the reverse strand (16) and
	# first in a pair (64): 64 + 16 = 80
	samtools view -b -f 80 ${BAM_PATH}/$FILE_NAME*.bam > ${OUTPUT_PATH}/a.fwd2.bam

	# combine the temporary files
	samtools merge -f ${OUTPUT_PATH}/$FILE_NAME.fwd.bam ${OUTPUT_PATH}/a.fwd1.bam ${OUTPUT_PATH}/a.fwd2.bam

	#sort the filtered BAM file
	samtools sort ${OUTPUT_PATH}/$FILE_NAME.fwd.bam -o ${OUTPUT_PATH}/$FILE_NAME.fwd.sorted.bam

	# index the filtered BAM file
	samtools index ${OUTPUT_PATH}/$FILE_NAME.fwd.sorted.bam

	# remove the temporary files
	rm ${OUTPUT_PATH}/a.fwd*.bam

	# include reads that map to the reverse strand (128)
	# and are second in a pair (16): 128 + 16 = 144
	samtools view -b -f 144 ${BAM_PATH}/$FILE_NAME*.bam > ${OUTPUT_PATH}/a.rev1.bam

	# include reads that are first in a pair (64), but
	# exclude those ones that map to the reverse strand (16)
	samtools view -b -f 64 -F 16 ${BAM_PATH}/$FILE_NAME*.bam > ${OUTPUT_PATH}/a.rev2.bam

	# merge the temporary files
	samtools merge -f ${OUTPUT_PATH}/$FILE_NAME.rev.bam ${OUTPUT_PATH}/a.rev1.bam ${OUTPUT_PATH}/a.rev2.bam

	#sort the filtered BAM file
	samtools sort ${OUTPUT_PATH}/$FILE_NAME.rev.bam -o ${OUTPUT_PATH}/$FILE_NAME.rev.sorted.bam

	# index the merged, filtered BAM file
	samtools index ${OUTPUT_PATH}/$FILE_NAME.rev.sorted.bam

	# remove temporary files
	rm ${OUTPUT_PATH}/a.rev*.bam

done

#Readme:
#https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html
#https://www.biostars.org/p/92935/
