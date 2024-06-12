#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --partition=normal_q
#SBATCH --mem=120G
#SBATCH -N 1
#SBATCH --account=computeomics


cat /projects/ciwars/CDC-WBS/TestRun/RawData_LanesConcatenated/CDC_Cocktail_zym_S82_CAT_R1_001.fastq.gz /projects/ciwars/CDC-WBS/TestRun/RawData_LanesConcatenated/CDC_Cocktail_zym_S82_CAT_R2_001.fastq.gz > cocktail_concat.fastq.gz

sourmash sketch dna cocktail_concat.fastq.gz -o cocktail.sig

samples=`ls /projects/ciwars/CDC-WBS/TestRun/RawData_LanesConcatenated/*_R1_001.fastq.gz | awk '{split($_,x,"_R1_001.fastq.gz"); print x[1]}' | sort | uniq`

for sample in $samples

do

name=$(echo ${sample} | cut -d'/' -f7)

cat ${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz > ${name}_combined.fastq.gz

sourmash sketch dna ${name}_combined.fastq.gz -o ${name}.sig

echo ${name}

sourmash search ${name}.sig cocktail.sig --containment -t 0
sourmash search cocktail.sig ${name}.sig --containment -t 0

done