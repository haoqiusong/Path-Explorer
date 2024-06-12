#!/bin/bash
#SBATCH --time=30:00:00
#SBATCH --partition=normal_q
#SBATCH --mem=240G
#SBATCH -N 1
#SBATCH --account=computeomics

minimap2 -d reference.mmi ../ref_genomes/Escherichia_coli_ATCC_11775.fasta

samples=`ls /projects/ciwars/haoqiu_all/cdc/sourmash/*_CAT_combined.fastq.gz | awk '{split($_,x,"_CAT_combined.fastq.gz"); print x[1]}' | sort | uniq`

for sample in $samples

do

name=$(echo ${sample} | cut -d'/' -f7)

minimap2 -ax map-ont reference.mmi ${sample}_CAT_combined.fastq.gz > ${name}_alignments.sam

samtools view -bS ${name}_alignments.sam | samtools sort -@8 -o ${name}_alignments.sorted.bam -

samtools depth ${name}_alignments.sorted.bam > ${name}_coverage.txt

done