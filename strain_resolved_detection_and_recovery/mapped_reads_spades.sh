#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --partition=normal_q
#SBATCH --mem=240G
#SBATCH -N 1
#SBATCH --account=computeomics

#minimap2 -d reference_genome.mmi ref.fasta

samples=`ls /projects/ciwars/CDC-WBS/TestRun/RawData_LanesConcatenated/*_R1_001.fastq.gz | awk '{split($_,x,"_R1_001.fastq.gz"); print x[1]}' | sort | uniq`

for sample in $samples

do

name=$(echo ${sample} | cut -d'/' -f7)

#minimap2 -a reference_genome.mmi ${sample}_R1_001.fastq.gz > ${name}_mapped_reads_R1.sam
#minimap2 -a reference_genome.mmi ${sample}_R2_001.fastq.gz > ${name}_mapped_reads_R2.sam

#samtools view -bS ${name}_mapped_reads_R1.sam > ${name}_mapped_reads_R1.bam
#samtools view -bS ${name}_mapped_reads_R2.sam > ${name}_mapped_reads_R2.bam

#samtools sort ${name}_mapped_reads_R1.bam -o sorted_${name}_mapped_reads_R1.bam
#samtools sort ${name}_mapped_reads_R2.bam -o sorted_${name}_mapped_reads_R2.bam

#samtools index sorted_${name}_mapped_reads_R1.bam
#samtools index sorted_${name}_mapped_reads_R2.bam

#samtools view -b -F 4 sorted_${name}_mapped_reads_R1.bam > ${name}_mapped_reads_only_R1.bam
#samtools view -b -F 4 sorted_${name}_mapped_reads_R2.bam > ${name}_mapped_reads_only_R2.bam

#samtools index ${name}_mapped_reads_only_R1.bam
#samtools index ${name}_mapped_reads_only_R2.bam

#samtools fastq -F 4 ${name}_mapped_reads_only_R1.bam > ${name}_R1.fastq
#samtools fastq -F 4 ${name}_mapped_reads_only_R2.bam > ${name}_R2.fastq

cat ${name}_R1.fastq ${name}_R2.fastq > ${name}_merged.fastq

mkdir ${name}
spades.py -s ${name}_merged.fastq -o ${name}

echo ${name}

done