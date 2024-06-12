#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --partition=normal_q
#SBATCH --mem=240G
#SBATCH -N 1
#SBATCH --account=computeomics


cd ../data1
samples=`ls *_R1_001.fastq | awk '{split($_,x,"_R1_001.fastq"); print x[1]}' | sort | uniq`
cd ../breadth
cd influenza

for sample in $samples

do

bowtie2-build ../../influenza_gens/influenza.fasta refgenome

bowtie2 -p 12 -x refgenome -1 ../../data1/${sample}_R1_001.fastq -2 ../../data1/${sample}_R2_001.fastq -S ${sample}.sam

samtools view -S -b ${sample}.sam > ${sample}.bam

samtools sort ${sample}.bam -o ${sample}_sorted.bam

samtools index ${sample}_sorted.bam

samtools mpileup ${sample}_sorted.bam | wc -l

bowtie2-inspect -s refgenome | awk '{ FS = "\t" } ; BEGIN{L=0}; {L=L+$3}; END{print L}'

echo ${sample}" done"

done