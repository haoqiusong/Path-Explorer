#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --partition=normal_q
#SBATCH --mem=240G
#SBATCH -N 1
#SBATCH --account=computeomics


cd ../data1
samples=`ls *_R1_001.fastq | awk '{split($_,x,"_R1_001.fastq"); print x[1]}' | sort | uniq`
cd ../example
cd norovirus

for sample in $samples 

do

bowtie2 -p 12 -x ../../norovirus/norovirus_ref -1 ../../data1/${sample}_R1_001.fastq -2 ../../data1/${sample}_R2_001.fastq -S ${sample}.sam

samtools view -S -b ${sample}.sam > ${sample}.bam

samtools sort ${sample}.bam -o ${sample}_sorted

bcftools mpileup -f *.fasta ${sample}_sorted | bcftools call -mv -Ov -o snp_${sample}.vcf

echo "done"

done