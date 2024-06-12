#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --partition=normal_q
#SBATCH --mem=60G
#SBATCH -N 1
#SBATCH --account=computeomics

for i in {1..14}; do
    ./themisto pseudoalign -q /projects/ciwars/haoqiu_all/cdc/extraction/8509-S${i}_S${i}_CAT_R1.fastq.paired.fq -o pseudoalignments_1.aln -i themisto_index/index --temp-dir themisto_index/tmp --rc --n-threads 16 --sort-output --gzip-output

    ./themisto pseudoalign -q /projects/ciwars/haoqiu_all/cdc/extraction/8509-S${i}_S${i}_CAT_R2.fastq.paired.fq -o pseudoalignments_2.aln -i themisto_index/index --temp-dir themisto_index/tmp --rc --n-threads 16 --sort-output --gzip-output

    mSWEEP --themisto-1 pseudoalignments_1.aln.gz --themisto-2 pseudoalignments_2.aln.gz -o mSWEEP -i reference_grouping.txt --write-probs

    mGEMS -r /projects/ciwars/haoqiu_all/cdc/extraction/8509-S${i}_S${i}_CAT_R1.fastq.paired.fq,/projects/ciwars/haoqiu_all/cdc/extraction/8509-S${i}_S${i}_CAT_R2.fastq.paired.fq -i reference_grouping.txt --themisto-alns pseudoalignments_1.aln.gz,pseudoalignments_2.aln.gz -o mGEMS-out --probs mSWEEP_probs.tsv -a mSWEEP_abundances.txt --index themisto_index

    cp mSWEEP_abundances.txt $i.txt

    echo $i
done

echo "Done"