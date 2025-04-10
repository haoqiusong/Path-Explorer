#!/usr/bin/env nextflow

// Enable DSL 2 syntax
nextflow.enable.dsl = 2

params.R1 = "Read1_file"
params.R2 = "Read2_file"
params.Ref = "Ref_file"

params.out_fname = "output"

process QC {

  input:
  path R1
  path R2

  output:
  path "${params.out_fname}_qc_R1.fastq.gz", emit: qc_R1
  path "${params.out_fname}_qc_R2.fastq.gz", emit: qc_R2

  """
  $projectDir/fastp -i $R1 -I $R2 -o ${params.out_fname}_qc_R1.fastq.gz -O ${params.out_fname}_qc_R2.fastq.gz --detect_adapter_for_pe --trim_poly_g --trim_poly_x --low_complexity_filter --average_qual 10 --thread 8
  """

}

process read_extraction {

  input:
  path qc_R1
  path qc_R2
  path Ref

  output:
  path "merged.fastq", emit: mapped_read
  path "R1.fastq", emit: extracted_R1
  path "R2.fastq", emit: extracted_R2

  """
    minimap2 -d reference_genome.mmi $Ref
    minimap2 -a reference_genome.mmi $qc_R1 > mapped_reads_R1.sam
    minimap2 -a reference_genome.mmi $qc_R2 > mapped_reads_R2.sam
    samtools view -bS mapped_reads_R1.sam > mapped_reads_R1.bam
    samtools view -bS mapped_reads_R2.sam > mapped_reads_R2.bam
    samtools sort mapped_reads_R1.bam -o sorted_mapped_reads_R1.bam
    samtools sort mapped_reads_R2.bam -o sorted_mapped_reads_R2.bam
    samtools index sorted_mapped_reads_R1.bam
    samtools index sorted_mapped_reads_R2.bam
    samtools view -b -F 4 sorted_mapped_reads_R1.bam > mapped_reads_only_R1.bam
    samtools view -b -F 4 sorted_mapped_reads_R2.bam > mapped_reads_only_R2.bam
    samtools index mapped_reads_only_R1.bam
    samtools index mapped_reads_only_R2.bam
    samtools fastq -F 4 mapped_reads_only_R1.bam > R1.fastq
    samtools fastq -F 4 mapped_reads_only_R2.bam > R2.fastq
    cat R1.fastq R2.fastq > merged.fastq
  """

}

process assembly {

  input:
  path mapped_read

  output:
  path "${params.out_fname}_assembly/scaffolds.fasta", emit: scaffolds

  """
  spades.py -s ${mapped_read} -o ${params.out_fname}_assembly
  """

}

process MapScaffoldsToRef {

  input:
  path Ref
  path scaffolds

  output:
  path "mapped_scaffolds.bam", emit: map_bam
  path "mapped_scaffolds.bam.bai", emit: map_bam_index

  script:
  """
  bowtie2-build $Ref ref_index
  bowtie2 -x ref_index -f ${scaffolds} -S scaffolds_to_reference.sam
  samtools view -S -b scaffolds_to_reference.sam > scaffolds_to_reference.bam
  samtools sort scaffolds_to_reference.bam -o mapped_scaffolds.bam
  samtools index mapped_scaffolds.bam
  """

}

process CleanAndAnnotateScaffolds {

  publishDir "$projectDir", mode: "copy"

  input:
  path Ref
  path map_bam
  path map_bam_index

  output:
  path "${params.out_fname}_recovered_genomes.fasta", emit: final_seq

  script:
  """
  python3 $projectDir/clean_and_annotate_scaffolds.py $Ref $map_bam ${params.out_fname}_recovered_genomes.fasta
  """

}

process CalculateRecoveryRate {

  publishDir "$projectDir", mode: "copy"

  input:
  path Ref
  path final_seq

  output:
  path "${params.out_fname}_recovery_rate.tsv"

  script:
  """
  python3 $projectDir/calculate_recovery_rate.py $Ref $final_seq ${params.out_fname}_recovery_rate.tsv
  """

}

process FastqPair {

  input:
  path extracted_R1
  path extracted_R2

  output:
  path "${extracted_R1}.paired.fq", emit: paired_R1
  path "${extracted_R2}.paired.fq", emit: paired_R2

  """
  fastq_pair $extracted_R1 $extracted_R2
  """

}

process RelativeAbundance {

  publishDir "$projectDir", mode: "copy"

  input:
  path paired_R1
  path paired_R2
  path Ref

  output:
  path "${params.out_fname}_1_abundances.txt", emit: abundance_table

  """
  python3 $projectDir/group.py $Ref reference_grouping.txt
  echo $Ref > ref.txt
  $projectDir/themisto build -k 31 -i ref.txt --index-prefix R --temp-dir ./ --n-threads 16
  awk 'BEGIN {OFS="\\n"; empty_seq = "N"; empty_qual = "F"} {header=\$0; getline seq; getline plus; getline qual; if (length(seq) > 0) print header, seq, plus, qual; else print header, empty_seq, plus, empty_qual}' $paired_R1 > R1_filtered.fastq
  awk 'BEGIN {OFS="\\n"; empty_seq = "N"; empty_qual = "F"} {header=\$0; getline seq; getline plus; getline qual; if (length(seq) > 0) print header, seq, plus, qual; else print header, empty_seq, plus, empty_qual}' $paired_R2 > R2_filtered.fastq
  $projectDir/themisto pseudoalign -q R1_filtered.fastq -o pseudoalignments_1.aln -i R --temp-dir ./ --rc --n-threads 16 --sort-output-lines --gzip-output
  $projectDir/themisto pseudoalign -q R2_filtered.fastq -o pseudoalignments_2.aln -i R --temp-dir ./ --rc --n-threads 16 --sort-output-lines --gzip-output
  mSWEEP --themisto-1 pseudoalignments_1.aln.gz --themisto-2 pseudoalignments_2.aln.gz -o ${params.out_fname} -i reference_grouping.txt
  """

}

process VisualizeAbundance {

  publishDir "$projectDir", mode: "copy"

  input:
  path abundance_table

  output:
  path "${params.out_fname}_relative_abundance.png"

  script:
  """
  python3 $projectDir/visualize_abundance.py ${abundance_table} ${params.out_fname}_relative_abundance.png
  """

}

workflow {
  
  fw_file_ch = Channel.from(params.R1)
  rev_file_ch = Channel.from(params.R2)
  ref_genome = Channel.from(params.Ref)
  
  qc_ch = QC(fw_file_ch, rev_file_ch)

  extracted_read = read_extraction(qc_ch.qc_R1, qc_ch.qc_R2, ref_genome)

  assemble = assembly(extracted_read.mapped_read)

  mb = MapScaffoldsToRef(ref_genome, assemble.scaffolds)

  sequences = CleanAndAnnotateScaffolds(ref_genome, mb.map_bam, mb.map_bam_index)

  CalculateRecoveryRate(ref_genome, sequences.final_seq)

  paired_reads = FastqPair(extracted_read.extracted_R1, extracted_read.extracted_R2)
  
  abundance = RelativeAbundance(paired_reads.paired_R1, paired_reads.paired_R2, ref_genome)

  VisualizeAbundance(abundance.abundance_table)

}