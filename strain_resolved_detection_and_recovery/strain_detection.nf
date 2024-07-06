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

process sourmash {

  input:
  path qc_R1
  path qc_R2
  path Ref

  output:
  path "${params.output_fname}_presence.csv", emit: presence

  """
  sourmash sketch dna $qc_R1 -o read1.sig
  sourmash sketch dna $qc_R2 -o read2.sig
  sourmash signature merge read1.sig read2.sig -o merged.sig
  sourmash sketch dna $Ref -o genome.sig --singleton
  sourmash gather merged.sig genome.sig -o ${params.output_fname}_presence.csv
  """

}

process FilterReferences {

  publishDir "$projectDir", mode: "copy"

  input:
  path presence
  path Ref

  output:
  path "${params.out_fname}_filtered_ref.fasta", emit: filtered_ref

  script:
  """
  python3 $projectDir/filter_references.py $presence $Ref ${params.out_fname}_filtered_ref.fasta
  """
}

process read_extraction {

  input:
  path qc_R1
  path qc_R2
  path filtered_ref

  output:
  path "merged.fastq", emit: mapped_read
  path "R1.fastq", emit: extracted_R1
  path "R2.fastq", emit: extracted_R2

  """
    minimap2 -d reference_genome.mmi $filtered_ref
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
  path filtered_ref
  path scaffolds

  output:
  path "mapped_scaffolds.bam", emit: map_bam
  path "mapped_scaffolds.bam.bai", emit: map_bam_index

  script:
  """
  bowtie2-build $filtered_ref ref_index
  bowtie2 -x ref_index -f ${scaffolds} -S scaffolds_to_reference.sam
  samtools view -S -b scaffolds_to_reference.sam > scaffolds_to_reference.bam
  samtools sort scaffolds_to_reference.bam -o mapped_scaffolds.bam
  samtools index mapped_scaffolds.bam
  """

}

process CleanAndAnnotateScaffolds {

  publishDir "$projectDir", mode: "copy"

  input:
  path filtered_ref
  path map_bam
  path map_bam_index

  output:
  path "${params.out_fname}_recovered_genomes.fasta", emit: final_seq

  script:
  """
  python3 $projectDir/clean_and_annotate_scaffolds.py $filtered_ref $map_bam ${params.out_fname}_recovered_genomes.fasta
  """

}

process CalculateRecoveryRate {

  publishDir "$projectDir", mode: "copy"

  input:
  path filtered_ref
  path final_seq

  output:
  path "${params.out_fname}_recovery_rate.tsv"

  script:
  """
  python3 $projectDir/calculate_recovery_rate.py $filtered_ref $final_seq ${params.out_fname}_recovery_rate.tsv
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
  path filtered_ref

  output:
  path "${params.out_fname}_1_abundances.txt", emit: abundance_table

  """
  python3 $projectDir/group.py $filtered_ref reference_grouping.txt
  echo $filtered_ref > ref.txt
  $projectDir/themisto build -k 31 -i ref.txt --index-prefix R --temp-dir $projectDir/themisto_index/tmp --mem-gigas 2 --n-threads 4 --file-colors
  awk 'BEGIN {OFS="\n"; empty_seq = "N"; empty_qual = "F"} {header=$0; getline seq; getline plus; getline qual; if (length(seq) > 0) print header, seq, plus, qual; else print header, empty_seq, plus, empty_qual}' $paired_R1 > R1_filtered.fastq
  awk 'BEGIN {OFS="\n"; empty_seq = "N"; empty_qual = "F"} {header=$0; getline seq; getline plus; getline qual; if (length(seq) > 0) print header, seq, plus, qual; else print header, empty_seq, plus, empty_qual}' $paired_R2 > R2_filtered.fastq
  $projectDir/themisto pseudoalign -q R1_filtered.fastq -o pseudoalignments_1.aln -i R --temp-dir $projectDir/themisto_index/tmp --rc --n-threads 16 --sort-output-lines --gzip-output
  $projectDir/themisto pseudoalign -q R2_filtered.fastq -o pseudoalignments_2.aln -i R --temp-dir $projectDir/themisto_index/tmp --rc --n-threads 16 --sort-output-lines --gzip-output
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

  pre = sourmash(qc_ch.qc_R1, qc_ch.qc_R2, ref_genome)

  filtered_refs = FilterReferences(pre.presence, ref_genome)

  extracted_read = read_extraction(qc_ch.qc_R1, qc_ch.qc_R2, filtered_refs.filtered_ref)

  assemble = assembly(extracted_read.mapped_read)

  mb = MapScaffoldsToRef(filtered_refs.filtered_ref, assemble.scaffolds)

  sequences = CleanAndAnnotateScaffolds(filtered_refs.filtered_ref, mb.map_bam, mb.map_bam_index)

  CalculateRecoveryRate(filtered_refs.filtered_ref, sequences.final_seq)

  paired_reads = FastqPair(extracted_read.extracted_R1, extracted_read.extracted_R2)
  
  abundance = RelativeAbundance(paired_reads.paired_R1, paired_reads.paired_R2, filtered_refs.filtered_ref)

  VisualizeAbundance(abundance.abundance_table)

}
