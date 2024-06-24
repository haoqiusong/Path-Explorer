"""
import pysam
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def clean_and_annotate_scaffolds(reference, bam_file, output_file):
    # Parse reference lengths
    ref_lengths = {}
    with pysam.FastxFile(reference) as ref:
        for entry in ref:
            ref_lengths[entry.name] = len(entry.sequence)

    # Initialize data structures
    scaffold_seqs = defaultdict(list)
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Extract aligned scaffolds and their positions
    for read in bam.fetch():
        if not read.is_unmapped:
            scaffold_seqs[read.reference_name].append((read.reference_start, read.query_sequence, read.query_qualities))
    
    bam.close()
    
    # Remove overlaps and merge sequences
    cleaned_scaffolds = {}
    for ref, sequences in scaffold_seqs.items():
        sequences.sort()  # Sort by start position
        merged_seq = sequences[0][1]
        last_end = sequences[0][0] + len(sequences[0][1])
        
        for start, seq, qual in sequences[1:]:
            if start < last_end:  # Overlapping sequence
                overlap_len = last_end - start
                non_overlapping_seq = seq[overlap_len:]
                non_overlapping_qual = qual[overlap_len:] if qual else None

                # Compare quality scores for overlapping region
                if qual:
                    overlapping_qual_merged = merged_seq[-overlap_len:]
                    overlapping_qual_new = seq[:overlap_len]

                    best_overlap_seq = ''.join(
                        merged_seq[-overlap_len:][i] if overlapping_qual_merged[i] >= overlapping_qual_new[i] else seq[i]
                        for i in range(overlap_len)
                    )

                    merged_seq = merged_seq[:-overlap_len] + best_overlap_seq + non_overlapping_seq
                else:
                    merged_seq = merged_seq[:-overlap_len] + seq

            else:
                gap_size = start - last_end
                merged_seq += "N" * gap_size + seq
            last_end = start + len(seq)
        
        cleaned_scaffolds[ref] = merged_seq
    
    # Write cleaned and annotated scaffolds to output
    with open(output_file, "w") as out:
        for ref, seq in cleaned_scaffolds.items():
            record = SeqRecord(Seq(seq), id=ref, description="")
            SeqIO.write(record, out, "fasta")

if __name__ == "__main__":
    reference = sys.argv[1]
    bam_file = sys.argv[2]
    output_file = sys.argv[3]
    clean_and_annotate_scaffolds(reference, bam_file, output_file)
"""


"""
import pysam
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def clean_and_annotate_scaffolds(reference, bam_file, output_file):
    # Parse reference sequences
    ref_sequences = {}
    with pysam.FastxFile(reference) as ref:
        for entry in ref:
            ref_sequences[entry.name] = list(entry.sequence)

    # Initialize data structures
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Map reads to the reference sequences
    for read in bam.fetch():
        if not read.is_unmapped:
            for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                if ref_pos is not None:
                    base = read.query_sequence[query_pos]
                    ref_name = read.reference_name
                    ref_sequences[ref_name][ref_pos] = base
    
    bam.close()
    
    # Convert lists back to strings for final output
    for ref in ref_sequences:
        ref_sequences[ref] = "".join(ref_sequences[ref])
    
    # Write cleaned and annotated scaffolds to output
    with open(output_file, "w") as out:
        for ref, seq in ref_sequences.items():
            record = SeqRecord(Seq(seq), id=ref, description="")
            SeqIO.write(record, out, "fasta")

if __name__ == "__main__":
    reference = sys.argv[1]
    bam_file = sys.argv[2]
    output_file = sys.argv[3]
    clean_and_annotate_scaffolds(reference, bam_file, output_file)
"""

"""
import pysam
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def clean_and_annotate_scaffolds(reference, bam_file, output_file):
    # Parse reference sequences and initialize with 'N'
    ref_sequences = {}
    with pysam.FastxFile(reference) as ref:
        for entry in ref:
            ref_sequences[entry.name] = ['N'] * len(entry.sequence)

    # Initialize data structures for storing mapped reads
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Map reads to the reference sequences and store bases and qualities
    scaffold_seqs = defaultdict(lambda: defaultdict(list))
    for read in bam.fetch():
        if not read.is_unmapped and read.query_sequence is not None and read.query_qualities is not None:
            for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                if ref_pos is not None:
                    base = read.query_sequence[query_pos]
                    qual = read.query_qualities[query_pos]
                    scaffold_seqs[read.reference_name][ref_pos].append((base, qual))

    bam.close()

    # Fill the reference sequences with the best base at each position
    for ref_name, positions in scaffold_seqs.items():
        for pos, bases in positions.items():
            if bases:
                base_counter = defaultdict(int)
                for base, qual in bases:
                    base_counter[base] += 1
                most_common_base = max(base_counter, key=base_counter.get)
                ref_sequences[ref_name][pos] = most_common_base

    # Convert lists back to strings for final output
    for ref_name in ref_sequences:
        ref_sequences[ref_name] = ''.join(ref_sequences[ref_name])

    # Write cleaned and annotated scaffolds to output
    with open(output_file, "w") as out:
        for ref_name, seq in ref_sequences.items():
            record = SeqRecord(Seq(seq), id=ref_name, description="")
            SeqIO.write(record, out, "fasta")

if __name__ == "__main__":
    reference = sys.argv[1]
    bam_file = sys.argv[2]
    output_file = sys.argv[3]
    clean_and_annotate_scaffolds(reference, bam_file, output_file)
"""

import sys
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def clean_and_annotate_scaffolds(reference, bam_file, output_file):
    # Load the reference genome
    reference_genome = SeqIO.to_dict(SeqIO.parse(reference, "fasta"))

    # Open the sorted BAM file
    bamfile = pysam.AlignmentFile(bam_file, "rb")

    # Prepare a dictionary to store the recovered sequences
    recovered_sequences = {contig: ["N"] * len(seq) for contig, seq in reference_genome.items()}

    # Iterate over each alignment
    for read in bamfile.fetch():
        if not read.is_unmapped:
            contig = read.reference_name
            start = read.reference_start
            query_sequence = read.query_sequence

            # Place the aligned sequence into the recovered_sequences dictionary
            recovered_sequences[contig][start:start+len(query_sequence)] = list(query_sequence)

    # Convert the recovered sequences to SeqRecord format
    recovered_seq_records = []
    for contig, seq_list in recovered_sequences.items():
        recovered_seq = "".join(seq_list)
        recovered_seq_record = SeqRecord(Seq(recovered_seq), id=contig, description="")
        recovered_seq_records.append(recovered_seq_record)

    # Write the recovered sequences to a FASTA file
    with open(output_file, "w") as output_handle:
        SeqIO.write(recovered_seq_records, output_handle, "fasta")

if __name__ == "__main__":
    reference = sys.argv[1]
    bam_file = sys.argv[2]
    #output_file = '/projects/ciwars/haoqiu_all/CDC_Final/Recovery_Pipeline/recovered_seqs.fasta'
    output_file = sys.argv[3]
    clean_and_annotate_scaffolds(reference, bam_file, output_file)