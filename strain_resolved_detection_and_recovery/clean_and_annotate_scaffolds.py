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
        # Remove all 'N's from the recovered sequence
        cleaned_recovered_seq = recovered_seq.replace('N', '')
        recovered_seq_record = SeqRecord(Seq(cleaned_recovered_seq), id=contig, description="")
        recovered_seq_records.append(recovered_seq_record)

    # Write the recovered sequences to a FASTA file
    with open(output_file, "w") as output_handle:
        SeqIO.write(recovered_seq_records, output_handle, "fasta")

if __name__ == "__main__":
    reference = sys.argv[1]
    bam_file = sys.argv[2]
    output_file = sys.argv[3]
    clean_and_annotate_scaffolds(reference, bam_file, output_file)
