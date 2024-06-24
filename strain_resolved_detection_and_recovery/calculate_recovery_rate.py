import pysam
import sys

def calculate_recovery_rate(reference, fasta_file, output_file):
    # Parse reference lengths
    ref_lengths = {}
    with pysam.FastxFile(reference) as ref:
        for entry in ref:
            ref_lengths[entry.name] = len(entry.sequence)

    # Initialize coverage dictionary
    coverage = {ref: 0 for ref in ref_lengths}

    # Calculate coverage by counting non-'N' bases
    with pysam.FastxFile(fasta_file) as fasta:
        for entry in fasta:
            if entry.name in coverage:
                coverage[entry.name] += sum(1 for base in entry.sequence if base != 'N')

    # Write recovery rates to output file
    with open(output_file, 'w') as out:
        out.write("Reference\tRecoveryRate(%)\n")
        for ref in ref_lengths:
            recovery_rate = (coverage[ref] / ref_lengths[ref]) * 100 if ref_lengths[ref] > 0 else 0
            out.write(f"{ref}\t{recovery_rate:.2f}\n")

if __name__ == "__main__":
    reference = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]
    calculate_recovery_rate(reference, fasta_file, output_file)