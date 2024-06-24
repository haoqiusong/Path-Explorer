import pandas as pd
import sys

def filter_references(presence_csv, ref_fasta, output_fasta):
    # Read the presence CSV file
    df = pd.read_csv(presence_csv)
    
    # Get the list of present genome names
    present_genomes_old = df['name'].tolist()

    present_genomes = []
    for pre in present_genomes_old:
        present_genomes.append(pre.split(' ')[0])

    # Read the reference fasta file
    with open(ref_fasta, 'r') as ref_file:
        lines = ref_file.readlines()
    
    # Filter out absent genomes
    write_genome = False
    with open(output_fasta, 'w') as out_file:
        for line in lines:
            if line.startswith('>'):
                genome_name = line.split()[0][1:]
                if genome_name in present_genomes:
                    write_genome = True
                else:
                    write_genome = False
            if write_genome:
                out_file.write(line)

if __name__ == '__main__':
    presence_csv = sys.argv[1]
    ref_fasta = sys.argv[2]
    output_fasta = sys.argv[3]
    
    filter_references(presence_csv, ref_fasta, output_fasta)