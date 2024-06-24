import re
import sys
import gzip
from collections import defaultdict

def parse_fasta(fasta_file):
    if fasta_file.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'
    
    with open_func(fasta_file, mode) as f:
        lines = f.readlines()
    
    genome_names = defaultdict()
    for line in lines:
        if line.startswith('>'):
            genome_name = re.search(r'>(\S+)', line).group(1)
            genome_names[genome_name] = line[1:-1]
    
    return genome_names

def create_reference_grouping_file(genome_names, output_file):
    with open(output_file, 'w') as f:
        for k, v in genome_names.items():
            f.write(v + '\t' + k + '\n')

if __name__ == '__main__':
    fasta_file = sys.argv[1]
    output_file = 'reference_grouping.txt'
    
    genome_names = parse_fasta(fasta_file)
    create_reference_grouping_file(genome_names, output_file)