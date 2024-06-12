import gzip
import pysam


# R1
for i in range(1, 15):
    sample = "/projects/ciwars/haoqiu_all/cdc/extraction/8509-S" + str(i) + "_S" + str(i) + "_CAT_mapped_reads_only_R1.bam"
    read_file = "/projects/ciwars/CDC-WBS/TestRun/RawData_LanesConcatenated/8509-S" + str(i) + "_S" + str(i) + "_CAT_R1_001.fastq.gz"

    mapped_read = []
    mapped_read_seq = []

    samfile = pysam.AlignmentFile(sample, "rb")

    references = samfile.references
    lengths = samfile.lengths

    for ref, lens in zip(references, lengths):
        contig_reads = samfile.fetch(ref)   
        for read in contig_reads:
            mapped_read.append(read.query_name)

    with gzip.open(read_file, 'rb') as s:
        while True:
            lines = [next(s) for _ in range(4)]
            if not lines:
                break
            lines[0] = lines[0].decode('ascii').strip()
            lines[1] = lines[1].decode('ascii').strip()
            lines[2] = lines[2].decode('ascii').strip()
            lines[3] = lines[3].decode('ascii').strip()
            tem_name = lines[0].split(' ')[0][1:]
            if tem_name in mapped_read:
                mapped_read_seq += lines
    
    with open('S'+str(i)+'_R1.fastq', 'w') as file:
        for item in mapped_read_seq:
            file.write("%s\n" % item)


# R2
for i in range(1, 15):
    sample = "/projects/ciwars/haoqiu_all/cdc/extraction/8509-S" + str(i) + "_S" + str(i) + "_CAT_mapped_reads_only_R2.bam"
    read_file = "/projects/ciwars/CDC-WBS/TestRun/RawData_LanesConcatenated/8509-S" + str(i) + "_S" + str(i) + "_CAT_R2_001.fastq.gz"

    mapped_read = []
    mapped_read_seq = []

    samfile = pysam.AlignmentFile(sample, "rb")

    references = samfile.references
    lengths = samfile.lengths

    for ref, lens in zip(references, lengths):
        contig_reads = samfile.fetch(ref)   
        for read in contig_reads:
            mapped_read.append(read.query_name)

    with gzip.open(read_file, 'rb') as s:
        while True:
            lines = [next(s) for _ in range(4)]
            if not lines:
                break
            lines[0] = lines[0].decode('ascii').strip()
            lines[1] = lines[1].decode('ascii').strip()
            lines[2] = lines[2].decode('ascii').strip()
            lines[3] = lines[3].decode('ascii').strip()
            tem_name = lines[0].split(' ')[0][1:]
            if tem_name in mapped_read:
                mapped_read_seq += lines
    
    with open('S'+str(i)+'_R2.fastq', 'w') as file:
        for item in mapped_read_seq:
            file.write("%s\n" % item)