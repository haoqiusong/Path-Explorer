import pysam

def calculate_breadth_of_coverage(bam_file, coverage_threshold):
    total_bases = 0
    covered_bases = 0

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for pileupcolumn in bam.pileup():
            total_bases += 1
            if pileupcolumn.n >= coverage_threshold:
                covered_bases += 1

    breadth_of_coverage = covered_bases / total_bases
    return breadth_of_coverage

# Usage example
bam_file = "/projects/ciwars/haoqiu_all/cdc/extraction/8509-S1_S1_CAT_mapped_reads_only_R1.bam"
coverage_threshold = 10
breadth = calculate_breadth_of_coverage(bam_file, coverage_threshold)
print(f"Breadth of coverage at {coverage_threshold}x: {breadth * 100:.2f}%")