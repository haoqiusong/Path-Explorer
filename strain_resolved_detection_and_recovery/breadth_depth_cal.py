import os

def calculate_breadth_of_coverage(coverage_file, coverage_threshold):
    covered_positions = 0
    total_positions = 0

    with open(coverage_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            coverage_depth = int(fields[2])  # Assuming depth is in the third column
            total_positions += 1
            if coverage_depth >= coverage_threshold:
                covered_positions += 1

    breadth_of_coverage = (covered_positions / total_positions) * 100
    return breadth_of_coverage


def calculate_depth_of_coverage(coverage_file, coverage_threshold):
    coverage_values = []
    with open(coverage_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            parts = line.split('\t')  # Split line by tab or space, depending on the file format
            coverage = int(parts[2])  # Assuming coverage value is in the second column
            coverage_values.append(coverage)
    total_coverage = sum(coverage_values)
    total_positions = len(coverage_values)
    average_depth = total_coverage / total_positions
    return average_depth


coverage_threshold = 1  # Set your desired coverage threshold
dir_path = '/projects/ciwars/haoqiu_all/cdc/breadth/scripts'

for file_path in os.listdir(dir_path):
    if os.path.isfile(os.path.join(dir_path, file_path)):
        if 'coverage' in file_path:
            sample_name = file_path.split('_')[1]
            breadth = calculate_breadth_of_coverage(file_path, coverage_threshold)
            depth = calculate_depth_of_coverage(file_path, coverage_threshold)
            print(sample_name)
            print(f"Breadth of Coverage: {breadth:.2f}%")
            print(f"Depth of Coverage: {depth:.2f}X\n")