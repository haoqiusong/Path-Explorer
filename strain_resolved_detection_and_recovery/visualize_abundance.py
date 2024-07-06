import pandas as pd
import matplotlib.pyplot as plt
import sys

def visualize_abundance(abundance_file, output_image):
    try:
        df = pd.read_csv(abundance_file, sep='\t', comment='#', header=None, names=['c_id', 'mean_theta'])
    except Exception as e:
        print(f"Error reading the file: {e}")
        return
    
    # Ensure the dataframe is not empty and has the correct columns
    if df.empty or 'c_id' not in df.columns or 'mean_theta' not in df.columns:
        print("The input file does not have the required format.")
        return

    # Plotting
    plt.figure(figsize=(14, 8))
    bars = plt.bar(df['c_id'], df['mean_theta'], color='skyblue')
    plt.xlabel('Genome Name', fontsize=14)
    plt.ylabel('Relative Abundance', fontsize=14)
    plt.title('Relative Abundance of Genomes in the Sample', fontsize=16)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)

    # Adding data labels
    for bar in bars:
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f'{bar.get_height():.4f}', 
                 ha='center', va='bottom', fontsize=12)

    # Adding grid lines
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()
    plt.savefig(output_image)
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python visualize_abundance.py <abundance_file> <output_image>")
        sys.exit(1)

    abundance_file = sys.argv[1]
    output_image = sys.argv[2]
    
    visualize_abundance(abundance_file, output_image)
