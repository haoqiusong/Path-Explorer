import pandas as pd
import matplotlib.pyplot as plt
import sys

def visualize_abundance(abundance_file, output_image):
    # Read the abundance data
    df = pd.read_csv(abundance_file, sep='\t')
    
    # Assuming the TSV file has columns "Strain" and "Relative_Abundance"
    strains = df['Strain']
    abundances = df['Relative_Abundance']
    
    # Create the bar plot
    plt.figure(figsize=(10, 6))
    plt.bar(strains, abundances, color='skyblue')
    
    # Add titles and labels
    plt.xlabel('Strain')
    plt.ylabel('Relative Abundance')
    plt.title('Relative Abundance of Strains')
    
    # Rotate x-axis labels if they are too long
    plt.xticks(rotation=90)
    
    # Save the plot
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