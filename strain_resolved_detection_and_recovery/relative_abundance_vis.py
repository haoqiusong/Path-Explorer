import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

#np.random.seed(0)
#data = np.random.rand(8, 14)
#data_normalized = data / data.sum(axis=0)

data_normalized = [[0.678697,0.213182,0.0854189,0.0227022,1.90229e-08,1.90229e-08,1.93185e-08,1.90229e-08],
[0.689384,0.236786,0.0614891,0.0123415,1.93944e-08,1.93944e-08,1.94733e-08,1.95521e-08],
[0.695296,0.20478,0.0778727,0.0220529,2.22541e-08,2.22541e-08,2.23547e-08,2.25562e-08],
[0.718129,0.240287,0.0343619,0.00722322,1.59817e-08,1.59817e-08,1.64285e-08,1.60375e-08],
[0.71765,0.24018,0.0348847,0.0072865,1.40617e-08,1.40617e-08,1.4279e-08,1.4366e-08],
[0.579552,0.160161,0.217168,0.0431191,4.25414e-08,4.25414e-08,4.32752e-08,4.25414e-08],
[0.791439,0.17266,0.0286915,0.00721316,5.46196e-08,5.46196e-08,5.75362e-08,5.6077e-08],
[0.583331,0.236531,0.134338,0.0458009,2.45488e-08,2.45488e-08,2.47873e-08,2.46681e-08],
[0.596753,0.260155,0.105548,0.0375474,2.91511e-08,2.91511e-08,2.94927e-08,2.91511e-08],
[0.617477,0.188447,0.138499,0.0555783,2.52124e-08,2.52124e-08,2.5972e-08,2.53389e-08],
[0.644113,0.233177,0.0858068,0.036905,3.06265e-08,3.06265e-08,3.12015e-08,3.06265e-08],
[0.638976,0.23798,0.0856756,0.0373701,2.8725e-08,2.8725e-08,2.97333e-08,2.8893e-08],
[0.531634,0.0837257,0.328465,0.0561748,3.26723e-08,3.26723e-08,3.35231e-08,3.26723e-08],
[0.79623,0.116533,0.0607777,0.0264587,3.61559e-08,3.61559e-08,3.81135e-08,3.61559e-08]]

# List transpose.
data_normalized = list(zip(*data_normalized))

data = pd.DataFrame(data_normalized,
                    columns=['S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12','S13','S14',],
                    index=['Escherichia_coli_ATCC_700973', 'Escherichia_coli_ATCC_11775', 'Enterococcus_faecium_ATCC_BAA_2320', 'Enterococcus_faecalis_ATCC_700802', 'Candida_auris_AR-0390', 'Candida_auris_AR-0385', 'blaOXA-1', 'blaCTX-M'])

# Create the heatmap
sns.set(font_scale=2)
plt.figure(figsize=(25, 15))  # Adjust the size as needed
sns.heatmap(data, annot=False, fmt=".1f", linewidths=.5, cmap='coolwarm', 
            cbar_kws={'label': 'Relative abundance of target genomes'})
#sns.set(font_scale=10)
# Add labels and a title if desired
plt.title('Relative Abundance of Target Genomes in Each Sample')
plt.xlabel('Sample')
plt.ylabel('Target Genome')
plt.tight_layout()
plt.savefig('z.png')