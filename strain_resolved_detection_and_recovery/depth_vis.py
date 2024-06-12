import numpy as np
import matplotlib.pyplot as plt
from tsmoothie.smoother import *

f = open("8509-S10_S10_coverage.txt", encoding = "utf-8")
s = f.readlines()
x = []
y = []
for i in range(len(s)):
    s[i] = s[i].strip()
    temp = s[i].split('\t')
    if temp[0].split('_')[1] == '1':
        x.append(int(temp[1]))
        if int(temp[2]) == 0:
            y.append(0)
        else:
            y.append(np.log2(int(temp[2])))

smoother = ConvolutionSmoother(window_len=10000, window_type='ones')
smoother.smooth(y)
#low, up = smoother.get_intervals('sigma_interval', n_sigma=5)

plt.figure(figsize=(10,7))
plt.xlabel("Genome locus (bp)")
plt.ylabel("Depth of coverage (log2-scaled)")
plt.title('Depth of coverage of Escherichia_coli_ATCC_11775 (contig 1) in Sample S10')
plt.plot(smoother.data[0], color='orange', label='original')
plt.plot(smoother.smooth_data[0], linewidth=1.5, color='blue', label='smoothed')
plt.legend()
#plt.fill_between(range(len(smoother.data[0])), low[0], up[0], alpha=0.3)
plt.savefig("Coverage1.png")