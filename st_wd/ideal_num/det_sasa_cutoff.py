import numpy as np
import pickle as pkl, pandas as pd, matplotlib.pyplot as plt
from sys import argv

script, ifg, csvdir, threshold = argv
# threshold is percent, cutoff is in angstroms

sasadf = pd.read_csv(csvdir+ifg+'_ifg_atom_density.csv', index_col = 'iFG_count')

hist, bin_edges = np.histogram(sasadf['iFG_sasa_CB_3A_probe'], bins=300)
cutoff = [c/sum(hist) for c in np.cumsum(hist) ]
cutoff = [i for i, e in enumerate(cutoff) if e > float(threshold)/100]
cutoff = cutoff[0] # this is now the index where to draw the line
cutoff = bin_edges[cutoff] # this is where to draw the line

# see what percentage fall under 10 anstroms
ten = len(bin_edges[bin_edges<10]) - 1
print('percentage of data that falls under 10 A**2')
print(np.cumsum(hist)[ten] / sum(hist)*100)


plt.hist(sasadf['iFG_sasa_CB_3A_probe'],bins=100)
plt.title('%s CB sasa, %s percent at %s A^2' %(ifg, threshold, cutoff)) 
plt.show()
