import pickle as pkl, pandas as pd, seaborn as sns, matplotlib.pyplot as plt
import numpy as np
from sys import argv

script, ifg = argv
corr, heatmap = pkl.load(open('../Lookups/noskip_%s_correlation.pkl'%ifg, 'rb'))

mask = np.zeros_like(heatmap, dtype=np.bool)
mask[np.triu_indices_from(mask,k=1)] = True

aas = list(heatmap.index)
sns.heatmap(heatmap, annot=True,mask=mask)
ticks = np.arange(20) + 0.5
plt.yticks(ticks,aas[::],rotation=0) 
plt.xticks(ticks,aas,rotation=0) 
plt.title('%s cooperativity score'%ifg)
plt.show()
