# commandline argument is fs, dssp, 3A, 4A, or 5A
import numpy as np
import pickle as pkl
from pprint import pprint
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

nrpdb_dict=pkl.load(open('sasa_database_from_07252017.pkl','rb'))


f, axarr = plt.subplots(4,5,figsize=(12,6.5), sharey=True, sharex=True)
fig_row = 0
fig_col = 0

for aa, resdict in nrpdb_dict.items():
    for m in ['4A', '3A']:
    #for m in ['4A', '3A', '1.8A', '1.4A']:
        x = resdict[m]
        x = [float(z) for z in x ]
        x = np.hstack(x)
        x[x==0] = 1 # bc don't want to take log(0)
        x = np.hstack(x)

        axarr[fig_row, fig_col].hist(x, bins=15)
        axarr[fig_row, fig_col].set_title(aa)

    if fig_col < 4:
        fig_col += 1
    else:
        fig_col = 0
        fig_row += 1

plt.xlim(0,200)
f.text(0.5, 0.95, 'nrPDB sasa cb probe: g=3A, b=4A', ha='center')
f.text(0.5, 0.04, 'sasa sq. angstroms', ha='center')
f.text(0.04, 0.5, 'counts', va='center', rotation='vertical')
#plt.tight_layout()
plt.yscale('log', nonposy='clip')
plt.show()
###
