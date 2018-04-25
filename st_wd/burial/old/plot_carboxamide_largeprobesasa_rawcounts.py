# commandline argument is fs, dssp, 3A, 4A, or 5A
import numpy as np
import pickle as pkl
from pprint import pprint
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sasa_dict = pkl.load(open('carboxamide_largeprobesasa_rawcountsdict.pkl','rb'))


f, axarr = plt.subplots(4,5,figsize=(12,6.5), sharey=True, sharex=True)
fig_row = 0
fig_col = 0

#c = ['blue', 'green', 'purple']
for aa, resdict in sasa_dict.items():
    for index, m in enumerate(['4A', '3A']):
        x = resdict[m]
        x = [float(z) for z in x ]
        x[x==0] = 1
        x = np.hstack(x)
        axarr[fig_row, fig_col].hist(x)
        axarr[fig_row, fig_col].set_title(aa)

    if fig_col < 4:
        fig_col += 1
    else:
        fig_col = 0
        fig_row += 1

plt.xlim(0,200)
f.text(0.5, 0.95, 'carboxamide sasa cb probe: g=3A, b=4A', ha='center')
f.text(0.5, 0.04, 'sasa sq. angstroms', ha='center')
f.text(0.04, 0.5, 'log(counts)', va='center', rotation='vertical')
#plt.tight_layout()
plt.yscale('log', nonposy='clip')
plt.show()
###
