import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
import seaborn as sns
from sys import argv

script, ifg = argv

db_df = pkl.load(open('/home/gpu/Sophia/combs/st_wd/Lookups/burial_scores/db_burial_scores.pkl','rb'))

vdm_df = pkl.load(open('/home/gpu/Sophia/combs/st_wd/Lookups/burial_scores/%s_ifg_burial_scores.pkl'%ifg,'rb'))

f, axarr = plt.subplots(4,5,figsize=(12,6.5), sharey=True, sharex=True)
fig_row = 0
fig_col = 0

index = [x for x in range(11)][1:]

for aa in vdm_df.columns.values:
    axarr[fig_row, fig_col].plot(index, vdm_df[aa],'g')
    #axarr[fig_row, fig_col].plot(index, db_df[[aa]], 'r--')
    axarr[fig_row, fig_col].set_title(aa)

    if fig_col < 4:
        fig_col += 1
    else:
        fig_col = 0
        fig_row += 1
plt.tight_layout()
plt.suptitle('burial scores for vdms of %s; red=db, green=ifg'%ifg)

plt.subplots_adjust(top=0.85) # so suptitle is compatible with tight layout
plt.show()
