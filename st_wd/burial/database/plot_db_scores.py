import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
import seaborn as sns
from sys import argv

ifg_df = pkl.load(open('/home/gpu/Sophia/combs/st_wd/Lookups/burial_scores/db_burial_scores.pkl','rb'))

f, axarr = plt.subplots(4,5,figsize=(12,6.5), sharey=True, sharex=True)
fig_row = 0
fig_col = 0

index = [x for x in range(33)][1:]

for aa in ifg_df.columns.values:
    #ifg = ifg_df[aa] 
    #db = db_df[aa] 
    #scores = ifg/db
    #axarr[fig_row, fig_col].plot(index, scores, lw=3)
    axarr[fig_row, fig_col].plot(index, ifg_df[aa],'r--')
    #axarr[fig_row, fig_col].plot(index, db, 'g--')
    axarr[fig_row, fig_col].set_title(aa)

    if fig_col < 4:
        fig_col += 1
    else:
        fig_col = 0
        fig_row += 1
plt.xlim(1,10)
plt.ylim(0,2.5)
plt.tight_layout()
plt.suptitle('database scores')

plt.subplots_adjust(top=0.85) # so suptitle is compatible with tight layout
plt.show()
###
