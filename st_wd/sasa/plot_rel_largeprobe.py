import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
import seaborn as sns
from sys import argv

script, dbscores, ifgscores, ifg = argv
ifg_df = pkl.load(open(ifgscores,'rb'))
db_df = pkl.load(open(dbscores,'rb'))

f, axarr = plt.subplots(4,5,figsize=(12,6.5), sharey=True, sharex=True)
fig_row = 0
fig_col = 0

index = [x*10 for x in range(20)][1:]

for aa in db_df.columns.values:
    #carb = ifg_df[0][aa] # carbdf has 3a and 4a. select 3A
    #db = db_df[2][aa] # db has 1.4, 1.8, 3, and 4. select 3A only
    carb = ifg_df[aa] 
    db = db_df[aa] 
    scores = carb/db
    axarr[fig_row, fig_col].plot(index, scores, lw=3)
    axarr[fig_row, fig_col].plot(index, carb,'r--')
    axarr[fig_row, fig_col].plot(index, db, 'g--')
    axarr[fig_row, fig_col].set_title(aa)

    if fig_col < 4:
        fig_col += 1
    else:
        fig_col = 0
        fig_row += 1

plt.tight_layout()
plt.suptitle('normalized 3A probe sasa at Cb: r=%s, g=db, b=ifg/db' %ifg)

plt.subplots_adjust(top=0.85) # so suptitle is compatible with tight layout
plt.show()
###
