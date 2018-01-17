import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
import seaborn as sns

# old method normalized by counts in bin / counts over all bins
#carbox_df = pkl.load(open('rel_largeprobe_sasa_carboxamide.pkl','rb'))
#db_df = pkl.load(open('rel_largeprobe_sasa_database_lookups.pkl','rb'))

carbox_df = pkl.load(open('scores_largeprobe_sasa_carboxylate_lookup.pkl','rb'))
db_df = pkl.load(open('scores_largeprobe_sasa_db_lookup.pkl','rb'))

f, axarr = plt.subplots(4,5,figsize=(12,6.5), sharey=True, sharex=True)
fig_row = 0
fig_col = 0

index = [x*10 for x in range(20)][1:]

for aa in db_df.columns.values:
    #carb = carbox_df[0][aa] # carbdf has 3a and 4a. select 3A
    #db = db_df[2][aa] # db has 1.4, 1.8, 3, and 4. select 3A only
    carb = carbox_df[aa] 
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
plt.suptitle('normalized 3A probe sasa at Cb: r=carboxamide, g=db, b=carboxamide/db')

plt.subplots_adjust(top=0.85) # so suptitle is compatible with tight layout
axes = plt.gca()
axes.set_ylim([0,3])

plt.show()
###
