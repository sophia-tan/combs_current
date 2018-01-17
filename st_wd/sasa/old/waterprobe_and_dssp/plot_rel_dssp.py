import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
import seaborn as sns

perc_aai = pkl.load(open('AAi_database_lookups.pkl','rb'))
perc_aai = perc_aai[1]

df = pkl.load(open('rsa_dssp_database_lookups.pkl','rb'))
print(df)
f, axarr = plt.subplots(4,5,figsize=(12,6.5), sharey=True)
fig_row = 0
fig_col = 0

index = [x/20 for x in range(20)][1:]

for aa in df.columns.values:
    x = df[aa]
    x = np.hstack(x)
    x = x/perc_aai[aa]
    axarr[fig_row, fig_col].plot(index, x)
    axarr[fig_row, fig_col].set_title(aa)

    if fig_col < 4:
        fig_col += 1
    else:
        fig_col = 0
        fig_row += 1

plt.tight_layout()
plt.suptitle('relative dssp sasa: f(aai) bin / f(total aa) bin normalized by f(aai) in db')
plt.subplots_adjust(top=0.85) # so suptitle is compatible with tight layout
plt.show()
###
