import pickle as pkl, pandas as pd, seaborn as sns, matplotlib.pyplot as plt
import numpy as np
from AAcodes import *
import itertools

sns.set(style="white")
freq = pkl.load(open('AAi_old3.5_combed_carboxamide_lookup.pkl','rb'))
totalcounts = sum(list(freq.values()))
print(totalcounts)
######df = pd.read_csv('carboxamide_vdm_pdb_info.csv')
######vdms = df[['iFG_count', 'resname']]
######
aas = sorted(list(one_to_three.keys()))
######heatmap = pd.DataFrame(index=aas, columns=aas)
######for aa1 in aas:
######    for aa2 in aas:
######        heatmap.ix[aa1,aa2] = 0
######
#######for ifgcount in range(1, 90):
######rowcount = 0
######for ifgcount in range(1, 99743):
######    vdmls = []
######    for ix, row in vdms[max(0,rowcount-10):min(rowcount+10,99843)].iterrows():
######        if row['iFG_count'] == ifgcount:
######            rowcount += 1
######            print(rowcount)
######            try:
######                res = three_to_one[row['resname']]
######                vdmls.append(row['resname'])
######            except:
######                pass
######    pairs = itertools.combinations(vdmls,2)
######    for p in pairs:
######        heatmap.ix[three_to_one[min(p)],three_to_one[max(p)]] += 1
######        heatmap.ix[three_to_one[max(p)],three_to_one[min(p)]] += 1
######heatmap = heatmap[heatmap.columns].astype(int)
######heatmap = np.triu(heatmap,k=1)
######pkl.dump(heatmap,open('heatmap.pkl','wb'))
heatmap = pkl.load(open('heatmap.pkl','rb'))
heatmap = heatmap.astype(float)

for ix1, aa1 in enumerate(aas): # aa1 = aai
    for ix2, aa2 in enumerate(aas):
        denominator = freq[one_to_three[aa2]] * freq[one_to_three[aa1]] / totalcounts / totalcounts * 48013
        #denominator = freq[one_to_three[aa2]] * freq[one_to_three[aa1]] / totalcounts / totalcounts * 190281
        print(heatmap[ix1][ix2], denominator)
        score = heatmap[ix1][ix2] / denominator
        if score == 0:
            score = np.log(1)
        else:
            score = np.log(score)
        heatmap[ix1][ix2] = np.round(score,1)

# generate mask for lower triangle
mask = np.zeros_like(heatmap, dtype=np.bool)
mask[np.tril_indices_from(mask)] = True
cmap = sns.diverging_palette(220, 10, as_cmap=True)

for rowix, row in enumerate(heatmap):
    for colix, col in enumerate(row):
        heatmap[rowix][colix] = heatmap[rowix][colix]

sns.heatmap(heatmap, mask=mask, cmap=cmap, annot=True,
            square=True, linewidths=.5, cbar_kws={"shrink": .5})

ticks = np.arange(20) + 0.5
plt.yticks(ticks,aas[::-1],rotation=0) 
plt.xticks(ticks,aas[::-1],rotation=0) 
plt.show()
