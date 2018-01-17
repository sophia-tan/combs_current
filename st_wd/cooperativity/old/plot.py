import pickle as pkl, pandas as pd, seaborn as sns, matplotlib.pyplot as plt
import numpy as np
from AAcodes import *
import itertools
from sys import argv

script, heatmap_pkl, vdm_pdb_info_csv, db_aa_lookup = argv

ifg = heatmap_pkl.split('_')[0]
df = pd.read_csv(vdm_pdb_info_csv)
vdms = df[['iFG_count', 'resname']]
def to_str(x):
    return ','.join(sorted(list(x)))
vdms_grouped_ifg = vdms.groupby('iFG_count')
vdms_grouped_ifg_str = vdms_grouped_ifg['resname'].agg([len, to_str])
vdms_grouped_ifg_str['to_str'].value_counts()
vdms_grouped_ifg_str_counts = vdms_grouped_ifg_str['to_str'].value_counts()
vdms = vdms_grouped_ifg_str
vdms = vdms[vdms['len']>1]
vdms_potentially_coop = sum(vdms['len'])# vdms with at least 1 other vdm in interaction sphere

aas = sorted(list(one_to_three.keys()))

freq = pkl.load(open(db_aa_lookup,'rb'))[0]
freque = pkl.load(open(db_aa_lookup,'rb'))[1]

heatmap = pkl.load(open(heatmap_pkl,'rb'))
# generate mask for lower triangle
mask = np.zeros_like(heatmap, dtype=np.bool)
mask[np.triu_indices_from(mask,k=1)] = True
        
for ix1, aa1 in enumerate(aas): # aa1 = aai
    for ix2, aa2 in enumerate(aas):
        
        # new method i don't believe
        faa1, faa2 = freq[one_to_three[aa1]], freq[on_to_three[aa2]]
        if faa1 > faa2: 
            denominator = faa2*freque[one_to_three[aa1]]
        else:
            denominator = faa1*freque[one_to_three[aa2]]



        #denominator = freq[one_to_three[aa2]] * freq[one_to_three[aa1]] * len(df)
        #denominator = freq[one_to_three[aa2]] * freq[one_to_three[aa1]] * vdms_potentially_coop
        score = heatmap[ix1][ix2] / denominator
        print(score)
        #print(aa1, aa2, heatmap[ix1][ix2], denominator)
        if score == 0:
            score = np.log10(0.00000000001)
        else:
            score = np.log10(score)
        heatmap[ix1][ix2] = np.round(score,1)

mx = np.ma.masked_array(heatmap, mask=mask)
print(mx)
vmax,vmin,cen = mx.max(), mx.min(), mx.mean()
#print(vmax,vmin,cen)

sns.heatmap(heatmap, annot=True,mask=mask, cmap="YlGnBu", vmax=vmax,center=cen,
        vmin=vmin,    square=True, linewidths=.5, cbar_kws={"shrink": .5})
ticks = np.arange(20) + 0.5
plt.yticks(ticks,aas[::-1],rotation=0) 
plt.xticks(ticks,aas,rotation=0) 
plt.title('cooperativity score for %s vdms'%ifg)
plt.show()
