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
#vdms_potentially_coop = sum(vdms['len'])# vdms with at least 1 other vdm in interaction sphere
total_list = [v.split(',') for v in vdms['to_str']]
total_list = [item for sublist in total_list for item in sublist]

aas = sorted(list(one_to_three.keys()))
# calculate how many of each AA makes more than one interaction
new_dict = {}
for amino in aas: 
    amino = one_to_three[amino]
    new_dict[amino] = total_list.count(amino)

# freq of that aa represented in vdm
freq = pkl.load(open(db_aa_lookup,'rb'))['bb_sc freq']

heatmap = pkl.load(open(heatmap_pkl,'rb'))
# generate mask for lower triangle
mask = np.zeros_like(heatmap, dtype=np.bool)
mask[np.triu_indices_from(mask,k=1)] = True
mx = np.ma.masked_array(heatmap, mask=mask)
print(np.sum(mx))
for ix1, aa1 in enumerate(aas): # aa1 = aai
    for ix2, aa2 in enumerate(aas):
        
        # new method i don't believe
        faa1, faa2 = new_dict[one_to_three[aa1]], new_dict[one_to_three[aa2]]
        if faa1 > faa2: 
            denominator = faa2*freq[one_to_three[aa1]]*2
        else:
            denominator = faa1*freq[one_to_three[aa2]]*2



        #denominator = freq[one_to_three[aa2]] * freq[one_to_three[aa1]] * len(df)
        #denominator = freq[one_to_three[aa2]] * freq[one_to_three[aa1]] * vdms_potentially_coop
        score = heatmap[ix1][ix2] / denominator
        print(aa1, aa2, new_dict[one_to_three[aa1]], new_dict[one_to_three[aa2]], freq[one_to_three[aa1]],freq[one_to_three[aa2]], heatmap[ix1][ix2], denominator)
        if score == 0:
            score = np.log10(0.00000000001)
        else:
            score = np.log10(score)
        heatmap[ix1][ix2] = np.round(score,1)

mx = np.ma.masked_array(heatmap, mask=mask)
vmax,vmin,cen = mx.max(), mx.min(), mx.mean()
#print(mx)
#print(vmax,vmin,cen)

sns.heatmap(heatmap, annot=True,mask=mask, cmap="YlGnBu", vmax=vmax,center=cen,
        vmin=vmin,    square=True, linewidths=.5, cbar_kws={"shrink": .5})
ticks = np.arange(20) + 0.5
plt.yticks(ticks,aas[::-1],rotation=0) 
plt.xticks(ticks,aas,rotation=0) 
plt.title('cooperativity score for %s vdms'%ifg)
plt.show()
