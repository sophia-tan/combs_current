import pickle as pkl, pandas as pd, seaborn as sns, matplotlib.pyplot as plt
import numpy as np
from AAcodes import *
import itertools
from sys import argv
import scipy

script, heatmap_pkl, vdm_pdb_info_csv, db_aa_lookup = argv

iFG = heatmap_pkl.split('_')[0]
df = pd.read_csv(vdm_pdb_info_csv)
vdms = df[['iFG_count', 'resname']]
def to_str(x):
    return ','.join(sorted(list(x)))
vdms_grouped_ifg = vdms.groupby('iFG_count')
num_ifgs = len(vdms_grouped_ifg)
vdms_grouped_ifg_str = vdms_grouped_ifg['resname'].agg([len, to_str])
vdms = vdms_grouped_ifg_str
vdms = vdms[vdms['len']>1]
num_ifgs_high = len(vdms)

aas = sorted(list(one_to_three.keys()))
heatmap = pd.DataFrame(index=aas, columns=aas)
heatmap = heatmap.replace(np.NaN, 0)
heatmap = heatmap.astype(float)
freq = pkl.load(open(db_aa_lookup,'rb'))[1]

# generate mask for lower triangle
mask = np.zeros_like(heatmap, dtype=np.bool)
mask[np.triu_indices_from(mask,k=1)] = True


###vdms_potentially_coop = sum(vdms['len'])# vdms with at least 1 other vdm in interaction sphere
#### to get total # of possible interactions (pairs), get the length for each vdm (which is the # 
#### of vdms each ifg interacts with), and calculate how many pairs it can make
###total_interactions = [scipy.special.comb(N,2) for N in vdms['len']]
###total_interactions = (sum(total_interactions))
###total_interactions = (scipy.special.comb(vdms_potentially_coop,2))


for ix1, AA1 in enumerate(aas): # aa1 = aai
    for ix2, AA2 in enumerate(aas):
        aa1, aa2 = one_to_three[AA1], one_to_three[AA2]
        count = 0
        for ifg in vdms['to_str']:
            if aa1 in ifg and aa2 in ifg:
                if aa1 == aa2:
                    if ifg.count(aa1) > 1:
                        count += 1
                else: 
                    count += 1
        
        denominator = freq[aa1] * freq[aa2]
        score = count/len(df)/denominator
        #if score == 0:
        #    score = np.log10(0.00001)
        #else:
        score = np.log10(score)
        print(np.round(score,1))
        heatmap[AA1][AA2] = np.round(score,1)
        print(heatmap)
        #print(aa1,aa2, count/len(df), denominator, score)

        mx = np.ma.masked_array(heatmap, mask=mask)
vmax,vmin,cen = mx.max(), mx.min(), mx.mean()
#print(vmax,vmin,cen)

sns.heatmap(heatmap, annot=True,mask=mask)
#sns.heatmap(heatmap, annot=True,mask=mask, cmap="YlGnBu", center=0,vmin=vmin,
#        square=True, linewidths=.5, cbar_kws={"shrink": .5})
#sns.heatmap(heatmap, annot=True,mask=mask, cmap="YlGnBu", vmax=vmax,center=cen,
#        vmin=vmin,    square=True, linewidths=.5, cbar_kws={"shrink": .5})
ticks = np.arange(20) + 0.5
plt.yticks(ticks,aas[::],rotation=0) 
plt.xticks(ticks,aas,rotation=0) 
plt.title('cooperativity score for %s vdms'%iFG)
plt.show()
