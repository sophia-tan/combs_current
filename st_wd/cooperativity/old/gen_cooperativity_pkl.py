import pickle as pkl, pandas as pd, seaborn as sns, matplotlib.pyplot as plt
import numpy as np
from AAcodes import *
import itertools
from sys import argv

script, vdm_pdb_info_csv, ifg = argv

aas = sorted(list(one_to_three.keys()))

df = pd.read_csv(vdm_pdb_info_csv)
vdms = df[['iFG_count', 'resname']]

heatmap = pd.DataFrame(index=aas, columns=aas)
heatmap = heatmap.replace(np.NaN, 0)

def to_str(x):
    return ','.join(sorted(list(x)))
vdms_grouped_ifg = vdms.groupby('iFG_count')
vdms_grouped_ifg_str = vdms_grouped_ifg['resname'].agg([len, to_str])
vdms_grouped_ifg_str['to_str'].value_counts()
vdms = vdms_grouped_ifg_str
vdms = vdms[vdms['len']>1]

for ix, row in vdms.iterrows():
    vdmls = row['to_str']
    vdmls = vdmls.split(',')

    pairs = itertools.combinations(vdmls,2)
    for p in pairs:
        if p[0] in three_to_one.keys() and p[1] in three_to_one.keys():
            heatmap.ix[three_to_one[p[0]],three_to_one[p[1]]] += 1
            heatmap.ix[three_to_one[p[1]],three_to_one[p[0]]] += 1

heatmap = heatmap[heatmap.columns].astype(int)
heatmap = np.tril(heatmap,k=1)
heatmap = heatmap.astype(float)

pkl.dump(heatmap, open('%s_cooperativity_counts.pkl'%ifg,'wb'))

