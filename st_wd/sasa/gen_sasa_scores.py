# commandline argument is fs, dssp, 3A, 4A, or 5A
import numpy as np
import pickle as pkl
from pprint import pprint
import pandas as pd
from sys import argv

script, ifg = argv

ifg_sasa_pkl = 'noskip_sasa_pkl/noskip_%s_sasa.pkl'%ifg # in wd 
ifg_sasa=pkl.load(open(ifg_sasa_pkl,'rb'))

sasa_dict = {}
AAs = ['ALA', 'ARG', 'VAL', 'MET', 'GLU', 'GLN', 'GLY', \
        'ASP', 'ASN', 'ILE', 'LEU', 'LYS', 'HIS', 'TRP', 'TYR', 'PHE', 'PRO', 'THR', 'SER', 'CYS']
for aa in AAs:
    sasa_dict[aa] = []

# fill sasa_dict with df values
for ix, row in ifg_sasa.iterrows():
    if row['resname_vdm'] in sasa_dict.keys():
        resname = row['resname_vdm']
        sasa_dict[resname].append(row['vdM_sasa_CB_3A_probe'])

# col = aa, row = sasa cutoffs
index = [x*10 for x in range(20)][1:]

all_sasa = [] # for every aa
for aa, x in sasa_dict.items():
    x = [float(z) for z in x ]
    for i in x:
        all_sasa.append(i)

lookup = pkl.load(open('../Lookups/AAi_freq/AAi_database_lookups.pkl','rb')) # lookup table for freq(AAi) in db!

df = pd.DataFrame(index = index, columns = AAs)
for _bin in index:
    bin_sasa = [z<_bin and z>(_bin-20) for z in all_sasa]
    total_aa_bin = sum(bin_sasa) # counts of all aa's in this bin
    for aa, x in sasa_dict.items():
        x = [float(z) for z in x ]
        x = [z<_bin and z>(_bin-20) for z in x]
        aai_counts = sum(x) # at whatever cutoff bin
        frac_aai_bin = aai_counts / total_aa_bin # f(AAi in binj / all aa in binj)
        freq_aai_db = lookup[1][aa]
        score = frac_aai_bin / freq_aai_db
        df.ix[_bin, aa] = round(score,2)
        print(aa, _bin, score)

pkl.dump(df, open('../Lookups/sasa_scores/scores_largeprobe_sasa_%s_lookup.pkl' % ifg,'wb'))
