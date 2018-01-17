import itertools
import numpy as np

resname_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                   'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                   'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                   'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
                   }

def correlation(df):
    scores = {}
    obs_exp = {}
    df = df[df.apply(get_resnames, axis=1)]
    vdm_gr = df.groupby('iFG_count')
    vdm_gr_agg_pairs = vdm_gr['resname'].agg([num_pairs, to_pairs])
    vdm_gr_agg_pairs_gte1 = vdm_gr_agg_pairs[vdm_gr_agg_pairs['num_pairs'] > 0]
    n_pairs_total = vdm_gr_agg_pairs_gte1['num_pairs'].sum()
    resns = set(df.groupby('resname').groups)
    print(n_pairs_total,'n pairs total')
    
    dataf=vdm_gr_agg_pairs_gte1
    for ix, ifg in dataf.iterrows():
        if ('CYS', 'CYS') in ifg['to_pairs']:
            print(ix)

def num_pairs(x):
    # return len(list(itertools.combinations(sorted(list(x)), 2)))
    # use permutations because
    # using pair[0] and pair[1] in lines 33,36
    return len(list(itertools.combinations(sorted(list(x)), 2)))

def to_pairs(x):
    return list(itertools.combinations(sorted(list(x)), 2))

def get_resnames(row):
    if row['resname'] in set(resname_dict.keys()):
        return True
    else:
        return False

######## START CODE HERE ###########
from sys import argv 
import pandas as pd, pickle as pkl

script, vdm_pdb_info_csv = argv
df = pd.read_csv(vdm_pdb_info_csv)

ifg = vdm_pdb_info_csv.split('_')[0]
corr = correlation(df)
print(counta)
