# see where the carboxamide trp-trp comes from
import pickle as pkl, pandas as pd
import numpy as np
from AAcodes import *
import itertools
from sys import argv

script, vdm_pdb_info_csv, ifg_pdb_info_csv, aa = argv

aas = sorted(list(one_to_three.keys()))

df = pd.read_csv(vdm_pdb_info_csv)
ifg_df = pd.read_csv(ifg_pdb_info_csv)
ifg_df = ifg_df[['iFG_count', 'pdb', 'resname', 'resnum', 'chid']]

vdms = df[['iFG_count', 'resname']] # resname is vdm, resname_y is ifg
df = df[['iFG_count', 'vdM_count', 'resname', 'resnum', 'chid', 'atom_names']]

def to_str(x):
    return ','.join(sorted(list(x)))
vdms_grouped_ifg = vdms.groupby('iFG_count')

totalcount = 0 # ex) how many trps there are if the vdms of interest are trps
hassccount = 0
bb = ['N', 'O', 'CA', 'C']
for ix, ifggroup in vdms_grouped_ifg: # for all the ifgs 
    if len(np.where(ifggroup['resname']==aa)[0]) > 1:
        iFG_count = list(set(ifggroup['iFG_count']))[0]
        interest = df[df['iFG_count']==iFG_count][df['resname']==aa]
        print(interest)
        print(ifg_df[ifg_df['iFG_count']==iFG_count])
        print('break')
        for v, vdmrow in interest.iterrows(): # for all the vdms of interest 
            totalcount += 1
            atoms = vdmrow['atom_names'].split(' ')
            # get all the bb atoms. this gets you intersection of atom_names list and bb
            inter = list(set(bb) & set(atoms))
            if len(inter) < len(atoms): # that means it has at least one sc atom
                hassccount += 1

#print(totalcount) 
#print(hassccount)
# this means out of all __totalcount__ vdms, __hassccount__ have possible sc interactions



