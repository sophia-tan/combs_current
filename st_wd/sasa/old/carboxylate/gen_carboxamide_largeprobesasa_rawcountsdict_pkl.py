# commandline argument is fs, dssp, 3A, 4A, or 5A
import numpy as np
import pickle as pkl
from pprint import pprint
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

carbox_sasa=pkl.load(open('carboxylate_all_sasa_methods.pkl','rb')) # df

sasa_dict = {}
AAs = ['ALA', 'ARG', 'VAL', 'MET', 'GLU', 'GLN', 'GLY', \
        'ASP', 'ASN', 'ILE', 'LEU', 'LYS', 'HIS', 'TRP', 'TYR', 'PHE', 'PRO', 'THR', 'SER', 'CYS']
for aa in AAs:
    sasa_dict[aa] = {'1.8A':[], '3A': [], '4A': []}


# fill aa_dict with df values
for ix, row in carbox_sasa.iterrows():
    if row['resname_vdm'] in sasa_dict.keys():
        resname = row['resname_vdm']
        sasa_dict[resname]['3A'].append(row['vdM_sasa_CB_3A_probe'])
        sasa_dict[resname]['4A'].append(row['vdM_sasa_CB_4A_probe'])


pkl.dump(sasa_dict, open('carboxylate_largeprobesasa_rawcountsdict.pkl','wb'))


