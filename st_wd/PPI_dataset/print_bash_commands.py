''' print bash commands to run calc_rmsds'''

import pandas as pd
from Functions import *

'''load dataset'''
zemu = pd.read_csv('kortemme_flexddg_dataset.csv')

'''group dataset by pdb'''
grouped = zemu.groupby('PDBFileID')
count = 0
for pdbID, rows in grouped:
    count += 1
    res_dict = generate_resdict(rows)
    # output each dict as a list. first element = pdbID, 
    # second element = resdict 
    for key,resname in res_dict.items():
        print('python calc_rmsds.py {} {} {}'.format(count,
        key[0],key[1:]), 'planar_group & > log')

