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
    # output each dict as a list. first element = index of 
    # PDB, second element is pdbID, third is resdict
    '''only get destabilizing and ala mutations?'''
    print('dict{} = [{},"{}",{}]'.format(count,count,pdbID,res_dict))


