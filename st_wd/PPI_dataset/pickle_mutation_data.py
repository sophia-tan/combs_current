# THIS SCRIPT GENERATES RESDICTS

import pandas as pd, pickle as pkl
from Functions import *

'''load dataset'''
zemu = pd.read_csv('kortemme_flexddg_dataset.csv')

'''group dataset by pdb'''
grouped = zemu.groupby('PDBFileID')
resdicts_list = []
datadf_list = []

for pdbID, rows in grouped:
    res_dict = generate_resdict(rows)
    data_df = generate_data_df(rows)
    resdicts_list.append([pdbID, res_dict])
    datadf_list.append(data_df)

pkl.dump([resdicts_list, datadf_list], 
    open('mutation_data.pkl','wb'))
