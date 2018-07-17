''' THIS SCRIPT PRINTS BASH COMMANDS TO RUN CALC_RMSDS
SO I DON'T HAVE TO TYPE IN ALL THE NUMBERS BY HAND '''

import pandas as pd
from PPI_Functions import *

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
    print('nohup python calc_rmsds.py {} planar_group_no_bb >> log'.format(count))
