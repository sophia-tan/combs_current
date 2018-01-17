# combine vdm sasa old info (3A, 4A, 5A) with 2A
import sys                              
import pandas as pd
import pickle as pkl

allvdms = pd.read_pickle(open('all_vdms_df.pkl','rb')) # contains info like vdm resname, chain etc.
olddf = pd.read_csv('carboxylate_vdm_sasa_info.csv')# old sasa df only has 3A, 4A, 5A

merged = pd.merge(allvdms, olddf, on=['iFG_count', 'vdM_count'])

pkl.dump(merged, open('carboxylate_all_sasa_methods.pkl','wb'))

