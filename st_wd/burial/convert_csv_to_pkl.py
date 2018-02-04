import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis
import pandas as pd
import pickle as pkl
from sys import argv

script, ifg = argv

ifgs_dir = '/home/gpu/Sophia/STcombs/20171118/'

if ifg == 'lonepair_imidazole': 
    df = analysis.refine_df(ifgs_dir, 0, False, ifg, partial_ifg='lonepair_imidazole')
else:
    df = analysis.refine_df(ifgs_dir, 0, False, ifg)
df = df[['iFG_count', 'pdb', 'resname_ifg', 'resnum_ifg', 'chid_ifg', 'resname_vdm', 'resnum_vdm', 'chid_vdm', 'vdM_count', 'bb', 'sc']]
pkl.dump(df, open('noskip_pkl/%s_noskip.pkl'%ifg,'wb'))
