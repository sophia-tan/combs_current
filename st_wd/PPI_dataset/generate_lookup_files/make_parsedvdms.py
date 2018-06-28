import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis
import pandas as pd, prody as pr
import pickle as pkl
from sys import argv
from combs.apps import *

script, ifg_res = argv

csv_dir = '/home/gpu/Sophia/combs/st_wd/20180207_combed_csvs/%s/'%ifg_res
db_pdbs = '/home/gpu/Sophia/combs/st_wd/20180207_db_molprobity_biolassem/'

an = analysis.Analyze(csv_dir)
df = an.get_distant_vdms(seq_distance=10)
df = df[df['resname_vdm'] != 'MSE'] # remove MSE, etc.

df = analysis.Analysis.remove_repeat_proteins(df)
df = df[['pdb', 'chid_ifg', 'resname_ifg', 'resnum_ifg', \
    'chid_vdm', 'resname_vdm', 'resnum_vdm', 'segi_ifg', 'segi_vdm']]

pkl.dump(df, open('/home/gpu/Sophia/combs/st_wd/Lookups/refinedvdms/vdms_of_{}.pkl'.format(ifg_res),'wb'))
