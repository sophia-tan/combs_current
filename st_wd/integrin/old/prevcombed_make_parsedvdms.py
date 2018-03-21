import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis
import pandas as pd, prody as pr
import pickle as pkl
from sys import argv

script, ifg_res = argv

csv_dir = '/home/gpu/Sophia/STcombs/20171118/%s/csv/'%ifg_res
db_pdbs = '/home/gpu/Sophia/STcombs/database/reduce'

an = analysis.Analyze(csv_dir)
df = an.get_distant_vdms(seq_distance=10)
print(len(df), 'dist')

df = analysis.Analysis.remove_repeat_proteins(df)
df = df[['pdb', 'chid_ifg', 'resname_ifg', 'resnum_ifg', \
    'chid_vdm', 'resname_vdm', 'resnum_vdm']]

pkl.dump(df, open('/home/gpu/Sophia/combs/st_wd/Lookups/refinedvdms/vdms_of_{}.pkl'.format(ifg_res),'wb'))
