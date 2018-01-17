import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis
import pandas as pd
import pickle as pkl
from sys import argv

script, ifg = argv

csv_dir = '/home/gpu/Sophia/STcombs/20171118/%s/csv'%ifg

an = analysis.Analyze(csv_dir)
df = an.get_distant_vdms(seq_distance=0)

# if imidazole, the only ifg relevant for apixaban is it if has at least one lone pair
# check for this by making sure there's no corresponding H 
if ifg == 'lonepair_imidazole':
    df = an.get_lonepair_imidazole_ifgs(dist_df = df)

df = df[['iFG_count', 'pdb', 'resname_ifg', 'resnum_ifg', 'chid_ifg', 'vdM_count', 'resname_vdm', 'resnum_vdm', 'chid_vdm', 'sequence_vdm', 'sec_struct_dssp_vdm']]
df = analysis.Analysis.remove_repeat_proteins(df)

# add sasa info
sasadf = an.vdm_sasa_info
merged = pd.merge(df,sasadf, on=['iFG_count', 'vdM_count'])

pkl.dump(merged, open('noskip_sasa_pkl/noskip_%s_sasa.pkl' %ifg, 'wb'))
