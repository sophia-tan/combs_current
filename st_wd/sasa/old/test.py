import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis
import pandas as pd
import pickle as pkl
from sys import argv

script, csv_dir, ifg = argv

an = analysis.Analyze(csv_dir)
df = an.get_distant_vdms(seq_distance=0)
#df = an.get_distant_vdms(seq_distance=10)
df = df[['iFG_count', 'pdb', 'resname_ifg', 'resnum_ifg', 'chid_ifg', 'vdM_count', 'resname_vdm', 'resnum_vdm', 'chid_vdm', 'sequence_vdm', 'sec_struct_dssp_vdm']]
df = df[:2000]
df = analysis.Analysis.remove_repeat_proteins(df)


pkl.dump(df, open('all_vdm_%s_sasa.pkl' %ifg, 'wb'))

## add sasa info
#sasadf = pd.read_csv('%s_vdm_sasa_info.csv' % ifg)
#merged = pd.merge(df,sasadf, on=['iFG_count', 'vdM_count'])


