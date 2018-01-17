import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis
import pandas as pd
import pickle as pkl
from sys import argv

script, csv_dir, path, ifg = argv

df = pkl.load(open(path, 'rb'))

an = analysis.Analyze(csv_dir)
sasadf = an.vdm_sasa_info
merged = pd.merge(df,sasadf, on=['iFG_count', 'vdM_count'])

pkl.dump(merged, open('sasa_all_vdm_%s_sasa.pkl' %ifg, 'wb'))


