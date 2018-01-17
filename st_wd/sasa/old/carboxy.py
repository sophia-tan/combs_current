# for merging multiple carboxylate files

import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
import objgraph
from combs import analysis
import pandas as pd
import pickle as pkl
from sys import argv

script = argv
def func():
    ifg = 'backbonCO'
    try:
        an = analysis.Analyze('/home/gpu/Sophia/STcombs/20171118/%s/csv/'%ifg)
    except:
        an = analysis.Analyze('/home/gpu/Sophia/STcombs/20171118/backboneCO/csv/')
    ls = []
    for i in ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v']:
        df = pkl.load(open('noskip_%s_sasa_%s.pkl'%(ifg,i), 'rb'))
        ls.append(df)
    df=pd.DataFrame()
    df = df.append(other=ls, ignore_index=True, verify_integrity=True)
    #df = df[['pdb','resname_ifg', 'resname_vdm', 'sequence_vdm', 'sec_struct_dssp_vdm']]
    df = analysis.Analysis.remove_repeat_proteins(df)
    # add sasa info
    sasadf = an.vdm_sasa_info
    merged = pd.merge(df,sasadf, on=['iFG_count', 'vdM_count'])
    
    pkl.dump(merged, open('noskip_%s_sasa.pkl' %ifg, 'wb'))
func()
