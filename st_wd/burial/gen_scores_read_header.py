# to calcualte FG burial, its c-beta is used...fix this! 


import sys                              
import pandas as pd
import pickle as pkl
from sys import argv
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis
import numpy as np

script, ifg = argv

df = pkl.load(open('noskip_pkl/%s_noskip.pkl'%ifg,'rb'))
database = pkl.load(open('database_mindistconvexhull.pkl','rb'))
# get burial dist from database
ifgs_df, vdms_df = analysis.Analysis.add_burial_dist(df, database)

AAs = ['ALA', 'ARG', 'VAL', 'MET', 'GLU', 'GLN', 'GLY', \
        'ASP', 'ASN', 'ILE', 'LEU', 'LYS', 'HIS', 'TRP', 'TYR', 'PHE', 'PRO', 'THR', 'SER', 'CYS']
db_aai_freq = pkl.load(open('../Lookups/AAi_freq/AAi_database_lookups.pkl','rb')) 
# col = aa, row = burial dist cutoffs
index = [x for x in range(11)][1:]
ifg_scores_df = pd.DataFrame(index = index, columns = ['score'])
vdm_scores_df = pd.DataFrame(index = index, columns = AAs)

### FOR IFGS
total_ifgs= sum([z<10 for z in ifgs_df['burialdist']]) # don't count info for > 10A bin
for _bin in index: #_bin is the upperbound (doesn't include)
    bin_counts = sum([z<_bin and z>=(_bin-1) for z in ifgs_df['burialdist']])
    frac_bin = bin_counts / total_ifgs
    score = -np.log10(frac_bin / 0.1) # bc 10 bins and 1/10 chance of randomly being in bin
    ifg_scores_df.ix[_bin, 'score'] = round(score,2)

### FOR VDMS
for _bin in index: #_bin is the upperbound (doesn't include)
    total_bin = sum([z<_bin and z>=(_bin-1) for z in vdms_df['burialdist']])
    if _bin==10:
        total_bin = sum([z>=(_bin-1) for z in vdms_df['burialdist']])
    for aa in AAs:
        aai = vdms_df[vdms_df['resname_vdm']==aa]['burialdist']
        aai_bin_counts = sum([z<_bin and z>=(_bin-1) for z in aai])
        if _bin==10:
            aai_bin_counts = sum([z>=(_bin-1) for z in aai])
        frac_aai_bin = aai_bin_counts / total_bin # f(AAi in binj / all aa in binj)
        freq_aai_db = db_aai_freq[1][aa]
        if frac_aai_bin == 0:
            frac_aai_bin = 0.001
        score = -np.log10(frac_aai_bin / freq_aai_db)
        vdm_scores_df.ix[_bin, aa] = round(score,2)
pkl.dump(ifg_scores_df, open('/home/gpu/Sophia/combs/st_wd/Lookups/burial_scores/%s_ifg_burial_scores.pkl'%ifg,'wb'))
pkl.dump(vdm_scores_df, open('/home/gpu/Sophia/combs/st_wd/Lookups/burial_scores/%s_vdm_burial_scores.pkl'%ifg,'wb'))
