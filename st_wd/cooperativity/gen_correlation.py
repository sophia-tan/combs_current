import numpy as np, pickle as pkl, pandas as pd
import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis, apps
from sys import argv 

script, ifg = argv

csv_dir = '/home/gpu/Sophia/STcombs/20171118/%s/csv'%ifg
# use noskip vdm's from sasa dir, that already have repeats removed, 
# and then get the vdms within 3.5A
combined = analysis.refine_df(csv_dir, seq_dist=0, threefive=False,ifg=ifg)
#dist_vdms, vdms_bb, vdms_sc, an = analysis.refine_df(csv_dir, seq_dist=0, threefive=True,ifg=ifg)

threefive_df = combined[['iFG_count', 'pdb', 'resname_ifg', 'resnum_ifg', 'chid_ifg', 'vdM_count', 'resname_vdm', 'resnum_vdm', 'chid_vdm']]
#threefive_df = threefive_df[['iFG_count', 'pdb', 'resname_ifg', 'resnum_ifg', 'chid_ifg', 'vdM_count', 'resname_vdm', 'resnum_vdm', 'chid_vdm', 'atom_names', 'dist_info']] # atom_names are atoms on vdm

AAifreq = pkl.load(open('/home/gpu/Sophia/combs/st_wd/Lookups/AAi_freq/AAi_freq_combed_{}.pkl'.format(ifg),'rb'))


corr = analysis.EnergyTerms.correlation(threefive_df,AAifreq)
scores,obs_exp = corr

# make heatmap
aas = sorted(list(apps.one_to_three.keys()))
heatmap = pd.DataFrame(index=aas, columns=aas)
heatmap = heatmap.replace(np.NaN, 0)
heatmap = heatmap.astype(float)

for ix1, AA1 in enumerate(aas): # aa1 = aai
    for ix2, AA2 in enumerate(aas):
        aa1, aa2 = apps.one_to_three[AA1], apps.one_to_three[AA2]
        try:
            heatmap[AA1][AA2] = np.round(scores[(aa1,aa2)],2)
        except:
            heatmap[AA1][AA2] = np.round(scores[(aa2,aa1)],2)
print(heatmap)
# format of corr is [scores, obs_exp]
pkl.dump([corr, heatmap], open('../Lookups/correlation/matan_noskip_%s_correlation.pkl'%ifg,'wb'))
