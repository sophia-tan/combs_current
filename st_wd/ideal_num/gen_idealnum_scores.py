# assumes that ifg_interaction_types.pkl is in working dir
import matplotlib.pyplot as plt, pickle as pkl, pandas as pd
from matplotlib import cm
import numpy as np
from sys import argv
from scipy import stats

script, ifg = argv
#sasa is the cutoff in A**2
sasa = 10
csv_dir = '/home/gpu/Sophia/STcombs/20171118/%s/csv/' %ifg
sasacsvpath=csv_dir+'%s_ifg_atom_density.csv'%ifg

sasacsv = pd.read_csv(sasacsvpath,index_col='iFG_count')
contactsdf, megalist = pkl.load(open('interaction_types_pkl/%s_interaction_types.pkl'%ifg,'rb'))

########## ONLY GET ELEMENTS IN LIST THAT ARE BELOW SASA CUTOFF ##################
# get iFGs that are below cutoff 
# poss_ifgs has more ifgs than is in the contacts pkl file bc contacts pkl file 
# has repeats removed
poss_ifgs = sasacsv.index[sasacsv['iFG_sasa_CB_3A_probe'] < float(sasa)].tolist()
# for every ifg in contactsdf (that has repeats removed), if it's also in poss_ifgs, then keep 
# the index number
contactsdf_ifgs = list(set(contactsdf['iFG_count']))
contactsdf_ifgs = np.array(contactsdf_ifgs)


indices = np.isin(contactsdf_ifgs, np.array(poss_ifgs))
indices = list(np.where(indices==True)[0])
num_ifgs_within_sasa_cutoff = len(indices)
newlist = []
for ls in megalist: 
    newlist.append([ls[i] for i in indices])
megalist = newlist
[water_contacts_ifglevel , hb_contacts_vdmlevel , vdw_contacts_vdmlevel,
    polar_contacts_vdmlevel, hb_contacts_atomlevel, ca_hb_contacts_atomlevel, 
    vdw_contacts_atomlevel, polar_contacts_atomlevel, totalhbond ] = megalist 

totalhbonds = np.histogram(totalhbond, bins=int(max(totalhbond))) ## list of counts starts 
# at the 0th # of bonds (so same as the bin edge). the LAST ONE is the sum of the 
# last 2 bin edges
vdws =np.histogram(vdw_contacts_vdmlevel, bins=int(max(vdw_contacts_vdmlevel)+1))
polars =np.histogram(polar_contacts_vdmlevel, bins=int(max(polar_contacts_vdmlevel)+1))
scores = {'hbonds':{}, 'vdws':{}, 'polars':{}}
lbl = ['hbonds', 'vdws', 'polars']
for ix, hist in enumerate([totalhbonds]):
#for ix, hist in enumerate([totalhbonds, vdws, polars]):
    score_dict = scores[lbl[ix]]
    N = np.cumsum(hist[0])
    N = sum([x<= 0.95*sum(hist[0]) for x in N])
    for count, bond_num in zip(hist[0], hist[1][:-1]):
        score_dict[bond_num] = -np.log10(count/num_ifgs_within_sasa_cutoff / (1/N))

    pkl.dump(scores, open('../Lookups/ideal_num_interactions/'+ifg+'_num_interactions_scores.pkl','wb'))

