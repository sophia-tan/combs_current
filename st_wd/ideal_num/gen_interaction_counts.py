import sys                              
import pandas as pd
import pickle as pkl
from sys import argv
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis
import numpy as np

script, ifg = argv

ifgs_dir = '/home/gpu/Sophia/STcombs/20171118/'

df = pkl.load(open('../burial/noskip_pkl/%s_noskip.pkl'%ifg,'rb'))
database = pkl.load(open('../burial/database_mindistconvexhull.pkl','rb'))
# get burial dist from database
ifgs_df, vdms_df = analysis.Analysis.add_burial_dist(df, database)

# add hbond info. merge on pdb, resnum_ifg, resname_ifg, chid_ifg.
ifgs_df = analysis.Analysis.add_interactions(ifgs_dir,ifg, ifgs_df)

#total = 0
#index = [x for x in range(11)][1:]
#for _bin in index: #_bin is the upperbound (doesn't include)
#    ifgs = ifgs_df[ifgs_df['burialdist'] >=(_bin-1)]
#    if _bin != 10:
#        ifgs = ifgs[ifgs['burialdist']<_bin]
#

#
#merged = add_hbond(df, 'hbond')
#merged = add_hbond(merged, 'ca_hbond')
#
## merged is the df that has the contacts info
#
#### start to count 
#df = merged
#df['number_hbonds'].fillna(0,inplace=True)
#df['ca_hbonds'].fillna(0,inplace=True)
#
## load in water info
#waterdf = pd.read_csv(csv_dir+ifg+'_ifg_hbond_water.csv',index_col = 'iFG_count')
## load ifg sasa info
## sasadf = blah
#sasa_cutoff = 1000000000000
#
## here, vdw means any heavy atom within 3.5A OR carbon within 4.8A of ifg and
## excludes hbonds!!! But includes polar atoms
#
## keep lists of how many times the ifg interacts on the VDM level and ATOMIC level
#water_contacts_ifglevel = [] # list of the num of water contacts an ifg makes
#hb_contacts_vdmlevel = [] # list of the num of hbond contacts an ifg makes
#vdw_contacts_vdmlevel = [] # list of the num of vdw contacts an ifg makes
#polar_contacts_vdmlevel = [] # list of the num of polar contacts an ifg makes
#ca_hb_contacts_atomlevel = []
#hb_contacts_atomlevel = []
#vdw_contacts_atomlevel = []
#polar_contacts_atomlevel = []
#totalhbond = [] # hbonds w. water AND vdm. do this on IFG LEVEL
## only makes sense to do water on ifg level, and and ca_hb atom level (since it's rare 
## to see >1 ca interaction per vdm)
#
#'''iterate through every ifg'''
#vdms_grouped_ifg = df.groupby('iFG_count')
#for ifgcount, ifg_ in vdms_grouped_ifg:
#    print(ifg_)
#    hb_vdmlevel = 0
#    vdw_vdmlevel = 0
#    polar_vdmlevel = 0
#    ca_hb_atomlevel = []
#    hb_atomlevel = []
#    vdw_atomlevel = []
#    polar_atomlevel = []
#    total_hbonds = 0
#    
##    sasa = database[
###sasadf.loc[ifgcount, 'iFG_sasa_CB_3A_probe']
##    if sasa < sasa_cutoff: 
##        # see if it interacts with water
##        if ifgcount in waterdf.index:
##            water_contacts_ifglevel.append(waterdf.loc[ifgcount, 'number_hbonds'])
##            total_hbonds += waterdf.loc[ifgcount, 'number_hbonds']
##        else:
##            water_contacts_ifglevel.append(0)
##
##        ''' iterate through every vdm'''
##        for ix, vdm in ifg_.iterrows():
##            # count hb
##            # store number of hbonds (atomlevel) as list, you can sum up the list later 
##            hb_atomlevel.append(vdm['number_hbonds'])
##            ca_hb_atomlevel.append(vdm['ca_hbonds'])
##            if vdm['number_hbonds'] > 0:
##                hb_vdmlevel += 1 # since this is on vdmlevel, can't count each vdm more than once
##                total_hbonds += vdm['number_hbonds']
##            else: # if it's not a hbond... then count it as vdw (for sure) and polar (maybe)
##                atoms = vdm['atom_names'].split(' ')
##                atoms = [a[0] for a in atoms]
##                
##                vdw_vdmlevel += 1
##                vdw_atomlevel.append(len(atoms))
##                
##                polar_atomlevel.append(atoms.count('N') + atoms.count('O'))
##                if 'O' in atoms or 'N' in atoms: # if polar 
##                    # then this vdm makes polar contacts with the ifg, so +1 vdmlevel
##                    polar_vdmlevel += 1
##        
##
##        hb_contacts_vdmlevel.append(hb_vdmlevel)
##        hb_contacts_atomlevel.append(hb_atomlevel)
##        ca_hb_contacts_atomlevel.append(ca_hb_atomlevel)
##        vdw_contacts_vdmlevel.append(vdw_vdmlevel)
##        vdw_contacts_atomlevel.append(vdw_atomlevel)
##        polar_contacts_vdmlevel.append(polar_vdmlevel)
##        polar_contacts_atomlevel.append(polar_atomlevel)
##        totalhbond.append(total_hbonds)
##
##megalist = []
##for ls in [water_contacts_ifglevel , hb_contacts_vdmlevel , vdw_contacts_vdmlevel,
##    polar_contacts_vdmlevel, hb_contacts_atomlevel, ca_hb_contacts_atomlevel, 
##    vdw_contacts_atomlevel, polar_contacts_atomlevel, totalhbond ]:
##    megalist.append(ls)
##
##pkl.dump([merged,megalist], open('interaction_types_pkl/%s_interaction_types.pkl'%ifg,'wb'))
