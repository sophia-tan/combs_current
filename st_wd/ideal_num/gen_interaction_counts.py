import sys                              
import pandas as pd
import pickle as pkl
from sys import argv
import numpy as np

script, ifg = argv

csv_dir = '/home/gpu/Sophia/STcombs/20171118/%s/csv/'%ifg

df = pkl.load(open('../sasa/noskip_sasa_pkl/noskip_%s_sasa.pkl'%ifg,'rb'))
df = df[['iFG_count', 'pdb', 'resname_ifg', 'resnum_ifg', 'chid_ifg', 'vdM_count', 'resname_vdm', 'resnum_vdm', 'chid_vdm', 'vdM_sasa_dssp_full_residue', 'vdM_sasa_CB_3A_probe']]

# add atom names and hbond info
atoms = pd.read_csv(csv_dir+'%s_vdm_pdb_info.csv'%ifg)
atoms = atoms[['iFG_count', 'vdM_count', 'atom_names']]
merged = pd.merge(df,atoms, on=['iFG_count', 'vdM_count'])

# add vdm hbond info
def add_hbond(merged, hb):
    ''' hb can be 'hbond' or 'ca_hbond' '''
    vdmhbond = pd.read_csv(csv_dir+'%s_ifg_%s_vdm.csv'%(ifg, hb))[['iFG_count', 'vdM_count', 'rel_resnums', 'number_hbonds']] # number_hbonds is only hbonds with vdm
    if hb == 'ca_hbond':
        vdmhbond = vdmhbond.rename(columns={'number_hbonds':'ca_hbonds'})
        vdmhbond = vdmhbond.loc[vdmhbond['rel_resnums'].isin(['0', '00', '000', '0000'])][['iFG_count', 'vdM_count', 'ca_hbonds']]
    else:
        vdmhbond = vdmhbond.loc[vdmhbond['rel_resnums'].isin(['0', '00', '000', '0000'])][['iFG_count', 'vdM_count', 'number_hbonds']]
    merged = pd.merge(merged, vdmhbond, on=['iFG_count','vdM_count'], how='left')
    return merged

merged = add_hbond(merged, 'hbond')
merged = add_hbond(merged, 'ca_hbond')

# merged is the df that has the contacts info

### start to count 
df = merged
df['number_hbonds'].fillna(0,inplace=True)
df['ca_hbonds'].fillna(0,inplace=True)

# load in water info
waterdf = pd.read_csv(csv_dir+ifg+'_ifg_hbond_water.csv',index_col = 'iFG_count')
# load ifg sasa info
sasadf = pd.read_csv(csv_dir+ifg+'_ifg_atom_density.csv', index_col = 'iFG_count')
sasa_cutoff = 1000000000000

# here, vdw means any heavy atom within 3.5A OR carbon within 4.8A of ifg and
# excludes hbonds!!! But includes polar atoms

# keep lists of how many times the ifg interacts on the VDM level and ATOMIC level
water_contacts_ifglevel = [] # list of the num of water contacts an ifg makes
hb_contacts_vdmlevel = [] # list of the num of hbond contacts an ifg makes
vdw_contacts_vdmlevel = [] # list of the num of vdw contacts an ifg makes
polar_contacts_vdmlevel = [] # list of the num of polar contacts an ifg makes
ca_hb_contacts_atomlevel = []
hb_contacts_atomlevel = []
vdw_contacts_atomlevel = []
polar_contacts_atomlevel = []
totalhbond = [] # hbonds w. water AND vdm. do this on IFG LEVEL
# only makes sense to do water on ifg level, and and ca_hb atom level (since it's rare 
# to see >1 ca interaction per vdm)

'''iterate through every ifg'''
vdms_grouped_ifg = df.groupby('iFG_count')
for ifgcount, ifg_ in vdms_grouped_ifg:
    hb_vdmlevel = 0
    vdw_vdmlevel = 0
    polar_vdmlevel = 0
    ca_hb_atomlevel = []
    hb_atomlevel = []
    vdw_atomlevel = []
    polar_atomlevel = []
    total_hbonds = 0
    
    sasa = sasadf.loc[ifgcount, 'iFG_sasa_CB_3A_probe']
    if sasa < sasa_cutoff: 
        # see if it interacts with water
        if ifgcount in waterdf.index:
            water_contacts_ifglevel.append(waterdf.loc[ifgcount, 'number_hbonds'])
            total_hbonds += waterdf.loc[ifgcount, 'number_hbonds']
        else:
            water_contacts_ifglevel.append(0)

        ''' iterate through every vdm'''
        for ix, vdm in ifg_.iterrows():
            # count hb
            # store number of hbonds (atomlevel) as list, you can sum up the list later 
            hb_atomlevel.append(vdm['number_hbonds'])
            ca_hb_atomlevel.append(vdm['ca_hbonds'])
            if vdm['number_hbonds'] > 0:
                hb_vdmlevel += 1 # since this is on vdmlevel, can't count each vdm more than once
                total_hbonds += vdm['number_hbonds']
            else: # if it's not a hbond... then count it as vdw (for sure) and polar (maybe)
                atoms = vdm['atom_names'].split(' ')
                atoms = [a[0] for a in atoms]
                
                vdw_vdmlevel += 1
                vdw_atomlevel.append(len(atoms))
                
                polar_atomlevel.append(atoms.count('N') + atoms.count('O'))
                if 'O' in atoms or 'N' in atoms: # if polar 
                    # then this vdm makes polar contacts with the ifg, so +1 vdmlevel
                    polar_vdmlevel += 1
        

        hb_contacts_vdmlevel.append(hb_vdmlevel)
        hb_contacts_atomlevel.append(hb_atomlevel)
        ca_hb_contacts_atomlevel.append(ca_hb_atomlevel)
        vdw_contacts_vdmlevel.append(vdw_vdmlevel)
        vdw_contacts_atomlevel.append(vdw_atomlevel)
        polar_contacts_vdmlevel.append(polar_vdmlevel)
        polar_contacts_atomlevel.append(polar_atomlevel)
        totalhbond.append(total_hbonds)

megalist = []
for ls in [water_contacts_ifglevel , hb_contacts_vdmlevel , vdw_contacts_vdmlevel,
    polar_contacts_vdmlevel, hb_contacts_atomlevel, ca_hb_contacts_atomlevel, 
    vdw_contacts_atomlevel, polar_contacts_atomlevel, totalhbond ]:
    megalist.append(ls)

pkl.dump([merged,megalist], open('interaction_types_pkl/%s_interaction_types.pkl'%ifg,'wb'))
