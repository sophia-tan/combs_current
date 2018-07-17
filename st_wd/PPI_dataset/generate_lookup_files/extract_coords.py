import os, sys, traceback
import pickle as pkl, prody as pr
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
import numpy as np

script, aa, multiple_of_20 = sys.argv
'''split up by every multiple of 20. 
Ex) 20, 40, etc.'''

look_dir = '/home/gpu/Sophia/combs/st_wd/Lookups/refinedvdms/'

coords = []
lookup = pkl.load(open(look_dir+'vdms_of_{}.pkl'.format(aa),'rb'))
multiple = int(multiple_of_20) * 1000

# check to see if file exists
try:
    x=pkl.load(open(look_dir+'coords_of_{}_multiple_{}.pkl'.format(aa, multiple_of_20),'rb'))
    print('exists')
except:
    for ix, row in lookup.iterrows():
        iloc = lookup.index.get_loc(row.name) # complicated way of getting iloc
        if iloc < multiple and iloc >= multiple - 20000: 
            try:
                ifg = []
                vdm = []
                ifgname,vdmname = row['resname_ifg'],row['resname_vdm']
                if vdmname in constants.atoms_dict.keys(): # ignore MSE, etc.
                    db_dir = '/home/gpu/Sophia/combs/st_wd/20180626_db_molprobity_biolassem/'
                    par = pr.parsePDB(db_dir + row['pdb'] + 'H.pdb')
                    ifgchid,ifgresnum = row['chid_ifg'],row['resnum_ifg']
                    vdmchid,vdmresnum = row['chid_vdm'],row['resnum_vdm']
                    ifgseg,vdmseg = row['segi_ifg'],row['segi_vdm']
                    ifg_all_atoms = constants.atoms_dict[ifgname]
                    vdm_all_atoms = constants.atoms_dict[vdmname]
                    for atom in ifg_all_atoms:
                        sele = 'chain {} and segment {} and resnum {} and name {}'.format(ifgchid,ifgseg,ifgresnum,atom)
                        ifg.append(par.select(sele).getCoords()[0])
                    for atom in vdm_all_atoms:
                        sele = 'chain {} and segment {} and resnum {} and name {}'.format(vdmchid,vdmseg,vdmresnum,atom)
                        vdm.append(par.select(sele).getCoords()[0])
                    coords.append([row.name,np.array(ifg),np.array(vdm)])
                else:
                    coords.append([row.name,None])
            except:
                traceback.print_exc()
                coords.append([row.name,None])
    pkl.dump(coords, open(look_dir+'coords_of_{}_multiple_{}.pkl'.format(aa, multiple_of_20),'wb'))
