import os, sys, traceback
import pickle as pkl, prody as pr
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
import numpy as np

script, aa = sys.argv

aas = ['alanine', 'aspartate', 'glutamine', 'lysine', 'proline', \
    'asparagine', 'glutamate', 'glycine', 'phenylalanine','leucine']
look_dir = '/home/gpu/Sophia/combs/st_wd/Lookups/refinedvdms/'

coords = []
lookup = pkl.load(open(look_dir+'vdms_of_{}.pkl'.format(aa),'rb'))
for ix, row in lookup.iterrows():
    try:
        ifg = []
        vdm = []
        ifgname,vdmname = row['resname_ifg'],row['resname_vdm']
        if aa in aas:
            db_dir = '/home/gpu/Sophia/STcombs/20171118/database/reduce/'
            #db_dir = '/home/gpu/Sophia/combs/st_wd/20180207_db_molprobity_biolassem/'
            par = pr.parsePDB(db_dir + row['pdb'] + 'H.pdb')
            ifgchid,ifgresnum = row['chid_ifg'],row['resnum_ifg']
            vdmchid,vdmresnum = row['chid_vdm'],row['resnum_vdm']
            #ifgseg,vdmseg = row['segi_ifg'],row['segi_vdm']
            for atom in constants.atoms_dict[ifgname]:
                sele = 'chain {} and resnum {} and name {}'.format(ifgchid,ifgresnum,atom)
                #sele = 'chain {} and segment {} and resnum {} and name {}'.format(ifgchid,ifgseg,ifgresnum,atom)
                ifg.append(par.select(sele).getCoords()[0])
            for atom in constants.atoms_dict[vdmname]:
                sele = 'chain {} and resnum {} and name {}'.format(vdmchid,vdmresnum,atom)
                #sele = 'chain {} and segment {} and resnum {} and name {}'.format(vdmchid,vdmseg,vdmresnum,atom)
                vdm.append(par.select(sele).getCoords()[0])
        else:
            db_dir = '/home/gpu/Sophia/STcombs/20171118/database/reduce/'
            par = pr.parsePDB(db_dir + row['pdb'] + 'H.pdb')
            ifgchid,ifgresnum = row['chid_ifg'],row['resnum_ifg']
            vdmchid,vdmresnum = row['chid_vdm'],row['resnum_vdm']
            for atom in constants.atoms_dict[ifgname]:
                sele = 'chain {} and resnum {} and name {}'.format(ifgchid,ifgresnum,atom)
                ifg.append(par.select(sele).getCoords()[0])
            for atom in constants.atoms_dict[vdmname]:
                sele = 'chain {} and resnum {} and name {}'.format(vdmchid,vdmresnum,atom)
                vdm.append(par.select(sele).getCoords()[0])
        coords.append([row.name,np.array(ifg),np.array(vdm)])
    except:
                traceback.print_exc()
                coords.append([row.name,None])
pkl.dump(coords, open(look_dir+'coords_of_{}.pkl'.format(aa),'wb'))