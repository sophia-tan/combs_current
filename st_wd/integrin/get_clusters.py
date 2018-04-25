import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
#from ScoringFunctions import *
#from PyrosettaScores import *
from residues_integrin import *
import prody as pr
import numpy as np, pickle as pkl, pandas as pd
from pprint import pprint
from itertools import *
from Scoring import *

#script, inputres = sys.argv # inputres is resnum
# from residues_integrin.py so we can do one AA at a time
# for getting the coords of combed vdms, which is slow

lookup_dir = '/home/gpu/Sophia/combs/st_wd/Lookups/'
db_dir = '/home/gpu/Sophia/combs/st_wd/20180207_db_molprobity_biolassem/'

bb = ['C', 'O', 'OXT', 'CA', 'N']
parsed = pr.parsePDB('integrin.pdb')
geomdict = {}
for resi, resn in integrin_res.items():
    if 1==1:
    #if resi[1:] == inputres:
        # for each target res, find the residues w/in 3.5 and 4.8A
        interacting_atoms = get_interacting_atoms(parsed, resi, resn)
        # list of list where inner list is [targetatomindex, int_resatomindex]
        interacting_resindices = set(list([parsed.select('index %s'%x[1]).getResindices()[0] for x in interacting_atoms]))
        
        if len(interacting_resindices) > 0:
            geomdict[(resi,resn)] = {}
            print('--------------------')
            print('Target Res', resi, resn)
            for int_resindex in interacting_resindices: 
                int_res = [int_resindex,parsed.select('resindex {}'.format(int_resindex)).getResnames()[0]]
                print('Interacting Res', int_res)
                # interacting_atoms has all the indices for all the atoms the target resn interacts with, 
                # but we only want to look at the ones for this int_res
                int_res_atoms = [x[1] for x in interacting_atoms if parsed.select('index %s'%x[1]).getResindices()[0] == int_resindex]
                target_res_atoms = [x[0] for x in interacting_atoms if parsed.select('index %s'%x[1]).getResindices()[0] == int_resindex]
                
                ######### first, treat targetres (chA) as ifg, and then treat int_res as ifg ########
                # 1) targetres as ifg
                # int_res is resname
                ifgtype, vdmtype, ifginfo, vdminfo = get_ifg_vdm(parsed,constants.three_letter_code[resn], \
                    int_res, target_res_atoms, int_res_atoms)
                chA_as_ifg = get_clusters(resi,parsed,ifgtype,vdmtype,ifginfo,vdminfo,lookup_dir, db_dir, \
                    ifgs=True)

                # 2) int_res as ifg
                ifgtype, vdmtype, ifginfo, vdminfo = get_ifg_vdm(parsed,int_res, \
                    constants.three_letter_code[resn], int_res_atoms, target_res_atoms)
                chB_as_ifg = get_clusters(resi,parsed,ifgtype,vdmtype,ifginfo,vdminfo,lookup_dir, \
                    db_dir, ifgs=True)
                geomdict[(resi,resn)][int_res[0]] = [int_res[1],chA_as_ifg, chB_as_ifg]

#pkl.dump(geomdict, open('ifgs_integrin_geom_dict.pkl','wb')) # each dict val is list of lists...look in plot.py for more details
