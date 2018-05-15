# for getting clusters for figures 

import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
from residues_integrin import *
import prody as pr
import numpy as np, pickle as pkl, pandas as pd
from pprint import pprint
from itertools import *
from Scoring import *

script, inputres = sys.argv
#script, inputres, method = sys.argv # inputres is resnum
# from residues_integrin.py so we can do one AA at a time
# for getting the coords of combed vdms, which is slow

lookup_dir = '/home/gpu/Sophia/combs/st_wd/Lookups/'
db_dir = '/home/gpu/Sophia/STcombs/20171118/database/reduce/'

bb = ['C', 'O', 'OXT', 'CA', 'N']
parsed = pr.parsePDB('pymolintegrin.pdb')
geomdict = {}
for resi, resn in integrin_res.items():
    #if 1==1:
    if resi[1:] == inputres:
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
                #print('Interacting Res', int_res)
                # interacting_atoms has all the indices for all the atoms the target resn interacts with, 
                # but we only want to look at the ones for this int_res
                int_res_atoms = [x[1] for x in interacting_atoms if parsed.select('index %s'%x[1]).getResindices()[0] == int_resindex]
                target_res_atoms = [x[0] for x in interacting_atoms if parsed.select('index %s'%x[1]).getResindices()[0] == int_resindex]
                ######### treat targetres (chA) as ifg ############
                # int_res is resname
                output_pdbs(parsed,constants.three_letter_code[resn], \
                    int_res, target_res_atoms, int_res_atoms,method='planar_group',targetresi=resi)
