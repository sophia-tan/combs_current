import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
from ResDicts import *
import prody as pr
import numpy as np, pickle as pkl, pandas as pd
from itertools import *
from Scoring import *
from Functions import *

script, pdb_num, inputchain, inputres, method = sys.argv
''' @pdb_num is the num of the dict from  ResDicts.py
    that we'll refer to. 
    @inputchain is the chID of the mutated residue 
    @inputres is the resnum of the mutated residue
    *the point of doing one res at a time is bc getting 
    the coords of combed vdms is slow'''

lookup_dir = '/home/gpu/Sophia/combs/st_wd/Lookups/'
db_dir = '/home/gpu/Sophia/STcombs/20171118/database/reduce/'
pdb_dict = get_pdbdict(pdb_num,list(locals().items()))

parsed = pr.parsePDB(pdb_dict[1])
for resi, resn_list in pdb_dict[2].items():
    #if 1==1:
    if resi[1:] == inputres:
        # for each target res, find the residues w/in 3.5 and 4.8A
        interacting_atoms = get_interacting_atoms(parsed, resi, resn_list[0])
        # list of list where inner list is [targetatomindex, int_resatomindex]
        interacting_resindices = set(list([parsed.select('index %s'%x[1]).getResindices()[0] for x in interacting_atoms]))
        print('--------------------')
        print('Target Res', resi, resn_list)
        for int_resindex in interacting_resindices: 
            int_res = [int_resindex,parsed.select('resindex {}'.format(int_resindex)).getResnames()[0]]
            print('Interacting Res', int_res)
            # interacting_atoms has all the indices for all the atoms the target resn
            # interacts with, but we only want to look at the ones for this int_res
            int_res_atoms = [x[1] for x in interacting_atoms if parsed.select('index %s'%x[1]).getResindices()[0] == int_resindex]
            target_res_atoms = [x[0] for x in interacting_atoms if parsed.select('index %s'%x[1]).getResindices()[0] == int_resindex]
            ######### treat targetres as ifg ############
            # int_res is resname
            print(int_res)
            if int_res[1] != 'HOH':
                score_interaction_and_dump(parsed,constants.three_letter_code[resn_list[0]], \
                    int_res, target_res_atoms, int_res_atoms,method=method,targetresi=resi, cutoff=0.5)
