import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
sys.path.append('/home/gpu/Sophia/combs/st_wd/Scoring_Function')
from ScoringInteractions import *
from combs.apps import *
from residues_integrin import *
import prody as pr, numpy as np, pickle as pkl, pandas as pd

script, inputres, method = sys.argv 
''' @inputres: the resnum from residues_integrin.py '''

parsed = pr.parsePDB('pymolintegrin.pdb')
# For each target res, find the residues w/in 3.5 and 4.8A
for resi, resn in integrin_res.items():
    #if 1==1:
    if resi[1:] == inputres:
        interacting_atoms = get_interacting_atoms(parsed, resi, resn)
        # list of list where inner list is [targetatomindex, int_resatomindex]
        interacting_resindices = set(list([parsed.select('index %s'%x[1]).
            getResindices()[0] for x in interacting_atoms]))
        print('--------------------')
        print('Target Res', resi, resn)
        
        # keep only the interacting atoms that interact w/ this int_res
        # (instead of all interacting atoms that interact w/ whole target res)
        for int_resindex in interacting_resindices: 
            int_res = [int_resindex,parsed.select('resindex {}'.format(
                int_resindex)).getResnames()[0]]
            print('Interacting Res', int_res)
            int_res_atoms = [x[1] for x in interacting_atoms if parsed.select(
                'index %s'%x[1]).getResindices()[0] == int_resindex]
            target_res_atoms = [x[0] for x in interacting_atoms if parsed.select(
                'index %s'%x[1]).getResindices()[0] == int_resindex]
            
            # ======= treat targetres as ifg ==========
            # Reminder: int_res is resname
            for rmsd in [.3,.4,.5]:
                matches = score_interaction_and_dump(parsed,
                    constants.three_letter_code[resn], int_res, target_res_atoms, 
                    int_res_atoms,method=method,targetresi=resi, cutoff=rmsd)
