import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
import prody as pr
import numpy as np, pickle as pkl, pandas as pd
from itertools import *
from clusterScoring import *
from Functions import *

pdb_index, method = 3, 'planar_group'
mutation_data = pkl.load(open('mutation_data.pkl','rb'))
pdb_resdict = mutation_data[0][int(pdb_index)]
pdbname, pdb_resdict = pdb_resdict[0], pdb_resdict[1]
                
parsed = pr.parsePDB(pdbname)
print('#####################################')
print(pdb_index, pdbname)
for resi, resn_list in pdb_resdict.items():
    # for each target res, find the residues w/in 3.5 and 4.8A
    interacting_atoms = get_interacting_atoms(parsed, resi, resn_list[0])
    # list of list where inner list is [targetatomindex, int_resatomindex]
    interacting_resindices = set(list([parsed.select('index %s'%x[1]).getResindices()[0] 
        for x in interacting_atoms]))
    print('--------------------')
    print('Target Res', resi, resn_list)
    for int_resindex in interacting_resindices: 
        int_res = [int_resindex,parsed.select('resindex {}'.format(int_resindex)).getResnames()[0]]
        print('Interacting Res', int_res)
        # interacting_atoms has all the indices for all the atoms the target resn
        # interacts with, but we only want to look at the ones for this int_res
        int_res_atoms = [x[1] for x in interacting_atoms if parsed.select('index %s'%x[1]).
            getResindices()[0] == int_resindex]
        target_res_atoms = [x[0] for x in interacting_atoms if parsed.select('index %s'%x[1])
            .getResindices()[0] == int_resindex]
        ######### treat targetres as ifg ############
        # int_res is resname
        if int_res[1] != 'HOH' and resn_list[0] not in ['W','C','M'] \
            and int_res[1] in constants.planar_atoms.keys(): # need to fix for newly combed
        #if int_res[1] != 'HOH':
            for rmsd in [.4]:
                matches = score_interaction_and_dump(parsed,
                    constants.three_letter_code[resn_list[0]], int_res, target_res_atoms, 
                    int_res_atoms,method=method,targetresi=resi, cutoff=rmsd, 
                    pdbix=pdb_index, pdbname=pdbname)
