# troubleshoot by doing just one target resi

import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
import prody as pr
import numpy as np, pickle as pkl, pandas as pd
from itertools import *
from Scoring import *
from Functions import *

script, pdb_index, inputres, method = sys.argv
''' @pdb_index is to distinguish pdbs in the dataset
     (index of element in the pickle file lists that contain
     that pdb) '''

mutation_data = pkl.load(open('mutation_data.pkl','rb'))
lookup_dir = '/home/gpu/Sophia/combs/st_wd/Lookups/'
db_dir = '/home/gpu/Sophia/STcombs/20171118/database/reduce/'
pdb_resdict = mutation_data[0][int(pdb_index)]
pdb_datadf = mutation_data[1][int(pdb_index)]
pdbname, pdb_resdict = pdb_resdict[0], pdb_resdict[1]
                
list_for_df = [] # will ultimately dump this df in output dir
parsed = pr.parsePDB(pdbname)
for resi, resn_list in pdb_resdict.items():
    if resi== inputres:
        print('yas')
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
            if int_res[1] != 'HOH' and resn_list[0] not in ['W','C','M']: # need to fix for newly combed
            #if int_res[1] != 'HOH':
                for rmsd in [.4,.5]:
                    matches = score_interaction_and_dump(parsed,
                        constants.three_letter_code[resn_list[0]], int_res, target_res_atoms, 
                        int_res_atoms,method=method,targetresi=resi, cutoff=rmsd)
                    if matches != None:
                        ifgchid, ifgresi, ifgresn, vdmchid, vdmresi, vdmresn, ifgatoms, \
                            vdmatoms, num_nn, norm_metrics = matches
                        list_for_df.append(pd.Series([pdb_index, pdbname, resn_list[0], 
                            resn_list[1], resn_list[2], rmsd, ifgchid, ifgresi, ifgresn, 
                            vdmchid, vdmresi, vdmresn, ifgatoms, vdmatoms, num_nn, 
                            norm_metrics[0], norm_metrics[1], norm_metrics[2][0], 
                            norm_metrics[2][1], norm_metrics[2][2], norm_metrics[2][3],
                            method], index = ['pdb ix', 'pdb id', 'WT res', 
                            'mutated res', 'mut string', 'rmsd', 'ifg chid', 
                            'ifg resi', 
                            'ifg resn', 'vdm chid', 'vdm resi', 'vdm resn', 
                            'ifg atoms', 'vdm atoms', 'num NN', 'num vdms', 
                            'num directly interacting vdms', 'avg num NNs', 
                            'avg num NNs w/o singles', 'median num NNs', 
                            'median num NNs w/o singles', 'method']))
        if len(interacting_resindices) == 0:
            for rmsd in [.4,.5]:
                list_for_df.append(pd.Series([pdb_index, pdbname, resn_list[0], 
                    resn_list[1], resn_list[2], rmsd, resi[0], resi[1:], 
                    resn_list[0],
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, method], 
                    index = ['pdb ix', 'pdb id', 'WT res', 'mutated res', 
                    'mut string', 'rmsd', 
                    'ifg chid', 'ifg resi', 'ifg resn', 'vdm chid', 'vdm resi', 
                    'vdm resn', 'ifg atoms', 'vdm atoms', 'num NN', 'num vdms', 
                    'num directly interacting vdms', 'avg num NNs', 
                    'avg num NNs w/o singles', 'median num NNs', 
                    'median num NNs w/o singles', 'method']))
