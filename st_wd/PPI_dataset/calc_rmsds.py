import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
import prody as pr
import numpy as np, pickle as pkl, pandas as pd
from itertools import *
sys.path.append('/home/gpu/Sophia/combs/st_wd/Scoring_Function')
from ScoringInteractions import *
from PPI_Functions import *

script, pdb_index, method = sys.argv
''' @pdb_index is to distinguish pdbs in the dataset
     (index of element in the pickle file lists that contain
     that pdb) 
    @method could be 'sc' or 'planar_group_no_bb'. We want to ignore ifg bb 
     interactions bc we're looking at mutations to ala, and bb interactions 
     don't change '''

try:
    pkl.load(open('./output_data/pdb_ix_{}_matches.pkl'.
        format(pdb_index),'rb'))
except:
    mutation_data = pkl.load(open('mutation_data.pkl','rb'))
    pdb_resdict = mutation_data[0][int(pdb_index)]
    pdbname, pdb_resdict = pdb_resdict[0], pdb_resdict[1]
                    
    list_for_df = [] # will ultimately dump this df in output dir
    parsed = pr.parsePDB(pdbname)
    parsed = parsed.select('not water')
    print('#####################################')
    print(pdb_index, pdbname)
    # For each target res, find the residues w/in 3.5 and 4.8A
    for resi, resn_list in pdb_resdict.items():
        interacting_atoms = get_interacting_atoms(parsed, resi, resn_list[0])
        # list of list where inner list is [targetatomindex, int_resatomindex]
        if interacting_atoms != None:
            interacting_resindices = set(list([parsed.select('index %s'%
                x[1]).getResindices()[0] for x in interacting_atoms]))
            print('--------------------')
            print('Target Res', resi, resn_list)
            # Keep only the interacting atoms that interact w/ this int_res
            # (instead of all interacting atoms that interact w/ whole target res)
            for int_resindex in interacting_resindices: 
                int_res = [int_resindex, parsed.select('resindex {}'.format(
                    int_resindex)).getResnames()[0]]
                print('Interacting Res', int_res)
                int_res_atoms = [x[1] for x in interacting_atoms if parsed.select(
                    'index %s'%x[1]).getResindices()[0] == int_resindex]
                target_res_atoms = [x[0] for x in interacting_atoms if parsed.select(
                    'index %s'%x[1]).getResindices()[0] == int_resindex]
                # ======= treat targetres as ifg ==========
                # Reminder: int_res is resname
                if int_res[1] in constants.planar_atoms.keys(): 
                    rmsd_list = [.2,.3,.4,.5,.6,.7,.8]
                    for bool_one in [True, False]: # filtering membrane proteins from db
                        for bool_two in [True, False]: # filtering buried/exposed ifgs/vdms
                            mega_ls = score_interaction_and_dump(parsed,
                                constants.three_letter_code[resn_list[0]], int_res, target_res_atoms, 
                                int_res_atoms,method=method,targetresi=resi, cutoff_ls=rmsd_list, 
                                pdbix=pdb_index, pdbname=pdbname, non_membrane=bool_one,
                                buried=bool_two, output_pdb=False)
                            if mega_ls != None:
                                for matches, rmsd in zip(mega_ls, rmsd_list):
                                    matches = mega_ls[ls_index]
                                    ifgchid, ifgresi, ifgresn, vdmchid, vdmresi, vdmresn, ifgatoms, \
                                        vdmatoms, num_nn, norm_metrics = matches
                                    list_for_df.append(pd.Series([pdb_index, pdbname, resn_list[0], 
                                        resn_list[1], resn_list[2], rmsd, ifgchid, ifgresi, ifgresn, 
                                        vdmchid, vdmresi, vdmresn, ifgatoms, vdmatoms, num_nn, 
                                        norm_metrics[0], norm_metrics[1], norm_metrics[2][0], 
                                        norm_metrics[2][1], norm_metrics[2][2], norm_metrics[2][3],
                                        method, bool_one, bool_two], index = ['pdb ix', 'pdb id', 'WT res', 
                                        'mutated res', 'mut string', 'rmsd', 'ifg chid', 
                                        'ifg resi', 
                                        'ifg resn', 'vdm chid', 'vdm resi', 'vdm resn', 
                                        'ifg atoms', 'vdm atoms', 'num NN', 'num vdms', 
                                        'num directly interacting vdms', 'avg num NNs', 
                                        'avg num NNs w/o singles', 'median num NNs', 
                                        'median num NNs w/o singles', 'method', 'nonmembrane', 'buried']))
            if len(interacting_resindices) == 0:
                for rmsd in [.2,.3,.4,.5,.6,.7,.8]:
                    list_for_df.append(pd.Series([pdb_index, pdbname, resn_list[0], 
                        resn_list[1], resn_list[2], rmsd, resi[0], resi[1:], 
                        resn_list[0],
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, method, None, None], 
                        index = ['pdb ix', 'pdb id', 'WT res', 'mutated res', 
                        'mut string', 'rmsd', 
                        'ifg chid', 'ifg resi', 'ifg resn', 'vdm chid', 'vdm resi', 
                        'vdm resn', 'ifg atoms', 'vdm atoms', 'num NN', 'num vdms', 
                        'num directly interacting vdms', 'avg num NNs', 
                        'avg num NNs w/o singles', 'median num NNs', 
                        'median num NNs w/o singles', 'method', 'nonmembrane', 'buried']))
    outputdf = pd.DataFrame(list_for_df)
    pkl.dump(outputdf,open('./output_data/pdb_ix_{}_matches_{}.pkl'.
        format(pdb_index, method),'wb'))
