# need to uncomment last line and do [T,F]
import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
sys.path.append('/home/gpu/Sophia/combs/st_wd/Scoring_Function')
from ScoringInteractions import *
from combs.apps import *
from AiibB3_residues import *
import prody as pr, numpy as np, pickle as pkl, pandas as pd

script, inputres, method = sys.argv 
''' @inputres: the resnum from residues_integrin.py '''

parsed = pr.parsePDB('pymolintegrin.pdb')
# For each target res, find the residues w/in 3.5 and 4.8A
for resi, resn in integrin_res.items():
    list_for_df = [] # will ultimately dump this df in output dir
    #if 1==1:
    if resi[1:] == inputres:
        interacting_atoms = get_interacting_atoms(parsed, resi, resn)
        # list of list where inner list is [targetatomindex, int_resatomindex]
        targetresindex = interacting_atoms[0][0]
        interacting_resindices = set(list([parsed.select('index %s'%x[1]).
            getResindices()[0] for x in interacting_atoms]))
        print('--------------------')
        print('Target Res', resi, resn)
        
        # keep only the interacting atoms that interact w/ this int_res
        # (instead of all interacting atoms that interact w/ whole target res)
        for int_resindex in interacting_resindices: 
            int_res = [int_resindex, parsed.select('resindex {}'.format(
                int_resindex)).getResnames()[0]]
            print('Interacting Res', int_res)
            int_res_atoms = [x[1] for x in interacting_atoms if parsed.select(
                'index %s'%x[1]).getResindices()[0] == int_resindex]
            target_res_atoms = [x[0] for x in interacting_atoms if parsed.select(
                'index %s'%x[1]).getResindices()[0] == int_resindex]
            
            # ===== get # geometric matches to score ======
            # Reminder: int_res is resname
            rmsd_list = [.3,.4,.5,.6,.7,.8]
            ls = [[targetresindex, constants.three_letter_code[resn]], int_res, 
                target_res_atoms, int_res_atoms]

            # first, treat hotspot residue as ifg. then, switch and treat it as vdm
            assignment_of_ifg = [ls, [ls[1], ls[0], ls[3], ls[2]]]

            for bool_one in [True, False]: # for filtering non_membrane proteins from db
                for bool_two in [True, False]: # for filtering buried ifgs/vdms
                    if method == 'planar_group_ala_mutants':
                        # when switching assignment of ifg/vdm, switch bb consideration status
                        methods_ls = ['planar_group_no_bb_for_ifg', 'planar_group_no_bb_for_vdm']
                    else: 
                        methods_ls = [method, method]

                    # take turns treating hotspot/interacting residues as ifg
                    for order_of_ifg, specific_method, switch_status in zip(
                            assignment_of_ifg, methods_ls, ['not_switched', 'switched']):
                        mega_ls = score_interaction_and_dump(parsed, order_of_ifg[0][1], 
                            order_of_ifg[1], order_of_ifg[2], order_of_ifg[3], 
                            method=specific_method, targetresi=resi, 
                            cutoff_ls=rmsd_list, output_pdb=True, non_membrane=bool_one, 
                            buried=bool_two) # mega_ls is a list of all results for all rmsd cutoffs
                                
                        if mega_ls != None:
                            for matches, rmsd in zip(mega_ls, rmsd_list):
                                ifgchid, ifgresi, ifgresn, vdmchid, vdmresi, vdmresn, ifgatoms, \
                                    vdmatoms, num_nn, norm_metrics = matches
                                list_for_df.append(pd.Series([resi, rmsd, ifgresi, ifgresn, 
                                    vdmresi, vdmresn, ifgatoms, vdmatoms, num_nn, 
                                    norm_metrics[0], norm_metrics[1], norm_metrics[2][0], 
                                    norm_metrics[2][1], norm_metrics[2][2], norm_metrics[2][3],
                                    specific_method, bool_one, bool_two, switch_status], index = [
                                    'hotspot res', 'rmsd', 'ifg resi', 
                                    'ifg resn', 'vdm resi', 'vdm resn', 
                                    'ifg atoms', 'vdm atoms', 'num NN', 'num vdms', 
                                    'num directly interacting vdms', 'avg num NNs', 
                                    'avg num NNs w/o singles', 'median num NNs', 
                                    'median num NNs w/o singles', 'method', 'nonmembrane', 
                                    'buried', 'switched ifg/vdm?']))
        outputdf = pd.DataFrame(list_for_df)
        pkl.dump(outputdf,open('./output_data/{}_matches_{}.pkl'.
            format(resi, method),'wb'))
