import sys
sys.path.append('/home/gpu/Sophia/combs/st_wd/Scoring_Function')
from ScoringFunctions import *

########################################################################
''' PROVIDE BASIC INFORMATION '''
########################################################################
lookup_dir = '/home/gpu/Sophia/combs/st_wd/Lookups/'
poses_dir = '/home/gpu/Sophia/combs/st_wd/designed_gln_binding_protein/gln_bp_poses'
ifg_list = ['carboxamide', 'carboxylate', 'amino']
scaffold_pdb = '1wdn_noligandH_ala_0001.pdb'
ifgselection='chain Z and resnum 1'
ligand_ifgs = {'carboxamide':[['CG', 'CD', 'OE1', 'NE2']],'carboxylate':[['CA', 'C', 'OXT', 'O']], 'amino':[['N', 'CA']]}

########################################################################
''' RUN SCORING FX FOR EACH DESIGN '''
########################################################################
# FIRST, get sasa for each residue at all ala scaffold
sasadict = cb_sasas(scaffold_pdb)
# THEN, get pairwise scores for lig-vdm and vdm-vdm
for pose in os.listdir(poses_dir):
    print('==================================================')
    print('Pose name: '+pose)
    score_poses(pose, lookup_dir, poses_dir, ifg_list, ligand_ifgs, scaffold_pdb, ifgselection, sasadict, PyR=False)
