import sys
sys.path.append('/home/gpu/Sophia/combs/st_wd/Scoring_Function')
from ScoringFunctions import *

########################################################################
''' PROVIDE BASIC INFORMATION '''
########################################################################
lookup_dir = '/home/gpu/Sophia/combs/st_wd/Lookups/'
poses_dir = '/home/gpu/Sophia/combs/st_wd/apixaban/poses'
ifg_list = ['carboxamide', 'lonepair_imidazole', 'backboneCO']
scaffold_pdb = '1y28_bbH_0001.pdb'
ifgselection = 'chain X and resnum 1'
ligand_ifgs = {'carboxamide':[['C10', 'C11', 'N3', 'O1']],'lonepair_imidazole':[['N6']], 'backboneCO':[['C8','O3'],['C19','O2']]} 

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
