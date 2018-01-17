import sys
from ScoringFunctions import *
#from PyrosettaScores import *
import prody as pr
import freesasa, math
import numpy as np, pickle as pkl, pandas as pd
from pprint import pprint
from apixaban_info import *

script, lookup_dir = sys.argv
ligname = 'APX'
ifg_list = ['backboneCO']
#ifg_list = ['carboxamide', 'lonepair_imidazole', 'backboneCO']

####### 1) score naked backbone scaffold ####################
for scaffold_pdb in ['1y28_bbH_0001.pdb']:
    scaffold_scoredict = score_scaffold(scaffold_pdb, lookup_dir)
    # outputs dict where key=resnum, value=burial score

####### 2) score single and pairwise terms ##################

''' for every pdb in the temp_pdbs folder...'''
temp_pdbs = ['1_mem_177235_backboneCO_temp_188A_SC_iFG_631327_vdM_2_iFlip_1_backboneCO_oriented.pdb',
'1_mem_244774_carboxamide_temp_212A_SC_iFG_121240_vdM_3_iFlip_1_carboxamide_oriented.pdb',
'1_mem_254344_backboneCO_temp_199A_SC_iFG_782043_vdM_4_iFlip_1_backboneCO_oriented.pdb',
'1_mem_272364_backboneCO_temp_200A_SC_iFG_368916_vdM_2_iFlip_1_backboneCO_oriented.pdb',
'1_mem_301246_backboneCO_temp_202A_SC_iFG_707846_vdM_1_iFlip_1_backboneCO_oriented.pdb',
'1_mem_329742_backboneCO_temp_210A_SC_iFG_697140_vdM_1_iFlip_1_backboneCO_oriented.pdb',
'1_mem_65798_carboxamide_temp_176A_SC_iFG_145088_vdM_1_iFlip_1_carboxamide_oriented.pdb',
'1_mem_757085_backboneCO_temp_251A_SC_iFG_80686_vdM_3_iFlip_1_backboneCO_oriented.pdb',
'2_mem_66703_carboxamide_temp_176A_SC_iFG_104322_vdM_3_iFlip_1_carboxamide_oriented.pdb',
'3_mem_142023_backboneCO_temp_176A_SC_iFG_361743_vdM_1_iFlip_1_backboneCO_oriented.pdb']

for ligand_vdm_pair in temp_pdbs:
    for ifgname in ifg_list: ### hardcoding this. need a list called ifg_list
        ifgatomslist = ligand_ifgs[ifgname]
        for ifgatoms in ifgatomslist:
            # get freq aai score of vdm
            vdmresnum = int(ligand_vdm_pair.split('_')[5][:-1])
            bb_score, sc_score = freqaai(ligand_vdm_pair, ifgname, ifgatoms, lookup_dir)
            #phipsi_score, rotamer_score = pyrosetta_scores(ligand_vdm_pair, vdmresnum)
            if bb_score is not None or sc_score is not None:
                print(bb_score, sc_score, ligand_vdm_pair)
                interactamer_score = interactamer_geom_ligand(ligand_vdm_pair, ifgname, ligname, lookup_dir,\
                    is_bb=bb_score is not None, is_sc=sc_score is not None, ifgatoms=ifgatoms)
                #score_list = [bb_score, sc_score, phipsi_score, rotamer_score]

''' for every pdb in the temp_pdbs_pairs folder '''
temp_pdbs_pairs = ['0_0_mem_temp_pair_iFG_145088_vdM_1_iFlip_1_carboxamide_oriented.pdb_iFG_631327_vdM_2_iFlip_1_backboneCO_oriented.pdb']
for vdm_vdm_pair in temp_pdbs_pairs:
    ifgname = 'backboneCO' # which are the ifgs the vdms are interacting with?!
    coop_score = cooperativity(vdm_vdm_pair, lookup_dir, ifgname)
    score_list = [coop_score]

####### 3) score whole pose ####################################
poses_list = ['pose.pdb']
for pose in poses_list:
    for ifgname in ifg_list:
        for ifgatoms in ligand_ifgs[ifgname]:
            ideal_hbond_score = ideal_hbonds(pose, ifgname, ifgatoms, lookup_dir)

