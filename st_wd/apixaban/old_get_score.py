from ScoringFunctions import *
#from PyrosettaScores import *
import freesasa, math, os
import numpy as np, pickle as pkl, pandas as pd

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

    # 1) SCORE LIGAND-VDM PAIRWISE INTERACTIONS
    ''' key = (vdm resnum, id, ifg). value = [[first ifg interaction info], [second ifg interaction info], etc.]
    Order of ifgs is determined by ligand_info.py. Example...
    {(vdm resnum, vdm ID, ifg): [[interacts with ifg True/False, ifgatoms],[interacts with ifg True/False, ifgatoms]]} '''
    resdict = {}
    temp_pdbs_dir = poses_dir + '/' + pose + '/temp_pdbs/'
    for residue_dir in os.listdir(temp_pdbs_dir):
        residue_num = residue_dir.split('_')[1]
        for vdm in os.listdir(temp_pdbs_dir+residue_dir): 
            vdm_id = vdm.split('_')[0]
            for ifgname in ifg_list: 
                ifgatomslist = ligand_ifgs[ifgname]
                for ifgatoms in ifgatomslist: # ex) separates amino1 and amino2 if there are 2 aminos in ligand 
                    # get freq aai score of vdm
                    lig_vdm_pdb = temp_pdbs_dir+residue_dir+'/'+vdm
                    bb_score, sc_score = freqaai(lig_vdm_pdb, ifgname, ifgatoms, lookup_dir)
                    try:
                        resdict[(residue_num, vdm_id)].append([(ifgname, ifgatoms), bb_score is not None or sc_score is not None])
                    except:
                        resdict[(residue_num, vdm_id)] = [[(ifgname, ifgatoms), bb_score is not None or sc_score is not None]]
                        
                    rosettavdmresnum = get_rosettaresnum(scaffold_pdb, residue_num)
                    if bb_score is not None or sc_score is not None:
                        # Rosetta score includes x0.45 for phi/psi score, and x0.7 for dunbrack score
                        #####phipsi_score, rotamer_score = pyrosetta_scores(lig_vdm_pdb,scaffold_pdb, rosettavdmresnum)
                        interactamer_score = interactamer_geom_ligand(lig_vdm_pdb, ifgname, lookup_dir,\
                                is_bb=bb_score is not None, is_sc=sc_score is not None, ifgatoms=ifgatoms, ifgselection=ifgselection)
    # 2) SCORE VDM-VDM PAIRWISE INTERACTIONS
    temp_pdbs_pair_dir = poses_dir + '/' + pose + '/temp_pdbs_pair/'
    for vdm_vdm_dir in os.listdir(temp_pdbs_pair_dir):
        first_vdm_res = vdm_vdm_dir.split('_')[1]
        second_vdm_res = vdm_vdm_dir.split('_')[2]
        for pair in os.listdir(temp_pdbs_pair_dir+vdm_vdm_dir):
            pair_pdb = temp_pdbs_pair_dir+vdm_vdm_dir+'/'+pair
            coop_score = get_cooperativity_score(first_vdm_res, second_vdm_res, pair, resdict, pair_pdb, lookup_dir)
            sasa_score,freq_aai_score,interactamer_score = get_paired_vdm_scores(first_vdm_res, second_vdm_res,pair_pdb, sasadict, lookup_dir, ifgselection=ifgselection)
    
    
    # 3) SCORE WHOLE POSE
    pose = poses_dir + '/' + pose + '/pose.pdb'
    hbondscores_list = []
    for ifgname in ifg_list:
        for ifgatoms in ligand_ifgs[ifgname]:
            ideal_hbond_score = ideal_hbonds(pose, ifgname, ifgatoms, lookup_dir)
            hbondscores_list.append(ideal_hbond_score)
    hbondscore = np.mean(hbondscores_list)
