import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
from ScoringFunctions import *
#from PyrosettaScores import *
from residues_integrin import *
import prody as pr
import numpy as np, pickle as pkl, pandas as pd
from pprint import pprint
from itertools import *
script, lookup_dir = sys.argv


######## 1) score naked backbone scaffold ####################
#for scaffold_pdb in ['1y28_bbH_0001.pdb']:
#    scaffold_scoredict = score_scaffold(scaffold_pdb, lookup_dir)
#    # outputs dict where key=resnum, value=burial score
#
######## 2) score single and pairwise terms ##################
#
#''' for every pdb in the temp_pdbs folder...'''
#temp_pdbs = ['1_mem_177235_backboneCO_temp_188A_SC_iFG_631327_vdM_2_iFlip_1_backboneCO_oriented.pdb',
#'1_mem_244774_carboxamide_temp_212A_SC_iFG_121240_vdM_3_iFlip_1_carboxamide_oriented.pdb',
#'1_mem_254344_backboneCO_temp_199A_SC_iFG_782043_vdM_4_iFlip_1_backboneCO_oriented.pdb',
#'1_mem_272364_backboneCO_temp_200A_SC_iFG_368916_vdM_2_iFlip_1_backboneCO_oriented.pdb',
#'1_mem_301246_backboneCO_temp_202A_SC_iFG_707846_vdM_1_iFlip_1_backboneCO_oriented.pdb',
#'1_mem_329742_backboneCO_temp_210A_SC_iFG_697140_vdM_1_iFlip_1_backboneCO_oriented.pdb',
#'1_mem_65798_carboxamide_temp_176A_SC_iFG_145088_vdM_1_iFlip_1_carboxamide_oriented.pdb',
#'1_mem_757085_backboneCO_temp_251A_SC_iFG_80686_vdM_3_iFlip_1_backboneCO_oriented.pdb',
#'2_mem_66703_carboxamide_temp_176A_SC_iFG_104322_vdM_3_iFlip_1_carboxamide_oriented.pdb',
#'3_mem_142023_backboneCO_temp_176A_SC_iFG_361743_vdM_1_iFlip_1_backboneCO_oriented.pdb']
#
#for ligand_vdm_pair in temp_pdbs:
#    scoredict = {} # keys=ifg, values=scorelist
#    for ifgname in ['carboxamide']:
#    #for ifgname in ifg_list: ### hardcoding this. need a list called ifg_list
#        ifgatoms = ligand_ifgs[ifgname]
#        # get freq aai score of vdm
#        vdmresnum = int(ligand_vdm_pair.split('_')[5][:-1])
#        bb_score, sc_score = freqaai(ligand_vdm_pair, ifgname, ifgatoms, lookup_dir)
#        #phipsi_score, rotamer_score = pyrosetta_scores(ligand_vdm_pair, vdmresnum)
#        if bb_score is not None or sc_score is not None:
#            print(bb_score, sc_score, ligand_vdm_pair)
#            interactamer_score = interactamer_geom_ligand(ligand_vdm_pair, ifgname, ligname, lookup_dir,\
#                is_bb=bb_score is not None, is_sc=sc_score is not None)
#            print(interactamer_score)
#            #score_list = [bb_score, sc_score, phipsi_score, rotamer_score]
#            #scoredict[ifgname] = score_list
#
#''' for every pdb in the temp_pdbs_pairs folder '''
#temp_pdbs_pairs = ['0_0_mem_temp_pair_iFG_145088_vdM_1_iFlip_1_carboxamide_oriented.pdb_iFG_631327_vdM_2_iFlip_1_backboneCO_oriented.pdb']
#for vdm_vdm_pair in temp_pdbs_pairs:
#    ifgname= 'carboxamide' #### hardcoded this in, need to have a variable for this ####
#    coop_score = cooperativity(vdm_vdm_pair, lookup_dir, ifgname)
#    score_list = [coop_score]
#
######## 3) score whole pose ####################################
#poses_list = ['pose.pdb']
#for pose in poses_list:
#    ideal_hbond_score = ideal_hbonds(pose, ifgname, ligand_ifgs[ifgname], lookup_dir)
#

ifg_naming = {}
ifg_naming['Carboxamide'] = 'carboxamide'
ifg_naming['Carboxylate'] = 'carboxylate'
ifg_naming['PhenolOH'] = 'tyrCOH'
ifg_naming['Alcohol'] = 'thrCOH'
ifg_naming['Guano'] = 'guanidino'

def get_interacting_atoms(parsed, resi, resn):    
    target = parsed.select('chain %s and resnum %s'%(resi[0], resi[1:]))
    assert len(list(set(target.getResnames()))) == 1
    assert constants.one_letter_code[target.getResnames()[0]] == resn

    # find out if interaction is BB or SC
    other_chain = parsed.select('not chain %s'%resi[0])
    bb = ['C', 'O', 'OXT', 'CA', 'N']
    polar = ['O', 'N']
    interacting_atoms = []
    #print(resi, resn)
    for atom in target:
        radius = 3.5
        #if atom.getName()[0] in polar:
        #    radius = 3.5
        #else:
        #    radius = 4.8
        for nbr in pr.findNeighbors(atom, radius, other_chain):
            ifgatom, vdmatom, dist = nbr
            ifgindex, vdmindex = ifgatom.getIndex(), vdmatom.getIndex()
            ifgatomelem, vdmatomelem = ifgatom.getName(), vdmatom.getName()
            if dist <= 3.5:
                interacting_atoms.append((ifgindex,vdmindex))
            #else:
            #    if ifgatomelem[0]=='C' and vdmatomelem[0]=='C':
            #        interacting_atoms.append(vdmindex)
            #        print(ifgatom, vdmatom, vdmatom.getResname())
    return list(set(interacting_atoms))

def get_ifg(int_res, int_res_atoms, targetresn, target_res_atoms):
    '''logic for selecting where the ifg is, and which residue it came from.
    this logic only works for this selection! not strict enough for other things
    look at target residue/atoms and the interacting residue/atoms it interacts with 
    (target = mutated res)
    if Y594, vdm is phenyl and ifg is the int_res's iFG
    if target is making backboneCO interaction, take the interacting residue's iFG
    else take target's ifg 
    '''
    intresatoms = [parsed.select('index %s'%x).getNames()[0] for x in int_res_atoms]
    if resi == 'B594':
        vdmname = 'T'
        vdm = parsed.select('chain B and resnum 603')
        ifgname = 'carboxamide'
        ifg = parsed.select('chain A and resnum 753')
        return ifgname, vdmname, ifg, vdm, target_res_atoms, int_res_atoms
    if resi == 'A785':
        ifgname = 'carboxylate'
        ifg = parsed.select('chain A and resnum 785')
        vdmname = 'Y'
        vdm = parsed.select('chain B and resnum 594')
        return ifgname, vdmname, ifg, vdm, target_res_atoms, int_res_atoms
    if 'O' in intresatoms:
        ifgname = 'backboneCO'
        vdmname = targetresn
        ifgresnum = [parsed.select('index %s'%x).getResnums()[0] for x in int_res_atoms][0]
        ifg = parsed.select('not chain %s and resnum %s'%(resi[0], ifgresnum))
        vdmresnum = [parsed.select('index %s'%x).getResnums()[0] for x in target_res_atoms][0]
        vdm = parsed.select('chain %s and resnum %s'%(resi[0],vdmresnum))
        return ifgname, vdmname, ifg, vdm, int_res_atoms, target_res_atoms
    if 'O' in target_res_atoms:
        ifgname = 'backboneCO'
        vdmname = int_res
        ifgresnum = [parsed.select('index %s'%x).getResnums()[0] for x in target_res_atoms][0]
        ifg = parsed.select('chain %s and resnum %s'%(resi[0], ifgresnum))
        vdmresnum = [parsed.select('index %s'%x).getResnums()[0] for x in int_res_atoms][0]
        vdm = parsed.select('not chain %s and resnum %s'%(resi[0],vdmresnum))
        return ifgname, vdmname, ifg, vdm, target_res_atoms, int_res_atoms
    else:
        ifgname = list(constants.interactamer_atoms[int_res].keys())[0]
        vdmname = targetresn
        ifgresnum = [parsed.select('index %s'%x).getResnums()[0] for x in int_res_atoms][0]
        ifg = parsed.select('not chain %s and resnum %s'%(resi[0], ifgresnum))
        vdmresnum = [parsed.select('index %s'%x).getResnums()[0] for x in target_res_atoms][0]
        vdm = parsed.select('chain %s and resnum %s'%(resi[0],vdmresnum))
        return ifgname, vdmname, ifg, vdm, int_res_atoms, target_res_atoms

    #except:
    #    targetres = constants.three_letter_code[targetresn]
    #    ifgname = list(constants.interactamer_atoms[targetres].keys())[0]
    #    vdmname = constants.one_letter_code[int_res]
    #    ifgresnum = [parsed.select('index %s'%x).getResnums()[0] for x in target_res_atoms][0]
    #    ifg = parsed.select('chain %s and resnum %s'%(resi[0], ifgresnum))
    #    vdmresnum = [parsed.select('index %s'%x).getResnums()[0] for x in int_res_atoms][0]
    #    vdm = parsed.select('not chain %s and resnum %s'%(resi[0],vdmresnum))
    #    return ifgname, vdmname, ifg, vdm, target_res_atoms, int_res_atoms

#### START CODE ####
bb = ['C', 'O', 'OXT', 'CA', 'N']
parsed = pr.parsePDB('integrin.pdb')
# get database frequencies 
db_dict = analysis.EnergyTerms.AAi_db_lookup(lookup_dir)
geomdict = {}
for resi, resn in integrin_res.items():
    geomdict[resi] = []
    # for each target res, find the residues w/in 3.5 and 4.8A
    interacting_atoms = get_interacting_atoms(parsed, resi, resn)
    # list of list where inner list is [ifgindex, vdmindex]
    # can get unique interacting residues by resname bc no repeat resname/resnum
    interacting_res = set(list([parsed.select('index %s'%x[1]).getResnames()[0] for x in interacting_atoms]))
    print(resi, resn)
    for int_res in interacting_res:
        int_res_atoms = [x[1] for x in interacting_atoms if parsed.select('index %s'%x[1]).getResnames()[0] == int_res]
        target_res_atoms = [x[0] for x in interacting_atoms if parsed.select('index %s'%x[1]).getResnames()[0] == int_res]
        ifgname, vdmname, ifg, vdm, ifgatoms, vdmatoms = get_ifg(int_res, int_res_atoms, resn, target_res_atoms)
        ###if ifgname in ifg_naming.keys():
        ###    ifgname = ifg_naming[ifgname]
        ###lookup = pkl.load(open(lookup_dir+'AAi_freq_combed_%s.pkl'%ifgname,'rb'))  
        ###vdmatoms = ' '.join([str(x) for x in vdmatoms])
        ###vdmatoms = parsed.select('index %s'%vdmatoms).getNames()
        ###bbscore,scscore = None, None
        ###vdmatoms = list(set(vdmatoms))
        ###bbinteraction = sum([x in bb for x in vdmatoms])
        ###scinteraction = sum([x not in bb for x in vdmatoms])
        ###vdmname = constants.three_letter_code[vdmname]
        ###if bbinteraction > 0: 
        ###    score = lookup.ix[vdmname, 'vdm_freq_bb']
        ###    bbscore = -np.log10(score / db_dict[vdmname])
        ###if scinteraction > 0: 
        ###    score = lookup.ix[vdmname, 'vdm_freq_sc']
        ###    scscore = -np.log10(score / db_dict[vdmname])
        ####print(bbscore, scscore)
        ###vdminfo = [vdm.getChids()[0], vdm.getResnums()[0]]
        ###ifginfo = [ifg.getChids()[0], ifg.getResnums()[0]]
        ###nnlist = []
        ###for d in [1, 1.2, 1.4, 1.6, 1.8]:
        ###    num_nn,expnum,geomscore = interactamer_geom_ligand('integrin.pdb', ifgname, ifg.getResnames()[0], lookup_dir, \
        ###    is_bb=bbscore is not None, is_sc=scscore is not None, vdmselection='chain %s and resnum %s and name ' \
        ###    %(vdminfo[0], vdminfo[1]), ifgselection='chain %s and resnum %s'%(ifginfo[0], ifginfo[1]), dist=d)
        ###    nnlist.append(num_nn)
        ###geomdict[resi].append([int_res, nnlist])
#pkl.dump(geomdict, open('geomdict.pkl','wb'))
    ### get correlation score
    #if len(interacting_res) != 0:
    #    combo=combinations(interacting_res,2)
    #    for pair in combo:
    #        print(pair)
    #        print(ifgname)


