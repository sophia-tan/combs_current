import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
from PyrosettaScores import *
import prody as pr
import freesasa, math, os
import numpy as np, pickle as pkl, pandas as pd
from combs import analysis, parse, apps
from combs.cluster import *

############################################################################################
########################### SASAs at cb 3A probe ###########################################
############################################################################################

# 1) get dictionary with vdm AAs, resnums, and sasas at cb 3A for all ALA ################
def cb_sasas(ala_pdb):
    sasa_dict = {} # key is resnum, value is cbsasa
    
    # parse design and do freesasa calc
    prody_parsed = pr.parsePDB(ala_pdb, altloc='A', model=1)
    fs_struct = freesasa.Structure(ala_pdb) # more atoms 
    fs_result = freesasa_cb(prody_parsed, probe_radius=3) # less atoms bc this is Cb cutoff 

    # get sasa for each resnum
    chain_start = min(prody_parsed.getResnums())
    chain_end = max(prody_parsed.getResnums())
    for resnum in range(chain_start,chain_end+1):
        # get Cb atoms 
        prody_pdb_bb_cb_atom_ind = prody_parsed.select('protein and (backbone or name CB) and \
            not element H D').getIndices()
        # get Cb atoms for resnum
        sele = prody_parsed.select('protein and (backbone or name CB) and resnum ' + str(resnum) \
            + ' and not element H D')
        bb_cb_atom_ind = sele.getIndices() 
        sasa_3A_probe = '{0:.2f}'.format(sum(fs_result.atomArea(i) for i in \
            np.where(np.in1d(prody_pdb_bb_cb_atom_ind,bb_cb_atom_ind))[0]))
        sasa_dict[int(resnum)] = float(sasa_3A_probe)
    return sasa_dict

def freesasa_cb(prody_parsed, probe_radius=3):
    cb_sele = prody_parsed.select('protein and (backbone or name CB) and not element H D')
    coords = list(x for y in cb_sele.getCoords() for x in y)
    radii = list(freesasa.Classifier().radius(x, y) for x, y in zip(cb_sele.getResnames(), \
        cb_sele.getNames()))
    return freesasa.calcCoord(coords, radii, freesasa.Parameters({'probe-radius': probe_radius}))

# 2) 'find sasa bin for that vdm and score' ###########################
def scoring_sasa(sasa, aa, lookup_dir, ifg):
    'for this vdm\'s cb_sasa and ifg, get score'
    # load lookup and db tables
    lookup = pkl.load(open(lookup_dir+'sasa_scores/scores_largeprobe_sasa_%s_lookup.pkl'%ifg,'rb'))
    db = pkl.load(open(lookup_dir+'sasa_scores/scores_largeprobe_sasa_db_lookup.pkl','rb'))
    # drop rows 150-180
    lookup = lookup.drop([160, 170, 180, 190])
    db = db.drop([160, 170, 180, 190])
    if sasa == 0:
        sasabin = 10
    elif sasa > 150: # this site is exposed. take value at 150
        sasabin = 150
    else:
        sasabin = int(math.ceil(sasa / 10.0)) * 10 
    score = -np.log10(lookup.loc[sasabin,aa]/db.loc[sasabin,aa])
    return round(score,2)

############################################################################################
############################ get freq_aai score ############################################
############################################################################################i
def get_ifgatoms(ligand_text):
    ifgdict = {}
    with open(ligand_text) as inF:
        for line in inF:
            line = line.strip()
            line = line.split(' ')
            ifg, resnum, atoms= line[0], line[1], line[2]
            atoms = atoms.split('+')
            atoms = ' '.join(atoms)
            ifgdict[ifg] = [resnum, atoms]
    return ifgdict

def get_interacting_atoms(ifg, v):
    bb = ['C', 'O', 'OXT', 'CA', 'N']
    polar = ['O', 'N']
    vdmatoms = []
    ifgatoms = []
    for atom in ifg:
        if atom.getName()[0] in polar:
            radius = 3.5
        else:
            radius = 4.8
        for nbr in pr.findNeighbors(atom, radius, v):
            ifgatom, vdmatom, dist = nbr
            ifgatom, vdmatom = ifgatom.getName(), vdmatom.getName()
            if dist <= 3.5:
                vdmatoms.append(vdmatom)
                ifgatoms.append(ifgatom)
            else:
                if ifgatom[0]=='C' and vdmatom[0]=='C':
                    vdmatoms.append(vdmatom)
                    ifgatoms.append(ifgatom)
    return ifgatoms, vdmatoms

def freqaai(pdbfile, ifg_name, ifgatoms, lookup_dir, typ=None, flipvdmvdm=False):
    parsed = pr.parsePDB(pdbfile)
    
    # automatically select vdm
    v = parsed.select('chain X and resnum 10 and not element H D')
    ifg = parsed.select('resnum 1 and name %s'%(' '.join(ifgatoms)))
    
    # but if typ==vdmvdm (instead of ligvdm), use diff chain and resnums
    if typ == 'vdmvdm':
        v = parsed.select('chain X and resnum 10 and not element H D')
        ifg = parsed.select('chain W and resnum 10 and name %s'%(' '.join(ifgatoms)))
    if flipvdmvdm==True:
        v = parsed.select('chain W and resnum 10 and not element H D')
        ifg = parsed.select('chain X and resnum 10 and name %s'%(' '.join(ifgatoms)))

    # find out if interaction is BB or SC
    bb = ['C', 'O', 'OXT', 'CA', 'N']
    vdm_name = v.getResnames()[0]
    ifgatoms, vdmatoms = get_interacting_atoms(ifg, v)
    vdmatoms = list(set(vdmatoms))
    bbinteraction = sum([x in bb for x in vdmatoms])
    scinteraction = sum([x not in bb for x in vdmatoms])

    # get lookup info for ifg
    pklfile = lookup_dir + 'AAi_freq/AAi_freq_combed_%s.pkl'%ifg_name
    lookup = pkl.load(open(pklfile,'rb'))
    
    # get database frequencies 
    db_dict = analysis.EnergyTerms.AAi_db_lookup(lookup_dir)
    if bbinteraction > 0:
        score = lookup.ix[vdm_name, 'vdm_freq_bb']
        bbscore = -np.log10(score / db_dict[vdm_name])
    else:
        bbscore = None 
    if scinteraction > 0: 
        score = lookup.ix[vdm_name, 'vdm_freq_sc']
        scscore = -np.log10(score / db_dict[vdm_name])
    else:
        scscore = None
    return bbscore, scscore


############################################################################################
############################ get vdm-vdm scores ############################################
############################################################################################

def cooperativity(vdm_vdm_pair, lookup_dir, ifgname):
    parsed = pr.parsePDB(vdm_vdm_pair)
    resnames = list(set(parsed.getResnames()))
    if len(resnames) == 1:
        vdm_a, vdm_b = resnames[0], resnames[0]
    elif len(resnames) == 2:
        vdm_a, vdm_b = resnames[0], resnames[1]

    corr_pkl = pkl.load(open(lookup_dir+'correlation/noskip_%s_correlation.pkl'%ifgname, 'rb'))[0][0]
    try:
        score = corr_pkl[(vdm_a, vdm_b)]
    except:
        score = corr_pkl[(vdm_b, vdm_a)]
    return score

def get_transformed_ifgcoords(resn, parsed, comb, is_bb, is_sc, ifgatoms=None, ifgselection=None, vdmselection=None):
    '''helper function'''
    ifgcoords = [] # list of list. inner list is [ifgcoords, BBorSCtype]
    if is_sc:
        for key in apps.constants.interactamer_atoms[resn].keys():
            origin_atom, plane_atom1, plane_atom2 = apps.constants.interactamer_atoms[resn][key]
            coords = transform_vdmifg(resn, ifgcoords, parsed, comb, origin_atom, plane_atom1, plane_atom2, ifgatoms=ifgatoms, ifgselection=ifgselection, vdmselection=vdmselection)
            ifgcoords.append([coords, 'SC'])
    if is_bb:
        coords1 = transform_vdmifg(resn, ifgcoords, parsed, comb, 'C', 'O', 'CA', ifgatoms=ifgatoms,ifgselection=ifgselection,vdmselection=vdmselection)
        ifgcoords.append([coords1, 'C_O'])
        try:
            coords2 = transform_vdmifg(resn, ifgcoords, parsed, comb, 'N', 'H', 'CA', ifgatoms=ifgatoms, ifgselection=ifgselection,vdmselection=vdmselection)
            ifgcoords.append([coords2, 'N_CA'])
        except:
            pass
    return ifgcoords


def transform_vdmifg(resn, ifgcoords, parsed, comb, origin_atom, plane_atom1, plane_atom2,ifgatoms=None,vdmselection=None,ifgselection=None):
    #origin_atom, plane_atom1, plane_atom2 = 'C', 'O', 'CA'
    if resn != 'TYR' and resn not in apps.constants.flip_residues:
        bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs, vdmelem1avgdists = \
            make_rel_vdm_coords(parsed, comb, origin_atom, plane_atom1, plane_atom2, scoring=True, ifgatoms=ifgatoms,vdmselection=vdmselection,ifgselection=ifgselection)
    
    elif resn == 'TYR':
        x = make_rel_vdm_coords(parsed, comb, origin_atom, plane_atom1, plane_atom2, scoring=True, ifgatoms=ifgatoms,vdmselection=vdmselection,ifgselection=ifgselection)
        y = make_rel_vdm_coords(parsed, comb, origin_atom, 'CE2', plane_atom2, unflipped=False, scoring=True, ifgatoms=ifgatoms,vdmselection=vdmselection,ifgselection=ifgselection)
        xelem1dist = x[8]
        yelem1dist = y[8]
        if xelem1dist <= yelem1dist:
            bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs, vdmelem1avgdists = x
        else:
            bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs, vdmelem1avgdists = y

    elif resn in apps.constants.flip_residues:
        x = make_rel_vdm_coords(parsed, comb, origin_atom, plane_atom1, plane_atom2, scoring=True,ifgatoms=ifgatoms,ifgselection=ifgselection,vdmselection=vdmselection)
        y = make_rel_vdm_coords(parsed, comb, origin_atom, plane_atom2, plane_atom1, unflipped=False, scoring=True,ifgatoms=ifgatoms,ifgselection=ifgselection,vdmselection=vdmselection)
        xelem1dist = x[8]
        yelem1dist = y[8]
        if xelem1dist <= yelem1dist:
            bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs, vdmelem1avgdists = x
        else:
            bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs, vdmelem1avgdists = y
    return ifg_coords

def interactamer_geom_ligand(pdb, ifgname, lookup_dir, is_bb, is_sc, ifgatoms=None,ifgselection=None,vdmselection=None):
    ##### make sure to add the ligand's resname and atoms to interactamer_atoms in constants.py ####

    megaifgcoords = [] # for HIS, there will be 2 sets. 4 sets for APX
    parsed = pr.parsePDB(pdb)
    if vdmselection is None:
        resn = parsed.select('chain X and resnum 10').getResnames()[0] # of vdm
    else:
        resn = parsed.select(vdmselection[:21]).getResnames()[0]
    if type(apps.constants.ifg_sele_dict[ifgname])==dict:
        comb = parse.Comb(apps.constants.ifg_sele_dict[ifgname])
        #  this is to get the NN scores for sc AND bb (C_O and N_CA) if the vdm 
        # makes bb AND sc interactions. compare to see which has the best score 
        ifgcoords = get_transformed_ifgcoords(resn, parsed, comb, is_bb, is_sc, ifgatoms=ifgatoms, ifgselection=ifgselection,vdmselection=vdmselection)
        megaifgcoords += ifgcoords
    elif type(apps.constants.ifg_sele_dict[ifgname])==list:
    # APX has more than 1 ways its resnames bc more than 1 ifg selection orders
        for ifgsel_order in apps.constants.ifg_sele_dict[ifgname]:
            comb = parse.Comb(ifgsel_order)
            ifgcoords = get_transformed_ifgcoords(resn, parsed, comb, is_bb, is_sc, ifgatoms=ifgatoms,ifgselection=ifgselection,vdmselection=vdmselection)
            megaifgcoords += ifgcoords
    nnscore = []
    for coords in megaifgcoords:
        # will have 2 for HIS
        coord, typ = coords
        coord = coord[0].reshape(1,-1)
        nnfit = pkl.load(open(lookup_dir+'fitted_NN/%s/NNfit_%s_with_%s_%s.pkl'%(ifgname,ifgname,resn,typ),'rb'))
        number_nn = max(0.01,len(nnfit.radius_neighbors(coord,return_distance=False)[0]))
        expnnfile = lookup_dir+'expected_NN/'+'expNN_%s_ifg.pkl'%ifgname
        exp_nn = pkl.load(open(expnnfile,'rb'))[(typ,resn)]
        nnscore.append(-np.log10(number_nn/exp_nn))
    return min(nnscore)
############################################################################################
############################ get final pose scores ############################################
############################################################################################
def find_hbonds(parsed, ifgatoms, ifgselection):
    ifgsele = parsed.select(ifgselection+' and name %s'%(' '.join(ifgatoms)))
    # find any hbond donors
    acc_hyd_don_list = []
    ifgdon = ifgsele.select('element O N')
    if ifgdon is not None:
        for donor in ifgdon.iterAtoms():
            hyds = pr.Contacts(parsed)(1.1, donor.getCoords()).select('element H D')
            if hyds is not None:
                for hyd in hyds.iterAtoms():
                    acceptors = pr.Contacts(parsed)(2.5, hyd.getCoords()).select('element O N \
                        and not (chain X and resnum 1)')
                    if acceptors is not None:
                        for acceptor in acceptors.iterAtoms():
                            acc_hyd_don_list.append([acceptor, hyd, donor])
    
    # find any hbond acceptors
    ifgacc = ifgsele.select('element O or (resname HIS and sidechain and element N)')
    if ifgacc is not None:
        for acceptor in ifgacc.iterAtoms():
            hyds = pr.Contacts(parsed)(2.5, acceptor.getCoords()).select('element H D')
            if hyds is not None:
                for hyd in hyds.iterAtoms():
                    donors = pr.Contacts(parsed)(1.1, hyd.getCoords()).select('element O N \
                        and not (chain X and resnum 1)')
                    if donors is not None:
                        for donor in donors.iterAtoms():
                            acc_hyd_don_list.append([acceptor,hyd,donor])
    hbond_list = check_hbond_angles(acc_hyd_don_list)
    return hbond_list

def check_hbond_angles(acc_hyd_don_list):
    newlist = []
    for triplet in acc_hyd_don_list:
        acc, hyd, don = triplet
        angle = pr.calcAngle(acc, hyd, don)
        if angle > 90:
            newlist.append(triplet)
    return newlist

def ideal_hbonds(pose, ifgname, ifgatoms, lookup_dir,ifgselection):
    pklfile = pkl.load(open(lookup_dir+'ideal_num_interactions/%s_num_interactions_scores.pkl'%ifgname,'rb'))['hbonds']
    parsed = pr.parsePDB(pose)
    hbond_list = find_hbonds(parsed,ifgatoms, ifgselection)
    num_hbonds = len(hbond_list)
    return pklfile[num_hbonds]

def get_rosettaresnum(scaffold_pdb, residue_num):
    # vdmresnum might be offset, depending on if pdb starts numbering at 1 or not
    parsed = pr.parsePDB(scaffold_pdb)
    first_resnum = parsed.getResnums()[0]
    rosettavdmresnum = int(residue_num[:-1]) - first_resnum + 1
    return rosettavdmresnum

def get_cooperativity_score(first_vdm_res, second_vdm_res, pair, resdict, pair_pdb, lookup_dir):
    first_vdm_id = str(int(pair.split('_')[0]) + 1)
    second_vdm_id = str(int(pair.split('_')[1])  + 1)
    first_vdm_info = resdict[(first_vdm_res, first_vdm_id)]
    second_vdm_info = resdict[(second_vdm_res, second_vdm_id)]
    first_vdm_ifgs = [x[0] for x in first_vdm_info if x[1]==True]
    second_vdm_ifgs = [x[0] for x in second_vdm_info if x[1]==True]
    # see which ifgs are shared by both vdms
    cooperative_ifgs = set([x[0] for x in first_vdm_ifgs if x in second_vdm_ifgs])
    cooperativity_scores = [] # one for each ifg shared by both vdms
    for coop_ifg in cooperative_ifgs:
        coop_score = cooperativity(pair_pdb, lookup_dir, coop_ifg)
        cooperativity_scores.append(coop_score)
    if len(cooperativity_scores) > 0:
        return min(cooperativity_scores)
    else:
        return None

def get_ifgnames(ifgatoms, ifgres):
    bb = ['C', 'O', 'OXT', 'CA', 'N']
    interactions = [] # list of lists were sublist[0] = ifgname, sublist[1] = ifgatoms
    for ifgname, val in apps.constants.ifg_sele_dict.items():
        sc_int, bb_int = False, False
        if ifgres in val:
            for atom in ifgatoms:
                if atom in val[ifgres].split(' '):
                    sc_int = True
                if atom in bb:
                    bb_int = True
        if sc_int:
            interactions.append([ifgname, val[ifgres].split(' ')])
        if bb_int:
            interactions.append(['backboneCO', 'C O']) # change to bb when it's done running!
    return interactions

def get_paired_vdm_scores(first_vdm_res, second_vdm_res,pair_pdb,sasadict,lookup_dir):
    parsed = pr.parsePDB(pair_pdb)
    first = parsed.select('chain W and resnum 10 and not element H D') 
    second = parsed.select('chain X and resnum 10 and not element H D') 
    first_sel = 'chain W and resnum 10'
    second_sel = 'chain X and resnum 10'
    firstresnum = int(first_vdm_res[:-1])
    secondresnum = int(second_vdm_res[:-1])
    sasa_scores = []
    freq_aai_scores = []
    interactamer_scores = []
    # will have at least 2 scores because first is vdm1 as ifg, second is vdm2 as ifg.
    # more than 2 scores means that a vdm is making bb AND sc interactions.  
    # (max 4 scores bc bb/bb, bb/sc, sc/bb, sc/sc)
    for ix, order in enumerate([[first,second],[second,first]]): 
        ifgatoms, vdmatoms = get_interacting_atoms(order[0], order[1])
        if len(ifgatoms)==0:
            return None, None, None
        # determine ifgatoms and ifgname to feed into freqaai fx
        ifgatoms = list(set(ifgatoms))
        ifgres = order[0].getResnames()[0]
        interactions = get_ifgnames(ifgatoms,ifgres) # returns list bc could be both bb and sc
        for interxn in interactions:
            ifgname, atomslist = interxn # called atomoslist instead of ifgatoms bc confusing
            if ix == 0:
                bbscscores = freqaai(pair_pdb, ifgname, atomslist, lookup_dir, typ='vdmvdm')
                if bbscscores[0]==None and bbscscores[1] == None:
                    return None, None, None
                sasa_scores.append(scoring_sasa(sasadict[firstresnum], first.getResnames()[0], lookup_dir, ifgname))
                interactamer_scores.append(interactamer_geom_ligand(pair_pdb, ifgname, lookup_dir, bbscscores[0] is \
                    not None,bbscscores[1] is not None, ifgselection=first_sel, vdmselection=second_sel+' and name '))
            if ix == 1:
                bbscscores = freqaai(pair_pdb, ifgname, atomslist, lookup_dir, typ='vdmvdm', flipvdmvdm=True)
                if bbscscores[0]==None and bbscscores[1] == None:
                    return None, None, None
                sasa_scores.append(scoring_sasa(sasadict[secondresnum], second.getResnames()[0], lookup_dir, ifgname))
                interactamer_scores.append(interactamer_geom_ligand(pair_pdb, ifgname, lookup_dir, bbscscores[0] is \
                    not None,bbscscores[1] is not None, ifgselection=second_sel, vdmselection=first_sel+ ' and name '))
            for scoretype in bbscscores: 
                if scoretype is not None:
                    freq_aai_scores.append(scoretype)
    return np.mean(sasa_scores), np.mean(freq_aai_scores),np.mean(interactamer_scores)

############################################################################################
######################## use above functions to score poses ################################
############################################################################################
def score_poses(pose, lookup_dir, poses_dir, ifg_list, ligand_ifgs, scaffold_pdb, ifgselection, sasadict, PyR=False):
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
                        if PyR==True:
                            phipsi_score, rotamer_score = pyrosetta_scores(lig_vdm_pdb,scaffold_pdb, rosettavdmresnum)
                        interactamer_score = interactamer_geom_ligand(lig_vdm_pdb, ifgname, lookup_dir,\
                                is_bb=bb_score is not None, is_sc=sc_score is not None, ifgatoms=ifgatoms, ifgselection=ifgselection)
                        if residue_num == '157A':
                            print(interactamer_score, vdm, ifgname, ifgatoms)
    # 2) SCORE VDM-VDM PAIRWISE INTERACTIONS
    temp_pdbs_pair_dir = poses_dir + '/' + pose + '/temp_pdbs_pair/'
    for vdm_vdm_dir in os.listdir(temp_pdbs_pair_dir):
        first_vdm_res = vdm_vdm_dir.split('_')[1]
        second_vdm_res = vdm_vdm_dir.split('_')[2]
        for pair in os.listdir(temp_pdbs_pair_dir+vdm_vdm_dir):
            pair_pdb = temp_pdbs_pair_dir+vdm_vdm_dir+'/'+pair
            coop_score = get_cooperativity_score(first_vdm_res, second_vdm_res, pair, resdict, pair_pdb, lookup_dir)
            sasa_score,freq_aai_score,interactamer_score = get_paired_vdm_scores(first_vdm_res, second_vdm_res,pair_pdb, sasadict, lookup_dir)
    
    # 3) SCORE WHOLE POSE
    pose = poses_dir + '/' + pose + '/pose.pdb'
    hbondscores_list = []
    for ifgname in ifg_list:
        for ifgatoms in ligand_ifgs[ifgname]:
            ideal_hbond_score = ideal_hbonds(pose, ifgname, ifgatoms, lookup_dir, ifgselection)
            hbondscores_list.append(ideal_hbond_score)
    hbondscore = np.mean(hbondscores_list)
