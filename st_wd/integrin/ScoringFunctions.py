# gets cb sasa scores for each vdm of a FG for a design
# commandline arguments: pdb of design, txtfile of vdm_resnum vdm_aa 
# make sure there's a SINGLE SPACE b/n resnum and aa 
# and aa is 3 letters

import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
import prody as pr
import freesasa, math
import numpy as np, pickle as pkl, pandas as pd
from combs import analysis, parse, apps
from combs.cluster import *

#script, design_pdb, fg_vdm_txt, ligand_text= sys.argv

############################################################################################
################### make pandas df for vdms and their scores ###############################
############################################################################################
def make_df(fg_vdm_txt):
    vdms_aa = []
    vdms_resnum = []
    ifg = []
    with open(fg_vdm_txt) as inF:
        for line in inF:
            line = line.strip()
            line = line.split(' ')
            vdms_resnum.append(line[0])
            vdms_aa.append(line[1])
            ifg.append(line[2])
    columns = ['iFG', 'AA', 'BB or SC', 'sasa_cb_3a', 'burial score', 
                'dunbrack', 'rama', 'rosetta hbonds', 'f(AAi) bb', 'f(AAi) sc']
    df = pd.DataFrame(columns=columns)
    #df = pd.DataFrame(columns=columns, index=pd.Series(vdms_resnum, name='resnum'))
    # fill in df
    alphabet='abcdefg'
    for ix in range(len(vdms_aa)):
        if len(ifg[ix].split(','))>1:
            alphabet_ix = 0
            for fg in ifg[ix].split(','):
                index = str(vdms_resnum[ix])+alphabet[alphabet_ix]
                df = fill_df(df, vdms_resnum, index, ix, fg, vdms_aa)
                alphabet_ix += 1
        else: 
            index = str(vdms_resnum[ix])
            df = fill_df(df, vdms_resnum, index, ix, ifg[ix], vdms_aa)
    return df

def fill_df(df,vdms_resnum, index, ix, fg, vdms_aa):
    df.loc[index, 'iFG'] = fg
    df.loc[index, 'AA'] = vdms_aa[ix]
    return df
############################################################################################
########################### SASAs at cb 3A probe ###########################################
############################################################################################

# 1) get dictionary with vdm AAs, resnums, and sasas at cb 3A probe ################
def cb_sasas(design_pdb):
    sasa_dict = {} # key is resnum, value is list [aa, cbsasa]
    
    # parse design and do freesasa calc
    prody_parsed = pr.parsePDB(design_pdb, altloc='A', model=1)
    fs_struct = freesasa.Structure(design_pdb) # more atoms 
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
        resname = sele.getResnames()[0]
        bb_cb_atom_ind = sele.getIndices() 
        sasa_3A_probe = '{0:.2f}'.format(sum(fs_result.atomArea(i) for i in \
            np.where(np.in1d(prody_pdb_bb_cb_atom_ind,bb_cb_atom_ind))[0]))
        sasa_dict[int(resnum)] = [resname, float(sasa_3A_probe)]
    return sasa_dict

def freesasa_cb(prody_parsed, probe_radius=1.4):
    cb_sele = prody_parsed.select('protein and (backbone or name CB) and not element H D')
    coords = list(x for y in cb_sele.getCoords() for x in y)
    radii = list(freesasa.Classifier().radius(x, y) for x, y in zip(cb_sele.getResnames(), \
        cb_sele.getNames()))
    return freesasa.calcCoord(coords, radii, freesasa.Parameters({'probe-radius': probe_radius}))

# 2) 'find sasa bin for that vdm and score' ###########################

def scoring_sasa(sasadict, lookup_dir):
    score_dict = {}
    # load lookup table
    lookup = pkl.load(open(lookup_dir+'scores_largeprobe_sasa_db_lookup.pkl','rb'))
    # take -np.log10 for all values in df, and drop rows 150-180
    lookup = lookup.drop([160, 170, 180, 190])
    lookup = lookup.applymap(lambda g: -np.log10(g))
    for ix, value in sasadict.items():
        aa, sasa = value[0], value[1]
        if sasa > 150: # this site is exposed
            score = -1
        else:
            if sasa == 0:
                sasabin = 10
            else:
                sasabin = int(math.ceil(sasa / 10.0)) * 10 
            score = lookup.loc[sasabin,aa]
        score_dict[ix] = score
    rounded_dict = {k:round(v,2) for k, v in score_dict.items()} 
    return rounded_dict

# 3) wrap them together ############
def score_scaffold(design_pdb, lookup_dir):
    sasadict = cb_sasas(design_pdb)
    score_dict = scoring_sasa(sasadict, lookup_dir)
    return score_dict  
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

def freqaai(ligand_vdm_pair, ifg_name, ifgatoms, lookup_dir):
    parsed = pr.parsePDB(ligand_vdm_pair)
    
    # automatically select vdm
    v = parsed.select('chain X and resnum 10 and not element H D')
    vdm_name = v.getResnames()[0]

    # find out if interaction is BB or SC
    bb = ['C', 'O', 'OXT', 'CA', 'N']
    polar = ['O', 'N']
    ifg = parsed.select('resnum 1 and name %s'%(' '.join(ifgatoms)))
    vdmatoms = []
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
            else:
                if ifgatom[0]=='C' and vdmatom[0]=='C':
                    vdmatoms.append(vdmatom)
    vdmatoms = list(set(vdmatoms))
    bbinteraction = sum([x in bb for x in vdmatoms])
    scinteraction = sum([x not in bb for x in vdmatoms])

    # get lookup info for ifg
    pklfile = lookup_dir + 'AAi_freq_combed_%s.pkl'%ifg_name
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

    corr_pkl = pkl.load(open(lookup_dir+'noskip_%s_correlation.pkl'%ifgname, 'rb'))[0][0]
    try:
        score = corr_pkl[(vdm_a, vdm_b)]
    except:
        score = corr_pkl[(vdm_b, vdm_a)]
    return score

def exp_num_nn(ifgname, resn):
    resnpkl = pkl.load(open(''))

def get_transformed_ifgcoords(resn, ifgcoords, parsed, comb, is_bb, is_sc,vdmselection=None, ifgselection=None):
    '''helper function'''
    ifgcoords = [] # list of list. inner list is [ifgcoords, BBorSCtype]
    if is_sc:
        for key in apps.constants.interactamer_atoms[resn].keys():
            origin_atom, plane_atom1, plane_atom2 = apps.constants.interactamer_atoms[resn][key]
            coords = transform_vdmifg(resn, ifgcoords, parsed, comb, origin_atom, plane_atom1, plane_atom2, vdmselection=vdmselection, ifgselection=ifgselection)
            ifgcoords.append([coords, 'SC'])
    if is_bb:
        coords1 = transform_vdmifg(resn, ifgcoords, parsed, comb, 'C', 'O', 'CA', vdmselection=vdmselection, ifgselection=ifgselection)
        ifgcoords.append([coords1, 'C_O'])
        # turning off N H CA bc it's not named 'H' for all vdms
        #coords2 = transform_vdmifg(resn, ifgcoords, parsed, comb, 'N', 'H', 'CA')
        #ifgcoords.append([coords2, 'N_CA'])
    return ifgcoords


def transform_vdmifg(resn, ifgcoords, parsed, comb, origin_atom, plane_atom1, plane_atom2, vdmselection=None, ifgselection=None):
    #origin_atom, plane_atom1, plane_atom2 = 'C', 'O', 'CA'
    if resn != 'TYR' and resn not in apps.constants.flip_residues:
        bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs, vdmelem1avgdists = \
            make_rel_vdm_coords(parsed, comb, origin_atom, plane_atom1, plane_atom2, vdmselection=vdmselection, ifgselection=ifgselection, scoring=True)
    
    elif resn == 'TYR':
        x = make_rel_vdm_coords(parsed, comb, origin_atom, plane_atom1, plane_atom2, vdmselection=vdmselection, ifgselection=ifgselection,  scoring=True)
        y = make_rel_vdm_coords(parsed, comb, origin_atom, 'CE2', plane_atom2, vdmselection=vdmselection, ifgselection=ifgselection, unflipped=False, scoring=True)
        xelem1dist = x[8]
        yelem1dist = y[8]
        if xelem1dist <= yelem1dist:
            bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs, vdmelem1avgdists = x
        else:
            bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs, vdmelem1avgdists = y

    elif resn in apps.constants.flip_residues:
        x = make_rel_vdm_coords(parsed, comb, origin_atom, plane_atom1, plane_atom2, vdmselection=vdmselection, ifgselection=ifgselection, scoring=True)
        y = make_rel_vdm_coords(parsed, comb, origin_atom, plane_atom2, plane_atom1, vdmselection=vdmselection, ifgselection=ifgselection, unflipped=False, scoring=True)
        xelem1dist = x[8]
        yelem1dist = y[8]
        if xelem1dist <= yelem1dist:
            bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs, vdmelem1avgdists = x
        else:
            bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs, vdmelem1avgdists = y
    return ifg_coords

def get_num_nn(pdb, ifgname, ligname, lookup_dir, is_bb, is_sc, vdmselection=None, ifgselection=None, dist=None):
    ##### make sure to add the ligand's resname and atoms to interactamer_atoms in constants.py ####
    parsed = pr.parsePDB(pdb)
    comb = parse.Comb(apps.constants.ifg_sele_dict[ifgname])
    try:
        resn = parsed.select('chain X and resnum 10').getResnames()[0] # of vdm
    except:
        resn = parsed.select(vdmselection[:-10]).getResnames()[0]
    ifgcoords = [] # for HIS, there will be 2 sets
    ifgcoords = get_transformed_ifgcoords(resn, ifgcoords, parsed, comb, is_bb, is_sc,vdmselection=vdmselection, ifgselection=ifgselection)
    #  will also need to get the NN scores for sc AND bb (C_O and N_CA) if the vdm 
    # makes bb AND sc interactions. compare to see which has the best score 
    # didn't worry about any of this for integrin
    num_nn = []
    for coords in ifgcoords:
        # will have 2 for HIS
        coord, typ = coords
        coord = coord[0].reshape(1,-1)
        if dist is None:
            nnfit = pkl.load(open(lookup_dir+'fitted_NN/NNfit_%s_with_%s_%s.pkl'%(ifgname,resn,typ),'rb'))
        else:
            nnfit = pkl.load(open(lookup_dir+'fitted_NN/NNfit_%s_with_%s_%s_%s.pkl'%(ifgname,resn,typ,dist),'rb'))
        number_nn = nnfit.radius_neighbors(coord,return_distance=False)[0]
        number_nn = len(number_nn)
        num_nn.append(number_nn)
    return num_nn, resn, nnfit

def get_exp_num(typ, resn, ifgname, nnfit):
    all_coords = pkl.load(open('/home/gpu/Sophia/STcombs/20171118/%s/clusters/%s/pickle/%s_rel_vdms.pickle'%(ifgname,typ,resn),'rb'))[:,6]
    avg = []
    for coord in all_coords:
        coord = coord.reshape(1,-1)
        number_nn = nnfit.radius_neighbors(coord,return_distance=False)[0]
        number_nn = len(number_nn)
        avg.append(number_nn)
    avg = np.mean(avg)
    return avg


def interactamer_geom_ligand(pdb, ifgname, ligname, lookup_dir, is_bb, is_sc, vdmselection=None, ifgselection=None, dist=None):
    num_nn, resn, nnfit = get_num_nn(pdb, ifgname, ligname, lookup_dir, is_bb, is_sc, vdmselection=vdmselection, ifgselection=ifgselection, dist=dist)
    index = np.where(num_nn==max(num_nn))[0] # get index of the largest # in num_nn
    try:
        num_nn =num_nn[index]
    except:
        num_nn=num_nn[0]
    # determine exp_num
    if is_sc:
        typ = 'SC'
    if is_bb:
        typ='C_O'
    # !!!!!! not robust for if it'sboth SC and bb!!!! !!!!!! #
    exp = get_exp_num(typ, resn, ifgname, nnfit)
    if num_nn != 0:
        score = np.round(-np.log10(num_nn/exp),1)
    else:
        score=None
    return num_nn,np.round(exp,1),score

############################################################################################
############################ get final pose scores ############################################
############################################################################################
def find_hbonds(parsed, ifgatoms):
    ifgsele = parsed.select('chain X and resnum 1 and name %s'%(' '.join(ifgatoms)))
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

def ideal_hbonds(pose, ifgname, ifgatoms, lookup_dir):
    pklfile = pkl.load(open(lookup_dir+'%s_num_hbonds_scores.pkl'%ifgname,'rb'))
    parsed = pr.parsePDB(pose)
    hbond_list = find_hbonds(parsed,ifgatoms)
    num_hbonds = len(hbond_list)
    return pklfile[num_hbonds]

