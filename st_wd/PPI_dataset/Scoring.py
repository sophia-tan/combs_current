from sklearn.neighbors import NearestNeighbors
import sys
import traceback
import copy
from combs.apps import *
import prody as pr
import numpy as np
import pickle as pkl
import pandas as pd
from itertools import *
sys.path.append('/home/gpu/Sophia/combs/src/')

np.warnings.filterwarnings('ignore')

def get_interacting_atoms(parsed, resi, resn):
    if resi[1:].isdigit(): # removes all insertion codes
        target = parsed.select('chain %s and resnum %s_' % (resi[0], resi[1:]))
    else: # has insertion code
        target = parsed.select('chain %s and resnum %s' % (resi[0], resi[1:].upper()))

    assert len(list(set(target.getResnames()))) == 1
    assert constants.one_letter_code[target.getResnames()[0]] == resn

    # find out if interaction is BB or SC
    other_chain = parsed.select('not (chain {})'.format(resi[0], resi[1:]))
    polar = ['O', 'N']
    interacting_atoms = []
    for atom in target:
        radius = 3.5
        if atom.getName()[0] in polar:
            radius = 3.5
        else:
            radius = 4.8
        for nbr in pr.findNeighbors(atom, radius, other_chain):
            ifgatom, vdmatom, dist = nbr
            ifgindex, vdmindex = ifgatom.getIndex(), vdmatom.getIndex()
            ifgatomelem, vdmatomelem = ifgatom.getName(), vdmatom.getName()
            if ifgatom.getElement() != 'H' and vdmatom.getElement() != 'H':
                if dist <= 3.5:
                    interacting_atoms.append((ifgindex, vdmindex))
                else:
                    if ifgatomelem[0] == 'C' and vdmatomelem[0] == 'C':
                        interacting_atoms.append((ifgindex, vdmindex))
    return list(set(interacting_atoms))


def get_ifg_vdm(
        parsed,
        ifgresn,
        vdmresn,
        ifg_contact_atoms,
        vdm_contact_atoms,
        method):
    '''determine if the ifg is making bb or sc intrxns, and if the vdm is making
    bb or sc intrxns'''
    ifginfo = [
        parsed.select(
            'index %s' %
            ifg_contact_atoms[0]).getChids()[0],
        parsed.select(
            'index %s' %
            ifg_contact_atoms[0]).getResnums()[0]]
    vdminfo = [
        parsed.select(
            'index %s' %
            vdm_contact_atoms[0]).getChids()[0],
        parsed.select(
            'index %s' %
            vdm_contact_atoms[0]).getResnums()[0]]

    if isinstance(ifgresn, list):
        ifgresn = ifgresn[1]
    if isinstance(vdmresn, list):
        vdmresn = vdmresn[1]

    ifgcontactatoms = [
        parsed.select(
            'index %s' %
            x).getNames()[0] for x in ifg_contact_atoms]
    vdmcontactatoms = [
        parsed.select(
            'index %s' %
            x).getNames()[0] for x in vdm_contact_atoms]

    if method == 'planar_group':
        bb = ['N CA C', 'CA C O']
        cacbcg = ['CA CB CG']
        ifgatoms = ''
        vdmatoms = ''
        for typ in [bb, constants.planar_atoms[ifgresn], cacbcg]:
            for element in typ:
                element = element.split(' ')
                if len(
                    set(element).intersection(
                        set(ifgcontactatoms))) > len(
                    set(ifgatoms).intersection(
                        set(ifgcontactatoms))):
                    ifgatoms = element
        for typ in [bb, constants.planar_atoms[vdmresn]]:
            for element in typ:
                element = element.split(' ')
                if len(
                    set(element).intersection(
                        set(vdmcontactatoms))) > len(
                    set(vdmatoms).intersection(
                        set(vdmcontactatoms))):
                    vdmatoms = element

    elif method == 'whole_res':
        ifgatoms = constants.atoms_dict[ifgresn]
        vdmatoms = constants.atoms_dict[vdmresn]

    elif method == 'BBorSC':
        bb = ['C', 'O', 'CA', 'N']
        ifgbb = len(set(bb).intersection(set(ifgcontactatoms)))
        vdmbb = len(set(bb).intersection(set(vdmcontactatoms)))
        if (len(set(ifgcontactatoms)) - ifgbb) >= ifgbb:
            ifgatoms = constants.AA_sc_dict[ifgresn]
            if ifgresn == 'ALA':
                ifgatoms = ['CA', 'CB']
        else:
            ifgatoms = ['N', 'CA', 'C', 'O']
        if (len(set(vdmcontactatoms)) - vdmbb) >= vdmbb:
            vdmatoms = constants.AA_sc_dict[vdmresn]
            if vdmresn == 'ALA':
                vdmatoms = ['CA', 'CB']
        else:
            vdmatoms = ['N', 'CA', 'C', 'O']

    # if bb interactions, need to get combed GLY
    # only for iFG, not vdm bc no combed bbs
    if len(set(ifgatoms).intersection(['N', 'CA', 'C', 'O'])) > 1:
        ifgresn = 'GLY'

    def assign(atomlist, res):
        return [constants.AAname[res], atomlist]

    ifgtype = assign(ifgatoms, ifgresn)
    vdmtype = assign(vdmatoms, vdmresn)

    return ifgtype, vdmtype, ifginfo, vdminfo


def only_contacting(row, ifgatoms, vdmatoms):
    dist = row['dist_info']
    dist = dist.lstrip('(')
    dist = dist.rstrip(')')
    dist = dist.split(') (')
    dist = [x.split(' ') for x in dist]
    for elem in dist:
        if elem[0] in ifgatoms and elem[1] in vdmatoms:
            return True
    return False


def filter_contact(ifgresn, vdmresn, ifgatoms, vdmatoms):
    ''' keep only the vdM rows where the ifgs and vdms are
    interacting through ifgatoms and vdmatoms '''

    threecode = constants.AAname[ifgresn]
    try:
        csv = pd.read_csv(
            '/home/gpu/Sophia/STcombs/20171118/{}/csv/{}_ifg_contact_vdm.csv'.format(threecode, threecode))
    except BaseException:
        csv = pd.read_csv(
            '/home/gpu/Sophia/combs/st_wd/20180207db_combed_csvs/{}/{}_ifg_contact_vdm.csv'.format(
                threecode, threecode))
    vdmpdbinfo = pkl.load(
        open(
            '/home/gpu/Sophia/combs/st_wd/Lookups/refinedvdms/vdms_of_{}.pkl'.format(threecode),
            'rb'))
    csv = pd.merge(csv, vdmpdbinfo, left_index=True, right_index=True)
    csv = csv[csv['resname_vdm'] == vdmresn]
    num_all_vdms = len(csv) 
    csv['in_contact'] = csv.apply(
        only_contacting,
        ifgatoms=ifgatoms,
        vdmatoms=vdmatoms,
        axis=1)
    csv = csv[csv['in_contact']]
    return num_all_vdms, csv  # which is actually a pandas df

def flip(atom_ls, resn):
    # account for flipped residues
    ls = []
    ls.append(atom_ls)
    if resn in constants.flip_names:
        copyls = atom_ls.copy()
        for ix, atom in enumerate(atom_ls):
            for pair in constants.flip_names[resn]:
                if pair[0] not in atom_ls and pair[1] not in atom_ls:
                    pass
                elif pair[0] in atom_ls or pair[1] in atom_ls:
                    if atom in pair:
                        if atom == pair[0]:
                            copyls[ix] = pair[1]
                        elif atom == pair[1]:
                            copyls[ix] = pair[0]
                if copyls not in ls:
                    ls.append(copyls)
    return ls


def get_order_of_atoms(item, ifgresn, vdmresn, ifgls, vdmls):
    lookupatoms = []
    ifgorder = np.array(constants.atoms_dict[ifgresn])
    vdmorder = np.array(constants.atoms_dict[vdmresn])
    for atom in ifgls:
        index = np.where(ifgorder == atom)[0][0]
        lookupatoms.append(np.array(item[1][index]))
    for atom in vdmls:
        index = np.where(vdmorder == atom)[0][0]
        lookupatoms.append(np.array(item[2][index]))
    lookupatoms = np.array(lookupatoms)
    return lookupatoms

def score_interaction_and_dump(
        parsed,
        ifgresn,
        vdmresn,
        ifg_contact_atoms,
        vdm_contact_atoms,
        method,
        targetresi, cutoff, pdbix, pdbname):
    cutoff = float(cutoff)
    ifgtype, vdmtype, ifginfo, vdminfo = get_ifg_vdm(
        parsed, ifgresn, vdmresn, ifg_contact_atoms, vdm_contact_atoms, method)

    if ifgtype[1] != ['N', 'CA', 'C'] and ifgtype[1] != ['CA', 'C', 'O']:
        ifgresn = constants.AAname_rev[ifgtype[0]]
        vdmresn = constants.AAname_rev[vdmtype[0]]
        ifgatoms = ifgtype[1]
        vdmatoms = vdmtype[1]

        # filter for only vdmresn vdms of ifgresn with ifgatoms
        # and vdmatoms directly involved in interactions
        num_all_vdms, lookupdf = filter_contact(ifgresn, vdmresn, ifgatoms, vdmatoms)
        query = []
        for atom in ifgatoms:
            query.append(
                parsed.select(
                    'chain {} and resnum {} and name {}'.format(
                        ifginfo[0],
                        ifginfo[1],
                        atom)).getCoords()[0])
        for atom in vdmatoms:
            query.append(
                parsed.select(
                    'chain {} and resnum {} and name {}'.format(
                        vdminfo[0],
                        vdminfo[1],
                        atom)).getCoords()[0])

        query = np.array(query)
        lookupcoords = pkl.load(open(
            '/home/gpu/Sophia/combs/st_wd/Lookups/refinedvdms/coords_of_{}.pkl'.format(ifgtype[0]), 'rb'))
        #lookupcoords = lookupcoords[:50] # delete

        ifglists = flip(ifgatoms, ifgresn)
        vdmlists = flip(vdmatoms, vdmresn)
        rmsds = []
        num_atoms = len(query)
        coords_ls = [
            item for item in lookupcoords if item[0] in lookupdf.index]
        lookupatoms_to_clus = []
        counter = 0 # to keep count of how many pdbs are being output
        for item in coords_ls:
            if len(item) == 3:
                compare_rmsds = []
                ifg_vdm_ind = []
                for ifg_ind, ifgls in enumerate(ifglists):
                    for vdm_ind, vdmls in enumerate(vdmlists):
                        lookupatoms = get_order_of_atoms(
                            item, ifgresn, vdmresn, ifgls, vdmls)
                        moved, transf = pr.superpose(lookupatoms, query)
                        temp_rmsd = pr.calcRMSD(moved, query)
                        compare_rmsds.append(temp_rmsd)
                        ifg_vdm_ind.append([moved, temp_rmsd])
                # item[0] is df index
                rmsds.append([item[0], min(compare_rmsds)])
                # get index of which one had min rmsd
                for which_ind, each in enumerate(ifg_vdm_ind):
                    if each[1] == min(compare_rmsds):
                        lookupatoms_to_clus.append(each[0])
                        ######################################################################## 
                        #                   output pdb if low rmsd 
                        ######################################################################## 
                        if each[1] < cutoff and counter < 30 and which_ind==0: 
                            # this is to ensure rmsd is below cutoff when not flipped
                            # bc don't want to take care of that in prody to output pdb
                            row = lookupdf.loc[item[0]]
                            try:
                                db_dir = '/home/gpu/Sophia/STcombs/20171118/database/reduce/'
                                par = pr.parsePDB(db_dir+row['pdb']+'H.pdb')
                            except:
                                db_dir = '/home/gpu/Sophia/combs/st_wd/20180207_db_molprobity_biolassem/'
                                par = pr.parsePDB(db_dir+row['pdb']+'H.pdb')
    
                            ifgchid, ifgresnum = row['chid_ifg'], row['resnum_ifg']
                            vdmchid, vdmresnum = row['chid_vdm'], row['resnum_vdm']
                            printout = copy.deepcopy(par)
                            printout = printout.select(
                                '(chain {} and resnum {}) or (chain {} and resnum {})'.format(
                                ifgchid,ifgresnum,vdmchid,vdmresnum))
                            printout.select('chain {} and resnum {}'.format(ifgchid,ifgresnum)).setChids('Y')
                            printout.select('chain {} and resnum {}'.format(vdmchid,vdmresnum)).setChids('X')
                            printout.select('all').setResnums(10)
                            printout_interactamer = []
                            integrin_interactamer = []
                            try: # skip the ones that have segment ids. will prob need to update this 
                            # for the newly combed stuff
                                for atom in ifgatoms:
                                    integrin_interactamer.append(parsed.select('chain {} and resnum {} and name {}'.format(ifginfo[0],ifginfo[1],atom)))
                                    printout_interactamer.append(printout.select('chain Y and resnum 10 and name {}'.format(atom)))
                                for atom in vdmatoms:
                                    integrin_interactamer.append(parsed.select('chain {} and resnum {} and name {}'.format(vdminfo[0],vdminfo[1],atom)))
                                    printout_interactamer.append(printout.select('chain X and resnum 10 and name {}'.format(atom)))
                                integrin_interactamer_prody = []
    
                                integrin_interactamer = sum(integrin_interactamer[1:], integrin_interactamer[0])
                                printout_interactamer = sum(printout_interactamer[1:], printout_interactamer[0])
                                try:
                                    assert len(integrin_interactamer) == len(printout_interactamer)
    
                                    interact_res = printout.select('(chain X and resnum 10) or (chain Y and resnum 10)')
                                    interactamer_transf = pr.applyTransformation(transf, printout_interactamer)
                                    outdir = './output_data/pdbfiles/'
                                    
                                    threecode = constants.AAname[ifgresn]
    
                                    pr.writePDB(outdir+'{}_{}_{}_{}{}_{}{}_{}_{}'.format(
                                        pdbix, pdbname, targetresi,ifginfo[1],
                                        ifgresn,vdminfo[1],vdmresn,cutoff,
                                        row.name), interactamer_transf)
                                    counter += 1
                                except:
                                    pass
                            except:
                                traceback.print_exc()
                                pass

            else:
                rmsds.append([int(item[0]), 100000])

        # count how many NNs the query intrxn has
        num_nn, norm_metrics = get_NN(lookupatoms_to_clus, 
            num_atoms, rmsds, query, cutoff, num_all_vdms)
        return ifginfo[0], ifginfo[1], ifgresn, vdminfo[0], vdminfo[1],\
            vdmresn, ifgatoms, vdmatoms, num_nn, norm_metrics

def get_NN(lookupatoms_to_clus, num_atoms, rmsds, query, cutoff, num_all_vdms):
    flat = [x.reshape(-1,) for x in lookupatoms_to_clus]
    euclid_dist = cutoff * np.sqrt(num_atoms)
    nbrs = NearestNeighbors(
        n_neighbors=1,
        metric='euclidean',
        radius=euclid_dist)
    nnfit = nbrs.fit(flat)
    query = query.reshape(-1,)
    num_nn = len(nnfit.radius_neighbors(query, return_distance=False)[0])

    # See if num_nn is close to the # of geometric matches within 
    # the rmsd cutoff
    rmsds = rmsds[1:]
    geom_matches = [r for i, r in rmsds if r < cutoff]
    assert abs(len(geom_matches) - num_nn) <= 1

    ''' Metrics to normalize num. geom matches (or NNs) by: 
    1) number of vdms 
    2) number of vdms whose planar FG is making direct
       intrxns with the iFG's planar FG 
    3) expectation value for # of NNs each ifg-vdm interaction has...
       a) avg with singletons
       b) avg w/o singletons
       c) median with singletons
       d) median w/o singletons '''
    
    norm_metrics = [] 
    norm_metrics.append(num_all_vdms)
    norm_metrics.append(len(flat))
    NNs = []
    for coord in flat:
        coord = coord.reshape (1,-1)
        number_nn=len(nnfit.radius_neighbors(coord,return_distance=False)[0])
        NNs.append(number_nn)
    NNs_no_sing = [i for i in NNs if i > 1] # exclude singletons
    exp_list = [np.mean(NNs), np.mean(NNs_no_sing), 
        np.median(NNs), np.median(NNs_no_sing)]
    norm_metrics.append(exp_list)
    return num_nn, norm_metrics
