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
    parsed = parsed.select('not element H D and not name OXT')
    if resi[1:].isdigit(): # removes all insertion codes
        target = parsed.select('chain %s and resnum %s_' % (resi[0], resi[1:]))
    else: # has insertion code
        target = parsed.select('chain %s and resnum %s' % (resi[0], resi[1:].upper()))

    assert len(list(set(target.getResnames()))) == 1
    assert constants.one_letter_code[target.getResnames()[0]] == resn

    # make sure that the target res doesn't have a segment id
    # or anything funky or missing density
    if len(target) != len(constants.AA_sc_dict[constants.one_to_three[resn]]) + 4:
        print('Find out why there are missing or extra atoms in the prody res sele!')
    else: 
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

    if method == 'planar_group' or method == 'planar_group_no_bb':
        bb = ['N CA C', 'CA C O']
        cacbcg = ['CA CB CG']
        ifgatoms = ''
        vdmatoms = ''
        if method == 'planar_group':
            possible_ifgatoms = [bb, constants.planar_atoms[ifgresn], cacbcg]
        elif method == 'planar_group_no_bb':
            possible_ifgatoms = [constants.planar_atoms[ifgresn], cacbcg]
        # ifg intrxn cannot be bb for analyzing ala mutants, but vdm can be
        possible_vdmatoms = [bb, constants.planar_atoms[vdmresn], cacbcg]
        # Determine ifgatoms 
        for typ in possible_ifgatoms:
            for element in typ:
                element = element.split(' ')
                # if element has more atoms in ifgcontactatoms than current %ifgatoms%:
                if len(
                    set(element).intersection(
                        set(ifgcontactatoms))) > len(
                    set(ifgatoms).intersection(
                        set(ifgcontactatoms))):
                    ifgatoms = element
        # Determine vdmatoms 
        for typ in possible_vdmatoms:
            for element in typ:
                element = element.split(' ')
                if len(
                    set(element).intersection(
                        set(vdmcontactatoms))) > len(
                    set(vdmatoms).intersection(
                        set(vdmcontactatoms))):
                    vdmatoms = element

    # this restricts ifgatoms to just the sc (bc it will be mutated away in ala 
    # scanning), but it can interact w/ bb or sc in vdm
    elif method == 'sc_only_ifg': 
        ifgatoms = constants.AA_sc_dict[ifgresn]
        
        bb = ['C', 'O', 'CA', 'N']
        vdmbb = len(set(bb).intersection(set(vdmcontactatoms))) 
        if (len(set(vdmcontactatoms))-vdmbb) >= vdmbb:
            vdmatoms = constants.AA_sc_dict[vdmresn]
            if vdmresn == 'ALA':
                vdmatoms = ['CA', 'CB']
        else: 
            vdmatoms = ['N', 'CA', 'C', 'O']

    elif method == 'whole_res':
        ifgatoms = constants.atoms_dict[ifgresn]
        vdmatoms = constants.atoms_dict[vdmresn]

    def assign(atomlist, res):
        return [constants.AAname[res], atomlist]

    ifgtype = assign(ifgatoms, ifgresn)
    vdmtype = assign(vdmatoms, vdmresn)

    bb = ['C', 'O', 'CA', 'N']
    if ifgatoms == '' and set(ifgcontactatoms).issubset(bb) \
        and method == 'planar_group_no_bb':
        pass
    elif ifgatoms == '' or vdmatoms == '':
        print('Investigate why %ifgatoms% or %vdmatoms% is empty!')

    return ifgtype, vdmtype, ifginfo, vdminfo


def only_contacting(row, ifgatoms, vdmatoms):
    con = row['probe_contact_pairs']
    con = con.lstrip('(')
    con = con.rstrip(')')
    con = con.split(') (')
    con = [x.split(' ') for x in con]
    for elem in con:
        if elem[0] in ifgatoms and elem[1] in vdmatoms and elem[2] != 'wc':
            return True
    return False

def filter_buried(csv, threecode):
    ''' choosing arbitrary min. convex hull dist of 3A 
    as cutoff to define buried vs exposed ifgs/vdms'''
    ifgsasa = pd.read_csv(
            '/home/gpu/Sophia/combs/st_wd/20180626_combed_csvs/{}/{}_ifg_atom_density.csv'.
            format(threecode, threecode))[['iFG_count', 'min_hull_dist_CB_CA']]
    ifgsasa.rename(columns={'min_hull_dist_CB_CA': 'iFG_min_hull_dist'}, inplace=True)
    vdmsasa = pd.read_csv(
            '/home/gpu/Sophia/combs/st_wd/20180626_combed_csvs/{}/{}_vdm_sasa_info.csv'.
            format(threecode, threecode))[['iFG_count', 'vdM_count', 'min_hull_dist_CB_CA']]
    vdmsasa.rename(columns={'min_hull_dist_CB_CA': 'vdM_min_hull_dist'}, inplace=True)

    print(len(csv), 'with solvent exposed')
    csv = pd.merge(csv, ifgsasa, on=['iFG_count'], left_index=True)
    csv = pd.merge(csv, vdmsasa, on=['iFG_count', 'vdM_count'],left_index=True)
    csv = csv[csv['iFG_min_hull_dist'] > 3]
    csv = csv[csv['vdM_min_hull_dist'] > 3]
    print(len(csv), 'only buried')
    return csv

def filter_membrane(csv, threecode):
    print(len(csv), 'including membrane proteins')
    
    membrane = []
    with open('/home/gpu/Sophia/combs/st_wd/Lookups/opm_entries.txt') as inF:
        for line in inF:
            membrane.append(line.strip())

    csv = csv[~csv['pdb'].isin(membrane)]
    print(len(csv), 'without membrane proteins')
    return csv
    

def filter_contact(ifgresn, vdmresn, ifgatoms, vdmatoms, 
    exclude_membrane, exclude_exposed):
    '''Find out which are the vdM rows where the ifgs and vdms are interacting 
    through ifgatoms and vdmatoms. Exclude probe contact type 'wc'.
    Returns both filtered and unfiltered df'''

    threecode = constants.AAname[ifgresn]
    csv = pd.read_csv(
            '/home/gpu/Sophia/combs/st_wd/20180626_combed_csvs/{}/{}_ifg_contact_vdm.csv'.
            format(threecode, threecode))[['iFG_count', 'vdM_count', 'ifg_probe_contacts', 'probe_contact_pairs']]
    
    vdmpdbinfo = pkl.load(
        open(
        '/home/gpu/Sophia/combs/st_wd/Lookups/refinedvdms/vdms_of_{}.pkl'.
        format(threecode),'rb'))
    csv = pd.merge(csv, vdmpdbinfo, left_index=True, right_index=True)
    csv = csv[csv['resname_vdm'] == vdmresn]
    
    if exclude_membrane==True:
        csv = filter_membrane(csv, threecode)

    if exclude_exposed==True:
        csv = filter_buried(csv, threecode)
    
    num_all_vdms = len(csv) 
    unfiltered_df = csv
    csv['in_contact'] = csv.apply(
        only_contacting,
        ifgatoms=ifgatoms,
        vdmatoms=vdmatoms,
        axis=1)
    csv = csv[csv['in_contact']]
    num_direct = len(csv)
    return num_all_vdms, num_direct, unfiltered_df, \
        csv # csv is just filtered  pandas df

def flip(atom_ls, resn):
    # account for flipped residues
    ls = []
    ls.append(atom_ls)
    resn = constants.AAname[resn]
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
        parsed, ifgresn, vdmresn, ifg_contact_atoms, vdm_contact_atoms,
        method, targetresi, cutoff, pdbix=None, pdbname=None, output_pdb=False,
        non_membrane=False, buried=False):

    cutoff = float(cutoff)
    ifgtype, vdmtype, ifginfo, vdminfo = get_ifg_vdm(
        parsed, ifgresn, vdmresn, ifg_contact_atoms, vdm_contact_atoms, method)

    ifgresn = constants.AAname_rev[ifgtype[0]]
    vdmresn = constants.AAname_rev[vdmtype[0]]
    ifgatoms = ifgtype[1]
    vdmatoms = vdmtype[1]
    
    if ifgatoms == '' or vdmatoms == '':
        pass # probably has bb atoms only
    elif set(ifgatoms).issubset(set(constants.atoms_dict[ifgresn])) and set(
        vdmatoms).issubset(set(constants.atoms_dict[vdmresn])):

        # To normalize the # of geom matches, count only combed %vdmresn% 
        # vdms of ifgresn. And then keep only
        # the combed vdms whose %vdmatoms% are directly involved in the combed 
        # interaction. (Reminder: %vdmatoms% are determined from the 
        # target/query structure.) This process is to disregard, for example, 
        # combed vdms interacting through its bb residues if the target 
        # structure (and therefore, %vdmatoms%) don't have interacting bb residues
        num_all_vdms, num_direct, lookupdf, filtered_df = filter_contact(
            ifgresn, vdmresn, ifgatoms, vdmatoms, exclude_membrane=non_membrane,
            exclude_exposed=buried)
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
        # I forgot, but len(item)==2 if there's a problem w/ extract_coords?
        # Also, only keep the coords if the vdm is in lookupdf
        coords_ls = [
            item for item in lookupcoords if (item[0] in lookupdf.index 
                and len(item)==3)]
            
        lookupatoms_to_clus = []
        counter = 0 # to keep count of how many pdbs are being output
        for item in coords_ls:
            compare_rmsds = []
            ifg_vdm_ind = []
            for ifg_ind, ifgls in enumerate(ifglists):
                for vdm_ind, vdmls in enumerate(vdmlists):
                    lookupatoms = get_order_of_atoms(
                        item, ifgresn, vdmresn, ifgls, vdmls)
                    moved, transf = pr.superpose(lookupatoms, query)
                    temp_rmsd = pr.calcRMSD(moved, query)
                    compare_rmsds.append(temp_rmsd)
                    ifg_vdm_ind.append([moved, transf, temp_rmsd])
            # item[0] is df index
            rmsds.append([item[0], min(compare_rmsds)])
            # get index of which one had min rmsd
            for which_ind, each in enumerate(ifg_vdm_ind):
                if each[2] == min(compare_rmsds):
                    lookupatoms_to_clus.append(each[0])
                    ######################################################################## 
                    #                   output pdb if low rmsd 
                    ######################################################################## 
                    if each[2] < cutoff and counter < 30 and which_ind==0 and output_pdb==True: 
                        # this is to ensure rmsd is below cutoff when not flipped
                        # bc don't want to take care of that in prody to output pdb
                        row = lookupdf.loc[item[0]]
                        db_dir = '/home/gpu/Sophia/combs/st_wd/20180207_db_molprobity_biolassem/'
                        par = pr.parsePDB(db_dir+row['pdb']+'H.pdb')
    
                        ifgchid, ifgresnum = row['chid_ifg'], row['resnum_ifg']
                        vdmchid, vdmresnum = row['chid_vdm'], row['resnum_vdm']
                        ifgseg, vdmseg = row['segi_ifg'],row['segi_vdm']

                        printout = copy.deepcopy(par)
                        printout = printout.select(
                            '(chain {} and segment {} and resnum {}) or (chain {} and segment {} and resnum {})'.format(
                            ifgchid,ifgseg,ifgresnum,vdmchid,vdmseg,vdmresnum))
                        printout.select('chain {} and segment {} and resnum {}'.format(ifgchid,ifgseg,ifgresnum)).setChids('Y')
                        printout.select('chain {} and segment {} and resnum {}'.format(vdmchid,vdmseg,vdmresnum)).setChids('X')
                        printout.select('all').setResnums(10)
                        printout_interactamer = []
                        integrin_interactamer = []
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
                            interactamer_transf = pr.applyTransformation(each[1], printout_interactamer)
                            outdir = './output_data/pdbfiles/'
                            threecode = constants.AAname[ifgresn]
                            if pdbix != None and pdbname != None:
                                pr.writePDB(outdir+'{}_{}_{}_{}{}_{}{}_{}_{}_{}'.format(
                                    pdbix, pdbname, targetresi,ifginfo[1],
                                    ifgresn,vdminfo[1],vdmresn,cutoff,
                                    row.name, method), interactamer_transf)
                            else: 
                                pr.writePDB(outdir+'{}_{}{}_{}{}_{}_{}_{}'.format(
                                    targetresi,ifginfo[1],
                                    ifgresn,vdminfo[1],vdmresn,cutoff,
                                    row.name, method), interactamer_transf)
                            counter += 1
                        except:
                            traceback.print_exc()
                            pass

        # count how many NNs the query intrxn has
        num_nn, norm_metrics = get_NN(lookupatoms_to_clus, 
            num_atoms, rmsds, query, cutoff, num_all_vdms, num_direct)
        return ifginfo[0], ifginfo[1], ifgresn, vdminfo[0], vdminfo[1],\
            vdmresn, ifgatoms, vdmatoms, num_nn, norm_metrics

def get_NN(lookupatoms_to_clus, num_atoms, rmsds, query, cutoff, num_all_vdms, num_direct):
    flat = [x.reshape(-1,) for x in lookupatoms_to_clus]
    euclid_dist = cutoff * np.sqrt(num_atoms)
    nbrs = NearestNeighbors(
        n_neighbors=1,
        metric='euclidean',
        radius=euclid_dist)
    nnfit = nbrs.fit(flat)
    query = query.reshape(1,-1)
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
    norm_metrics.append(num_direct)
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
