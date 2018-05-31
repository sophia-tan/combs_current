from sklearn.neighbors import NearestNeighbors
import sys
import traceback
import copy
from combs.apps import *
from cluster_for_sophia import *
import prody as pr
import numpy as np
import pickle as pkl
import pandas as pd
from itertools import *
sys.path.append('/home/gpu/Sophia/combs/src/')


def get_interacting_atoms(parsed, resi, resn):
    target = parsed.select('chain %s and resnum %s' % (resi[0], resi[1:]))
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
        ifgatoms = ''
        vdmatoms = ''
        for typ in [bb, constants.planar_atoms[ifgresn]]:
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


def getcoords(row, info):
    ifgres, vdmres, ifgatoms, vdmatoms, query_atoms, lookup, targetresi, ignore = info
    try:
        try:
            db_dir = '/home/gpu/Sophia/STcombs/20171118/database/reduce/'
            par = pr.parsePDB(db_dir + row['pdb'] + 'H.pdb')
        except BaseException:
            db_dir = '/home/gpu/Sophia/combs/st_wd/20180207_db_molprobity_biolassem/'
            par = pr.parsePDB(db_dir + row['pdb'] + 'H.pdb')

        ifgchid, ifgresnum = row['chid_ifg'], row['resnum_ifg']
        vdmchid, vdmresnum = row['chid_vdm'], row['resnum_vdm']
        printout = copy.deepcopy(par)
        printout = printout.select(
            '(chain {} and resnum {}) or (chain {} and resnum {})'.format(
                ifgchid, ifgresnum, vdmchid, vdmresnum))
        printout.select(
            'chain {} and resnum {}'.format(
                ifgchid, ifgresnum)).setChids('Y')
        printout.select(
            'chain {} and resnum {}'.format(
                vdmchid, vdmresnum)).setChids('X')
        printout.select('all').setResnums(10)

        ifgatoms = list(set(ifgatoms))
        vdmatoms = list(set(vdmatoms))

        # account for flipped residues
        # 1) do ifg
        lookupifg_list = []
        lookupifg_list.append(ifgatoms)
        if ifgres in constants.flip_names:
            copyls = ifgatoms.copy()
            for ix, atom in enumerate(ifgatoms):
                for pair in constants.flip_names[ifgres]:
                    if pair[0] not in ifgatoms and pair[1] not in ifgatoms:
                        pass
                    elif pair[0] in ifgatoms or pair[1] in ifgatoms:
                        if atom in pair:
                            if atom == pair[0]:
                                copyls[ix] = pair[1]
                            elif atom == pair[1]:
                                copyls[ix] = pair[0]
                    if copyls not in lookupifg_list:
                        lookupifg_list.append(copyls)
        # 2) do vdm
        lookupvdm_list = []
        lookupvdm_list.append(vdmatoms)
        if vdmres in constants.flip_names:
            copyls = vdmatoms.copy()
            for ix, atom in enumerate(vdmatoms):
                for pair in constants.flip_names[vdmres]:
                    if pair[0] not in vdmatoms and pair[1] not in vdmatoms:
                        pass
                    elif pair[0] in vdmatoms or pair[1] in vdmatoms:
                        if atom in pair:
                            if atom == pair[0]:
                                copyls[ix] = pair[1]
                            elif atom == pair[1]:
                                copyls[ix] = pair[0]
                    if copyls not in lookupvdm_list:
                        lookupvdm_list.append(copyls)

        outdir = './output_data/pdbfiles/'
        pr.writePDB(
            outdir +
            '{}_{}_{}_{}.pdb'.format(
                targetresi,
                ifgres,
                vdmres,
                row.name),
            printout)

        def get_lookup_relevant_atoms(
                par,
                ifgatoms,
                ifgchid,
                ifgresnum,
                vdmatoms,
                vdmchid,
                vdmresnum):
            coords = []
            for atom in ifgatoms:
                selection = par.select(
                    'chain {} and resnum {} and name {}'.format(
                        ifgchid, ifgresnum, atom))
                coords.append(selection.getCoords()[0])
            for atom in vdmatoms:
                selection = par.select(
                    'chain {} and resnum {} and name {}'.format(
                        vdmchid, vdmresnum, atom))
                coords.append(selection.getCoords()[0])
            return np.array(coords)

        compare_rmsds = []
        query_atoms = np.array(query_atoms)
        for ifgls in lookupifg_list:
            for vdmls in lookupvdm_list:
                lookupatoms = get_lookup_relevant_atoms(
                    par, ifgls, ifgchid, ifgresnum, vdmls, vdmchid, vdmresnum)
                moved, transf = pr.superpose(lookupatoms, query_atoms)
                compare_rmsds.append(pr.calcRMSD(moved, query_atoms))
        return [min(compare_rmsds), len(query_atoms)]

    except Exception:
        traceback.print_exc()
        return np.nan

# def get_coordslistforcluster(targetresi,parsed,ifg, vdm, ifgatoms,vdmatoms,lookup_dir,db_dir,ifginfo,vdminfo,info,int_res):
#    try:
#        pkl.load(open('./output_data/{}_{}_{}_ifg_{}_vdm_{}_coords.pkl'.format(targetresi,int_res[0],int_res[1],ifg,vdm),'rb'))
#        return []
#    except:
#        query_atoms = []
#        ifgatoms = list(set(ifgatoms))
#        vdmatoms = list(set(vdmatoms))
#        for atom in ifgatoms:
#            selection = parsed.select('chain {} and resnum {} and name {}'.format(ifginfo[0],ifginfo[1],atom))
#            query_atoms.append(selection.getCoords()[0])
#        for atom in vdmatoms:
#            selection = parsed.select('chain {} and resnum {} and name {}'.format(vdminfo[0],vdminfo[1],atom))
#            query_atoms.append(selection.getCoords()[0])
#
#        # get coords for all combed vdms
#        lookup = pkl.load(open(lookup_dir+'refinedvdms/vdms_of_{}.pkl'.format(ifg), 'rb'))
#        lookup = lookup[lookup['resname_vdm']==constants.AAname_rev[vdm]]
#        #lookup = lookup[:5] ###delete###
#        info += [query_atoms, lookup, targetresi,int_res[0]]
#        lookup['rmsds'] = lookup.apply(getcoords,info=info,axis=1)
#        pkl.dump(lookup, open('./output_data/{}_{}_{}_ifg_{}_vdm_{}_coords.pkl'.format(targetresi,int_res[0],int_res[1],ifg,vdm),'wb'))
#        return lookup


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


def num_mems_query_clus(mems):
    '''find # of members in same cluster as query intrxn'''
    for clus in mems:
        for mem in clus:
            if mem == 0:
                return len(clus)


def score_interaction_and_dump(
        parsed,
        ifgresn,
        vdmresn,
        ifg_contact_atoms,
        vdm_contact_atoms,
        method,
        targetresi, cutoff):
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
        # lookupcoords = lookupcoords[:100] # delete

        ifglists = flip(ifgatoms, ifgresn)
        vdmlists = flip(vdmatoms, vdmresn)
        rmsds = []
        num_atoms = len(query)
        coords_ls = [
            item for item in lookupcoords if item[0] in lookupdf.index]
        lookupatoms_to_clus = []
        lookupatoms_to_clus.append(query)  # first element is always query

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
                for each in ifg_vdm_ind:
                    if each[1] == min(compare_rmsds):
                        lookupatoms_to_clus.append(moved)

            else:
                rmsds.append([int(item[0]), 100000])

        # count how many NNs the query intrxn has
        num_nn, norm_metrics = get_NN(lookupatoms_to_clus, 
            num_atoms, rmsds, query, cutoff, num_all_vdms)
        output_info = [num_nn, norm_metrics]
        pkl.dump(output_info, 
            open('./output_data/{}_{}{}_{}{}_matches_{}_{}.pkl'.format(
            targetresi,ifginfo[1],ifgresn,vdminfo[1],vdmresn,
            method,cutoff),'wb'))


def get_NN(lookupatoms_to_clus, num_atoms, rmsds, query, cutoff, num_all_vdms):
    flat = [x.reshape(-1,) for x in lookupatoms_to_clus]
    euclid_dist = cutoff * np.sqrt(num_atoms)
    print('euclidean distance: {} for # atoms: {} with rmsd: {}'
            .format(euclid_dist,num_atoms,cutoff))
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

#            try:
#            except:
#                np.savez_compressed('./output_data/{}_{}{}_{}{}_pairwisematrix_{}'.format(targetresi,\
#                    ifginfo[1],ifgresn,vdminfo[1],vdmresn,method),D)
#
#        adj_mat = make_adj_mat(D, 0.5)
#        mems,centroids = greedy(adj_mat)
#        avg_clus_size_wo_singles = np.mean([len(x) for x in mems if len(x) > 1])
#        avg_clus_size_with_singles = np.mean([len(x) for x in mems])
#        med_clus_size_wo_singles = np.median([len(x) for x in mems if len(x) > 1])
#        med_clus_size_with_singles = np.median([len(x) for x in mems])
#        num_clustered = len(lookupatoms_to_clus) # have direct intrxns w/ ifg
#
#        # find size of clus integrin is in (element 0 of coords array)
#        integrin_clus = num_mems_integrin_clus(mems)
#
#        rmsds.append([num_atoms,ifgatoms,vdmatoms,integrin_clus,num_clustered,\
#        avg_clus_size_wo_singles, avg_clus_size_with_singles, med_clus_size_wo_singles,\
#        med_clus_size_with_singles])
#        rmsds = np.array(rmsds)
#        pkl.dump(rmsds, open('./output_data/{}_{}{}_{}{}_rmsds_{}.pkl'.format(targetresi,\
#            ifginfo[1],ifgresn,vdminfo[1],vdmresn,method),'wb'))
#        return rmsds
#    else:
#        print('ERROR',ifginfo,ifgresn,vdminfo,vdmresn)
#