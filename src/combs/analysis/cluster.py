__all__ = ['cluster_rmsd_mat', 'cluster_adj_mat', 'make_ifg_sele', 'make_new_coord_sys', 'make_ifg_coords_array', 'make_bb_coords_array']

import numpy as np
import prody as pr
import traceback
from scipy.sparse import csr_matrix


def greedy_cluster(adj_mat):
    """Takes an adjacency matrix as input.
        All values of adj_mat are 1 or 0:  1 if <= to cutoff, 0 if > cutoff.
        Can generate adj_mat from data in column format with:
        sklearn.neighbors.NearestNeighbors(metric='euclidean', radius=cutoff).fit(data).radius_neighbors_graph(data)"""

    if not isinstance(adj_mat, csr_matrix):
        try:
            adj_mat = csr_matrix(adj_mat)
        except:
            print('adj_mat distance matrix must be scipy csr_matrix (or able to convert to one)')
            return

    assert adj_mat.shape[0] == adj_mat.shape[1], 'Distance matrix is not square.'

    all_mems = []
    cents = []
    indices = np.arange(adj_mat.shape[0])


    while adj_mat.shape[0] > 0:
        cent = adj_mat.sum(axis=1).argmax()
        cents.append(indices[cent])
        row = adj_mat.getrow(cent)
        tf = ~row.toarray().astype(bool)[0]
        mems = indices[~tf]
        all_mems.append(mems)
        indices = indices[tf]
        adj_mat = adj_mat[tf][:, tf]

    return all_mems, cents


def greedy_cluster_gt(adj_mat, gt=1):
    """Takes an adjacency matrix as input.
        All values of adj_mat are 1 or 0:  1 if <= to cutoff, 0 if > cutoff.
        Can generate adj_mat from data in column format with:
        sklearn.neighbors.NearestNeighbors(metric='euclidean', radius=cutoff).fit(data).radius_neighbors_graph(data)"""

    if not isinstance(adj_mat, csr_matrix):
        try:
            adj_mat = csr_matrix(adj_mat)
        except:
            print('adj_mat distance matrix must be scipy csr_matrix (or able to convert to one)')
            return

    assert adj_mat.shape[0] == adj_mat.shape[1], 'Distance matrix is not square.'

    all_mems = []
    cents = []
    indices = np.arange(adj_mat.shape[0])
    mems = np.ones(len(range(gt + 1)))

    while len(mems) > gt:
        cent = adj_mat.sum(axis=1).argmax()
        cents.append(indices[cent])
        row = adj_mat.getrow(cent)
        tf = ~row.toarray().astype(bool)[0]
        mems = indices[~tf]
        all_mems.append(mems)
        indices = indices[tf]
        adj_mat = adj_mat[tf][:, tf]

    return all_mems, cents


def greedy_cluster_pc(adj_mat, pc=0.5):
    """Takes an adjacency matrix as input.
        All values of adj_mat are 1 or 0:  1 if <= to cutoff, 0 if > cutoff.
        Can generate adj_mat from data in column format with:
        sklearn.neighbors.NearestNeighbors(metric='euclidean', radius=cutoff).fit(data).radius_neighbors_graph(data)"""

    if not isinstance(adj_mat, csr_matrix):
        try:
            adj_mat = csr_matrix(adj_mat)
        except:
            print('adj_mat distance matrix must be scipy csr_matrix (or able to convert to one)')
            return

    assert adj_mat.shape[0] == adj_mat.shape[1], 'Distance matrix is not square.'

    all_mems = []
    cents = []
    indices = np.arange(adj_mat.shape[0])
    totlen = adj_mat.shape[0]
    while adj_mat.shape[0] / totlen > pc:
        cent = adj_mat.sum(axis=1).argmax()
        cents.append(indices[cent])
        row = adj_mat.getrow(cent)
        tf = ~row.toarray().astype(bool)[0]
        mems = indices[~tf]
        all_mems.append(mems)
        indices = indices[tf]
        adj_mat = adj_mat[tf][:, tf]

    return all_mems, cents

def cluster_rmsd_mat(rmsd_matrix, cutoff=1, method='Greedy'):
    if method == 'Greedy':
        indices = np.arange(0, len(rmsd_matrix))
        boolean_rmsd_matrix = rmsd_matrix <= cutoff
        centroids = []
        members = []
        while boolean_rmsd_matrix.any():
            centroid = np.argmax(np.sum(boolean_rmsd_matrix, axis=0))
            mem_indices = boolean_rmsd_matrix[:, centroid]
            mems = indices[mem_indices]
            boolean_rmsd_matrix[mems, :] = False
            boolean_rmsd_matrix[:, mems] = False
            centroids.append(centroid)
            members.append(mems)
        return members, centroids


def cluster_adj_mat(adjacency_matrix):
    rows = adjacency_matrix.nonzero()[0]
    cols = adjacency_matrix.nonzero()[1]
    member_sets = []
    centroids = []
    while cols.any():
        un, inv, coun = np.unique(cols, return_counts=True, return_inverse=True)
        centroid = un[np.argmax(coun)]
        centroids.append(centroid)
        member_set = rows[cols == centroid]
        member_sets.append(member_set)
        ind_col = np.in1d(un, member_set, assume_unique=True, invert=True)
        ind_col = ind_col[inv]
        rows = rows[ind_col]
        cols = cols[ind_col]
        ind_row = np.in1d(rows, member_set, invert=True)
        cols = cols[ind_row]
        rows = rows[ind_row]
    return member_sets, centroids


def cluster_adj_mat_top_pc(adjacency_matrix, pc=0.5):
    rows = adjacency_matrix.nonzero()[0]
    cols = adjacency_matrix.nonzero()[1]
    member_sets = []
    centroids = []
    totlen = len(cols)
    while len(cols)/totlen > pc:
        un, inv, coun = np.unique(cols, return_counts=True, return_inverse=True)
        centroid = un[np.argmax(coun)]
        centroids.append(centroid)
        member_set = rows[cols == centroid]
        member_sets.append(member_set)
        ind_col = np.in1d(un, member_set, assume_unique=True, invert=True)
        ind_col = ind_col[inv]
        rows = rows[ind_col]
        cols = cols[ind_col]
        ind_row = np.in1d(rows, member_set, invert=True)
        cols = cols[ind_row]
        rows = rows[ind_row]
    return member_sets, centroids


def cluster_adj_mat_gt(adjacency_matrix, gt=1):
    rows = adjacency_matrix.nonzero()[0]
    cols = adjacency_matrix.nonzero()[1]
    member_sets = []
    centroids = []
    member_set = np.ones(len(range(gt + 1)))
    while len(member_set) > gt:
        un, inv, coun = np.unique(cols, return_counts=True, return_inverse=True)
        centroid = un[np.argmax(coun)]
        centroids.append(centroid)
        member_set = rows[cols == centroid]
        member_sets.append(member_set)
        ind_col = np.in1d(un, member_set, assume_unique=True, invert=True)
        ind_col = ind_col[inv]
        rows = rows[ind_col]
        cols = cols[ind_col]
        ind_row = np.in1d(rows, member_set, invert=True)
        cols = cols[ind_row]
        rows = rows[ind_row]
    return member_sets, centroids


def make_ifg_sele(pdb, comb):
    """uses iFG definitions in comb object tmeo select iFGs in the parsed protein object that have all atoms
    and occupancies = 1.
    """
    # There is a problem with this code: What if one wants to select atoms from a HEME as an iFG?
    # It is not represented in the one_letter_code dictionary...
    possible_ifgs = []
    if comb.num_res_ifg == 1:
        poss_ifg_sel = pdb.select('chain Y and resnum 10')
        if poss_ifg_sel is not None:
            ifg_resindices, indices = np.unique(poss_ifg_sel.getResindices(), return_index=True)
            ifg_resnames = poss_ifg_sel.getResnames()[indices]
            for ifg_resindex, ifg_resname in zip(ifg_resindices, ifg_resnames):
                    ifg_selection = pdb.select('resindex ' + str(ifg_resindex) + ' and name '
                                                          + comb.ifg_sele_dict[1][ifg_resname])
                    if ifg_selection is not None:
                        num_atoms = len(ifg_selection)
                        if num_atoms == len(comb.ifg_sele_dict[1][ifg_resname].split()):
                            if all(ifg_selection.getResnums() > 0):
                                possible_ifgs.append(ifg_selection)
        return possible_ifgs[0]
    else:
        poss_ifg_sel = pdb.select('chain Y and resnum 10to11')
        if poss_ifg_sel is not None:
            ifg_resindices_cat_list, indices = np.unique(poss_ifg_sel.getResindices(), return_index=True)
            ifg_resnames_cat_list = poss_ifg_sel.getResnames()[indices]
            ifg_resindex_pairs = [ifg_resindices_cat_list[i:i + 2] for i in range(0, len(ifg_resindices_cat_list), 2)]
            ifg_resname_pairs = [ifg_resnames_cat_list[i:i + 2] for i in range(0, len(ifg_resnames_cat_list), 2)]
            for ifg_resindex_pair, ifg_resname_pair in zip(ifg_resindex_pairs, ifg_resname_pairs):
                resind1, resind2 = ifg_resindex_pair
                resname1, resname2 = ifg_resname_pair
                try:
                    ifg_selection = pdb.select('(resindex ' + str(resind1) + ' and name '
                                                      + comb.ifg_sele_dict[1][resname1]+')'
                                                      + ' or (resindex ' + str(resind2) + ' and name '
                                                      + comb.ifg_sele_dict[2][resname2]+')')
                except KeyError:
                    print('Non-canonical residue in iFG, skipping.')
                    ifg_selection = None
                if ifg_selection is not None:
                    num_atoms = len(ifg_selection)
                    names = comb.ifg_sele_dict[1][resname1].split()
                    names.extend(comb.ifg_sele_dict[2][resname2].split())
                    if num_atoms == len(names):
                        if all(ifg_selection.getResnums() > 0):
                            possible_ifgs.append(ifg_selection)
        return possible_ifgs[0]


def make_new_coord_sys(pdb, comb, origin_atom, plane_atom1, plane_atom2):
    origin_coords = pdb.select('chain X and resnum 10 and name ' + origin_atom).getCoords()[0]
    pdb_coords = pdb.getCoords()
    pdb_coords_neworigin = pdb_coords - origin_coords
    plane_atom1_coords = pdb.select('chain X and resnum 10 and name ' + plane_atom1).getCoords()[0] - origin_coords
    plane_atom2_coords = pdb.select('chain X and resnum 10 and name ' + plane_atom2).getCoords()[0] - origin_coords
    x_norm = plane_atom1_coords / np.linalg.norm(plane_atom1_coords)
    orthvec = np.cross(plane_atom1_coords, plane_atom2_coords)
    z_norm = -1 * orthvec / np.linalg.norm(orthvec)
    orthvec2 = np.cross(plane_atom1_coords, orthvec)
    y_norm = orthvec2 / np.linalg.norm(orthvec2)
    R = np.array([x_norm, y_norm, z_norm])
    pdb_coords_neworigin_rot = np.dot(pdb_coords_neworigin, R.T)
    pdbcopy = pdb.copy()
    pdbcopy.setCoords(pdb_coords_neworigin_rot)
    ifg_sel = make_ifg_sele(pdbcopy, comb)
    ifg_coords = [y for x in ifg_sel.getCoords() for y in x]
    return ifg_coords, pdbcopy

def make_new_coord_sys_gen(pdb, comb, origin_atom, plane_atom1, plane_atom2, ifg_selection):
    origin_coords = pdb.select(origin_atom).getCoords()[0]
    pdb_coords = pdb.getCoords()
    pdb_coords_neworigin = pdb_coords - origin_coords
    plane_atom1_coords = pdb.select(plane_atom1).getCoords()[0] - origin_coords
    plane_atom2_coords = pdb.select(plane_atom2).getCoords()[0] - origin_coords
    x_norm = plane_atom1_coords / np.linalg.norm(plane_atom1_coords)
    orthvec = np.cross(plane_atom1_coords, plane_atom2_coords)
    z_norm = -1 * orthvec / np.linalg.norm(orthvec)
    orthvec2 = np.cross(plane_atom1_coords, orthvec)
    y_norm = orthvec2 / np.linalg.norm(orthvec2)
    R = np.array([x_norm, y_norm, z_norm])
    pdb_coords_neworigin_rot = np.dot(pdb_coords_neworigin, R.T)
    pdbcopy = pdb.copy()
    pdbcopy.setCoords(pdb_coords_neworigin_rot)
    ifg_sel = make_ifg_sele_gen(pdbcopy, ifg_selection, comb)
    ifg_coords = [y for x in ifg_sel.getCoords() for y in x]
    return ifg_coords, pdbcopy

def make_ifg_sele_gen(pdb, ifg_selection, comb):
    """uses iFG definitions in comb object to select iFGs in the parsed protein object that have all atoms
    and occupancies = 1.
    """
    # There is a problem with this code: What if one wants to select atoms from a HEME as an iFG?
    # It is not represented in the one_letter_code dictionary...
    possible_ifgs = []
    if comb.num_res_ifg == 1:
        poss_ifg_sel = pdb.select(ifg_selection)
        if poss_ifg_sel is not None:
            ifg_resindices, indices = np.unique(poss_ifg_sel.getResindices(), return_index=True)
            ifg_resnames = poss_ifg_sel.getResnames()[indices]
            for ifg_resindex, ifg_resname in zip(ifg_resindices, ifg_resnames):
                    ifg_selection = pdb.select('resindex ' + str(ifg_resindex) + ' and name '
                                                          + comb.ifg_sele_dict[1][ifg_resname])
                    if ifg_selection is not None:
                        num_atoms = len(ifg_selection)
                        if num_atoms == len(comb.ifg_sele_dict[1][ifg_resname].split()):
                            if all(ifg_selection.getResnums() > 0):
                                possible_ifgs.append(ifg_selection)
        return possible_ifgs[0]
    else:
        poss_ifg_sel = pdb.select(ifg_selection)
        if poss_ifg_sel is not None:
            ifg_resindices_cat_list, indices = np.unique(poss_ifg_sel.getResindices(), return_index=True)
            ifg_resnames_cat_list = poss_ifg_sel.getResnames()[indices]
            ifg_resindex_pairs = [ifg_resindices_cat_list[i:i + 2] for i in range(0, len(ifg_resindices_cat_list), 2)]
            ifg_resname_pairs = [ifg_resnames_cat_list[i:i + 2] for i in range(0, len(ifg_resnames_cat_list), 2)]
            for ifg_resindex_pair, ifg_resname_pair in zip(ifg_resindex_pairs, ifg_resname_pairs):
                resind1, resind2 = ifg_resindex_pair
                resname1, resname2 = ifg_resname_pair
                try:
                    ifg_selection = pdb.select('(resindex ' + str(resind1) + ' and name '
                                                      + comb.ifg_sele_dict[1][resname1]+')'
                                                      + ' or (resindex ' + str(resind2) + ' and name '
                                                      + comb.ifg_sele_dict[2][resname2]+')')
                except KeyError:
                    print('Non-canonical residue in iFG, skipping.')
                    ifg_selection = None
                if ifg_selection is not None:
                    num_atoms = len(ifg_selection)
                    names = comb.ifg_sele_dict[1][resname1].split()
                    names.extend(comb.ifg_sele_dict[2][resname2].split())
                    if num_atoms == len(names):
                        if all(ifg_selection.getResnums() > 0):
                            possible_ifgs.append(ifg_selection)
        return possible_ifgs[0]

def make_ifg_coords_array(pdbs, comb, origin_atom, plane_atom1, plane_atom2):
    if isinstance(pdbs, pr.atomic.atomgroup.AtomGroup):
        pdbs = [pdbs]
    ifg_coords_matrix = []
    pdbs_rot = []
    for pdb in pdbs:
        try:
            ifg_coords, pdb_rot = make_new_coord_sys(pdb, comb, origin_atom, plane_atom1, plane_atom2)
            ifg_coords_matrix.append(ifg_coords)
            pdbs_rot.append(pdb_rot)
        except Exception:
            traceback.print_exc()
    return np.array(ifg_coords_matrix), pdbs_rot

def superpose(mob_obj, mob_sel, targ_sel, flips=None):
    mob_names = mob_sel.getNames()
    targ_names = targ_sel.getNames()
    mob_sel = mob_sel.getCoords()
    targ_sel = targ_sel.getCoords()
    if any(mob_names != targ_names):
        mob_ind = np.argsort(mob_names)
        targ_ind = np.argsort(targ_names)
        mob_names = mob_names[mob_ind]
        targ_names = targ_names[targ_ind]
        mob_sel = mob_sel[mob_ind]
        targ_sel = targ_sel[targ_ind]
    if flips:
        for names in flips:
            i, = np.where(mob_names == names[0])
            j, = np.where(mob_names == names[1])
            temp = mob_sel[i]
            mob_sel[i] = mob_sel[j]
            mob_sel[j] = temp
    mob_sel_com = mob_sel.mean(0)
    targ_sel_com = targ_sel.mean(0)
    mob_sel_cen = mob_sel - mob_sel_com
    targ_sel_cen = targ_sel - targ_sel_com
    cov_matrix = np.dot(mob_sel_cen.T, targ_sel_cen)
    U, S, Wt = np.linalg.svd(cov_matrix)
    R = np.dot(U, Wt)
    if np.linalg.det(R) < 0.:
        Wt[-1] *= -1
        R = np.dot(U, Wt)
    mob_obj_copy = mob_obj.copy()
    mob_obj = mob_obj.getCoords()
    mob_obj_cen_rot = np.dot((mob_obj-mob_sel_com), R)
    mob_obj_cen_rot_targ = mob_obj_cen_rot + targ_sel_com
    mob_obj_copy.setCoords(mob_obj_cen_rot_targ)
    return mob_obj_copy

def make_new_coord_sys_ifg_align(pdb, targ_sel, comb):
    ifg_sel = make_ifg_sele(pdb, comb)
    pdb_copy_super = superpose(pdb, ifg_sel, targ_sel)
    bb_coords = [y for x in pdb_copy_super.select('chain X and name CA').getCoords() for y in x]
    return bb_coords, pdb_copy_super

def make_bb_coords_array(pdbs, comb):
    assert isinstance(pdbs, list)
    bb_coords_matrix = []
    pdbs_rot = []
    align_to_pdb = pdbs[0]
    for pdb in pdbs:
        try:
            bb_coords, pdb_rot = make_new_coord_sys_ifg_align(pdb, make_ifg_sele(align_to_pdb, comb), comb)
            bb_coords_matrix.append(bb_coords)
            pdbs_rot.append(pdb_rot)
        except Exception:
            traceback.print_exc()
    return np.array(bb_coords_matrix), pdbs_rot

# pair_dist = scipy.spatial.distance.pdist(iFGcoords_array)
#
# pair_dist_matrix = scipy.spatial.distance.squareform(pair_dist)
#
# ifg_dict = {'ASN': 'CB CG ND2 OD1', 'GLN': 'CG CD NE2 OE1'}
#    ...: cb = combs.Comb(ifg_dict)
#    ...:
#    ...: ifg_coords_asn,asn_vdms =combs.analysis.make_ifg_coords_array(vdms['ASN'], cb, 'CG', 'ND2', 'OD
#    ...: 1')
#    ...: from scipy import spatial
#    ...: pair_dist = spatial.distance.pdist(ifg_coords_asn)
#    ...: pair_dist_mat = spatial.distance.squareform(pair_dist)
#    ...: mem,cen = combs.cluster(pair_dist_mat)
#    ...: mem[0]

# import os
# from prody import *
# for i,clust in enumerate(mem[:10]):
#     ...:     os.mkdir('/Users/npolizzi/Desktop/clusters_asp/' + str(i))
#     ...:     for m in clust:
#     ...:         writePDB('/Users/npolizzi/Desktop/clusters_asp/' + str(i) + '/' + repr(asp_vdms[m]).sp
#     ...: lit()[1] + '.pdb.gz', asp_vdms[m])
