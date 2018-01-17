__all__ = ['get_rois', 'make_relative_interactamer_coords', 'find_hotspot', 'print_hotspot_pdbs', 'make_interactamers',
           'load_interactamers', 'get_rois_phi_psi', 'make_rel_vdms', 'load_rel_vdms', 'find_hotspot_bb',
           'print_hotspot_pdbs_bb']

# get residues of interest (roi) in protein of interest (poi)
# do translation and rotation of vdms onto rois in poi
# make a graph of ifgs which are grouped by roi:
#   ifg location is a node. edge between nodes if distance between ifgs is < some cutoff and ifgs are from different roi
# find longest path in graph
from scipy import spatial
import networkx as nx
import itertools
import numpy as np
import traceback
import os
import pickle
import prody as pr
import pandas as pd
from ..apps.constants import interactamer_atoms, flip_residues, flip_sets, flip_names
from ..analysis.Cluster import superpose
from ..analysis.Cluster import make_ifg_sele
from ..analysis.Analysis import Analyze
from sklearn.neighbors import NearestNeighbors
import functools
# think about using multipool

#make interactamers and print to file:
#fmt should be: ifg_count, vmd_count, vdm_interact_coords, ifg_interact_coords
# for each iFG of interest:
    # 1. load the interactamers
    # 2. find transformation of vdm_interact_coords onto roi
    # 3. move ifg_interact_coords via that transformation.
    # 4. remove clashing ifg_interact_coords.  Need to also keep track of ifg_count and vdm_count numbers of remaining interactamers.
    # 5. find groups of similar location ifg_interact_coords
# 6. find groups with distances between ifg_interact_coords that satisfy those on the ligand of interest (loi)
# 7. find poses of loi that superimpose onto an ifg from each group and loi does not clash with protein of interst (poi)

# sample = Sample()
# sample.poi=
# sample.bs_residues=
# sample.rois=
# sample.ligand_conformers=
# sample.load_rel_vdms(ifg_info)

# ifg_info1=IntFG_info()
# ifg_info1.name = 'carboxamide'
# ifg_info1.ifg_dict = ifg_dict_carboxamide
# ifg_info1.rel_vdm_path = rel_vdm_path
# ifg_info1.rel_vdms = {}
# ifg_info1.hotspots = []
#
# ifg_info2=IntFG_info()
# ifg_info2.name = 'amino'
# ifg_info2.ifg_dict = ifg_dict_amino
# ifg_info2.rel_vdm_path = rel_vdm_path
# ifg_info2.rel_vdms = {}
# ifg_info2.hotspots = []
#
# ifg_info3=IntFG_info()
# ifg_info3.name = 'carboxylate'
# ifg_info3.ifg_dict = ifg_dict_carboxylate
# ifg_info3.rel_vdm_path = rel_vdm_path
# ifg_info3.rel_vdms = {}
# ifg_info3.hotspots = []
#
# ifg_infos = [ifg_info1, ifg_info2, ifg_info3]
# for ifg_info in ifg_infos:
#     sample.load_rel_vdms(ifg_info) #delete clashing ones
#     ifg_info.find_hotspots(rmsd_cutoff=) #find and store hotspots
#     ifg_info.print_hotspots(outdir)  #optionally print hotspots
# Poses = []
# for conformer in sample.ligand_conformers:
#     sample.find_ligand_cst_hotspots(ifg_infos) #sets hotspot_sets attribute
#     for hotspot_set in sample.hotspot_sets:
#         pose = Pose(conformer, hotspot_set)
#         pose.dock_ligand()
#         pose.check_bb_clash()
#         pose.get_sequences()
#         Poses.append(pose)
# sample.store_poses(Poses)  #store poses in a pickle

#Pose class should have a print_pose() method.





def get_rotation_translation(mob_coords, targ_coords):
    mob_coords_com = mob_coords.mean(0)
    targ_coords_com = targ_coords.mean(0)
    mob_coords_cen = mob_coords - mob_coords_com
    targ_coords_cen = targ_coords - targ_coords_com
    cov_matrix = np.dot(mob_coords_cen.T, targ_coords_cen)
    U, S, Wt = np.linalg.svd(cov_matrix)
    R = np.dot(U, Wt)
    if np.linalg.det(R) < 0.:
        Wt[-1] *= -1
        R = np.dot(U, Wt)
    return R, mob_coords_com, targ_coords_com


def move_interactamers(apply_coords, rotation, int_com, targ_com):
    apply_coords_cen_rot = np.dot((apply_coords-int_com), rotation)
    return apply_coords_cen_rot + targ_com


def get_rois(poi, resn_chid_pairs):
    return [poi.select('resnum ' + str(resn) + ' and chain ' + chid) for resn, chid in resn_chid_pairs]


def get_rois_phi_psi(poi, resn_chid_pairs):
    phi_psi = []
    for resn, chid in resn_chid_pairs:
        try:
            phi = pr.calcPhi(poi[chid, resn])
        except:
            phi = None
            # traceback.print_exc()
        try:
            psi = pr.calcPsi(poi[chid, resn])
        except:
            psi = None
        phi_psi.append((phi, psi))
    return phi_psi


def get_poi_coords_for_clashsc(roi, poi, include_sc=False):
    if include_sc:
        coords_sel = poi.select('(not water and not element H D) exwithin 10 of sel',
                                   sel=roi.select('not element H D'))
    else:
        coords_sel = poi.select('(backbone and not element H D) exwithin 10 of sel',
                                    sel=roi.select('not element H D'))
    poi_coords = coords_sel.getCoords()
    return poi_coords


def calc_clash_score(poi_coords, ifg_coords, clash_cutoff=2.5):
    shape = ifg_coords.shape
    return np.sum(np.sum(spatial.distance.cdist(poi_coords, ifg_coords.flatten().reshape(shape[0]*shape[1],
                                                 shape[2])).T.reshape(shape[0], shape[1], len(poi_coords)) < clash_cutoff,
                         axis=1), axis=1)


def load_rel_vdms(path_to_pickle):
    interactamers = {}
    # interactamers['BB_CO'] = {}
    # interactamers['BB_NH'] = {}
    for file in [file for file in os.listdir(path_to_pickle) if file[0] is not '.']:
        splitfile = file.split('_')
        if splitfile[0] == 'BB':
            if splitfile[2] == 'N':
                with open(path_to_pickle + '/' + file, 'rb') as infile:
                    interactamers[splitfile[0] + '_'
                                  + splitfile[1]]['N_not_O'] = pickle.load(infile)
            elif splitfile[2] == 'O':
                with open(path_to_pickle + '/' + file, 'rb') as infile:
                    interactamers[splitfile[0] + '_'
                                  + splitfile[1]]['O_not_N'] = pickle.load(infile)
            else:
                with open(path_to_pickle + '/' + file, 'rb') as infile:
                    interactamers[splitfile[0] + '_'
                                  + splitfile[1]][tuple([int(n) for n in splitfile[2:-1]])] = pickle.load(infile)
        else:
            with open(path_to_pickle + '/' + file, 'rb') as infile:
                interactamers[file.split('rel_vdms')[0][:-1]] = pickle.load(infile)
    return interactamers


def get_roi_coords(roi, resn=None):
    resn = resn or roi.getResnames()[0]
    return np.array([roi.select('name ' + interactamer_atoms[resn][i]).getCoords()[0] for i in range(3)])


def find_hotspot(rois, phi_psis, interactamers, poi, rmsd_cutoff=1, include_sc_clash=False, clash_cutoff=2.5):
    """rois should be a list of prody selections of the residues of interest in the binding site. vdms should be in
    a dictionary format with keys are Resname and values as list of van der Mers of that resname"""

    bs_num_resn = {}
    for i, roi in enumerate(rois):
        bs_num_resn[i] = roi.getResnames()[0]

    bs_num_phi_psi = {}
    for i, phi_psi in enumerate(phi_psis):
        bs_num_phi_psi[i] = phi_psi
    # print(bs_num_phi_psi)

    bs_num_ifg_coords = {}
    bs_num_vdm_tags = {}
    for key in bs_num_resn.keys():
        poi_coords = get_poi_coords_for_clashsc(rois[key], poi, include_sc=include_sc_clash)

        # for sidechain
        roi_coords = get_roi_coords(rois[key])
        r, int_coords_com, roi_coords_com = get_rotation_translation(interactamers[bs_num_resn[key]][0, 4], roi_coords)
        ifg_coords = np.array([coords for coords in interactamers[bs_num_resn[key]][:, 5]])
        ifg_coords = move_interactamers(ifg_coords, r, int_coords_com, roi_coords_com)
        filter_index = np.where(calc_clash_score(poi_coords, ifg_coords, clash_cutoff=clash_cutoff) < 1)[0]
        if filter_index.size > 0:
            ifg_coords = ifg_coords[filter_index]
            # print(interactamers[bs_num_resn[key]][filter_index, :4])
            # print(np.array([bs_num_resn[key]]*len(filter_index)).T)
            # print(interactamers[bs_num_resn[key]][filter_index, :4].shape)
            # print(np.array([bs_num_resn[key]] * len(filter_index)).T.shape)
            # bs_num_vdm_tags[key] = np.concatenate((interactamers[bs_num_resn[key]][filter_index, :4],
            #                                        np.array([[bs_num_resn[key]]*len(filter_index)]).T), axis=1)
            # a[row_idx[:, None], col_idx]
            bs_num_vdm_tags[key] = interactamers[bs_num_resn[key]][filter_index[:, None], np.array([0, 1, 2, 3, 6])]

            bs_num_ifg_coords[key] = [coords.flatten() for coords in ifg_coords]
            # print(key, bs_num_vdm_tags[key].shape)
        # for backbone CO
        # need to load backbone interactamers in a special phipsi way
        for phi_low, phi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
            for psi_low, psi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
                if (phi_low <= bs_num_phi_psi[key][0] < phi_high) and (psi_low <= bs_num_phi_psi[key][1] < psi_high):
                    phi_psi_key = (phi_low, phi_high, psi_low, psi_high)

        # print(phi_psi_key)
        roi_coords = get_roi_coords(rois[key], resn='BB_CO')
        # print(roi_coords)
        # print(interactamers['BB_CO'][phi_psi_key][:,:4])
        try:
            r, int_coords_com, roi_coords_com = get_rotation_translation(interactamers['BB_CO'][phi_psi_key][0, 4], roi_coords)
            ifg_coords = np.array([coords for coords in interactamers['BB_CO'][phi_psi_key][:, 5]])
        # print(ifg_coords)
            ifg_coords = move_interactamers(ifg_coords, r, int_coords_com, roi_coords_com)
            filter_index = np.where(calc_clash_score(poi_coords, ifg_coords, clash_cutoff=clash_cutoff) < 1)[0]
        # print(ifg_coords)
        # print(filter_index)
        # print(bs_num_vdm_tags[key])
            if filter_index.size > 0:
                # print(filter_index[0].size)
                ifg_coords = ifg_coords[filter_index]

                bs_num_vdm_tags[key] = np.concatenate((bs_num_vdm_tags[key],
                                                       interactamers['BB_CO'][phi_psi_key][filter_index[:, None],
                                                                                           np.array([0, 1, 2, 3, 6])]), axis=0)
                                                # np.array([['BB_CO_' + '_'.join(str(nn) for nn in phi_psi_key)]*len(filter_index)]).T), axis=1)), axis=0)
                bs_num_ifg_coords[key].extend(coords.flatten() for coords in ifg_coords)
        except Exception:
            print('Number of interactamers in BB_CO phi/psi set ' + ','.join(str(x) for x in phi_psi_key) + ': ',
                  len(interactamers['BB_CO'][phi_psi_key]))
            traceback.print_exc()
            r, int_coords_com, roi_coords_com = get_rotation_translation(interactamers['BB_CO']['O_not_N'][0, 4],
                                                                         roi_coords)

        ifg_coords = np.array([coords for coords in interactamers['BB_CO']['O_not_N'][:, 5]])
        # print(ifg_coords)
        ifg_coords = move_interactamers(ifg_coords, r, int_coords_com, roi_coords_com)
        filter_index = np.where(calc_clash_score(poi_coords, ifg_coords, clash_cutoff=clash_cutoff) < 1)[0]
        # print(ifg_coords)
        # print(filter_index)
        # print(bs_num_vdm_tags[key])
        if filter_index.size > 0:
            # print(filter_index[0].size)
            ifg_coords = ifg_coords[filter_index]

            bs_num_vdm_tags[key] = np.concatenate((bs_num_vdm_tags[key],
                                                  interactamers['BB_CO']['O_not_N'][filter_index[:, None],
                                                                                    np.array([0, 1, 2, 3, 6])]), axis=0)
            bs_num_ifg_coords[key].extend(coords.flatten() for coords in ifg_coords)

            # print(key, bs_num_vdm_tags[key].shape)
        # for backbone NH
        roi_coords = get_roi_coords(rois[key], resn='BB_NH')
        try:
            r, int_coords_com, roi_coords_com = get_rotation_translation(interactamers['BB_NH'][phi_psi_key][0, 4],
                                                                     roi_coords)
            ifg_coords = np.array([coords for coords in interactamers['BB_NH'][phi_psi_key][:, 5]])
            ifg_coords = move_interactamers(ifg_coords, r, int_coords_com, roi_coords_com)
            filter_index = np.where(calc_clash_score(poi_coords, ifg_coords, clash_cutoff=clash_cutoff) < 1)[0]
            if filter_index.size > 0:
                ifg_coords = ifg_coords[filter_index]
                bs_num_vdm_tags[key] = np.concatenate((bs_num_vdm_tags[key],
                                                       interactamers['BB_NH'][phi_psi_key][filter_index[:, None],
                                                                                           np.array([0, 1, 2, 3, 6])]),
                                                      axis=0)
                                                # np.array([['BB_NH_' + '_'.join(str(nn) for nn in phi_psi_key)] * len(filter_index)]).T), axis=1)),axis=0)
                bs_num_ifg_coords[key].extend(coords.flatten() for coords in ifg_coords)
                # print(key, bs_num_vdm_tags[key].shape)
        except Exception:
            print('Number of interactamers in BB_NH phi/psi set ' + ','.join(str(x) for x in phi_psi_key) + ': ',
                  len(interactamers['BB_NH'][phi_psi_key]))
            traceback.print_exc()
            r, int_coords_com, roi_coords_com = get_rotation_translation(interactamers['BB_NH']['N_not_O'][0, 4],
                                                                         roi_coords)


        ifg_coords = np.array([coords for coords in interactamers['BB_NH']['N_not_O'][:, 5]])
        ifg_coords = move_interactamers(ifg_coords, r, int_coords_com, roi_coords_com)
        filter_index = np.where(calc_clash_score(poi_coords, ifg_coords, clash_cutoff=clash_cutoff) < 1)[0]
        if filter_index.size > 0:
            ifg_coords = ifg_coords[filter_index]
            bs_num_vdm_tags[key] = np.concatenate(
                (bs_num_vdm_tags[key], interactamers['BB_NH']['N_not_O'][filter_index[:, None],
                                                                         np.array([0, 1, 2, 3, 6])]), axis=0)
            # np.array([['BB_NH_' + '_'.join(str(nn) for nn in phi_psi_key)] * len(filter_index)]).T), axis=1)),axis=0)
            bs_num_ifg_coords[key].extend(coords.flatten() for coords in ifg_coords)

    g = nx.Graph()
    for i, j in itertools.combinations(bs_num_resn.keys(), 2):
        # print(i, len(bs_num_ifg_coords[i]))
        # print(j, len(bs_num_ifg_coords[j]))
        # print(i, bs_num_vdm_tags[i].shape)
        # print(j, bs_num_vdm_tags[j].shape)
        try:
            distpairs = spatial.distance.cdist(bs_num_ifg_coords[i], bs_num_ifg_coords[j])
            for m, n in np.array(np.where(distpairs < rmsd_cutoff)).T:
                g.add_nodes_from([(i, m), (j, n)])
                g.add_edge((i, m), (j, n))
        except Exception:
            traceback.print_exc()
    return g, nx.connected_components(g), bs_num_vdm_tags



def find_hotspot_bb(rois, interactamers, poi, rmsd_cutoff=1, include_sc_clash=False, clash_cutoff=2.5):
    """rois should be a list of prody selections of the residues of interest in the binding site. vdms should be in
    a dictionary format with keys are Resname and values as list of van der Mers of that resname"""

    # bs_num_resn = {}
    # for i, roi in enumerate(rois):
    #     bs_num_resn[i] = roi.getResnames()[0]

    # bs_num_phi_psi = {}
    # for i, phi_psi in enumerate(phi_psis):
    #     bs_num_phi_psi[i] = phi_psi
    # # print(bs_num_phi_psi)
    # print([val for val in interactamers.values()])
    interactamers = np.vstack((interactamers.values()))

    bs_num_ifg_coords = {}
    bs_num_ifg_coords_noflat = {}
    bs_num_vdm_tags = {}
    # for key in bs_num_resn.keys():
    for k in range(len(rois)):
        poi_coords = get_poi_coords_for_clashsc(rois[k], poi, include_sc=include_sc_clash)

        # for sidechain
        roi_coords = get_roi_coords(rois[k], resn='BB_NH')
        r, int_coords_com, roi_coords_com = get_rotation_translation(interactamers[0, 4], roi_coords)
        ifg_coords = np.array([coords for coords in interactamers[:, 5]])
        ifg_coords = move_interactamers(ifg_coords, r, int_coords_com, roi_coords_com)
        print('removing clash at pos', k)
        filter_index = np.where(calc_clash_score(poi_coords, ifg_coords, clash_cutoff=clash_cutoff) < 1)[0]
        if filter_index.size > 0:
            ifg_coords = ifg_coords[filter_index]
            # print(interactamers[bs_num_resn[key]][filter_index, :4])
            # print(np.array([bs_num_resn[key]]*len(filter_index)).T)
            # print(interactamers[bs_num_resn[key]][filter_index, :4].shape)
            # print(np.array([bs_num_resn[key]] * len(filter_index)).T.shape)
            # bs_num_vdm_tags[key] = np.concatenate((interactamers[bs_num_resn[key]][filter_index, :4],
            #                                        np.array([[bs_num_resn[key]]*len(filter_index)]).T), axis=1)
            # a[row_idx[:, None], col_idx]
            bs_num_vdm_tags[k] = interactamers[filter_index[:, None], np.array([0, 1, 2, 3, 6])]

            bs_num_ifg_coords[k] = [coords.flatten() for coords in ifg_coords]

            bs_num_ifg_coords_noflat[k] = ifg_coords
            # print(key, bs_num_vdm_tags[key].shape)

    g = nx.Graph()
    for i, j in itertools.combinations(list(range(len(rois))), 2):
        # print(i, len(bs_num_ifg_coords[i]))
        # print(j, len(bs_num_ifg_coords[j]))
        # print(i, bs_num_vdm_tags[i].shape)
        # print(j, bs_num_vdm_tags[j].shape)
        print('graph ', i, j)
        try:
            # distpairs = spatial.distance.cdist(bs_num_ifg_coords[i], bs_num_ifg_coords[j])
            nbrs = NearestNeighbors(n_neighbors=300, metric='euclidean', radius=rmsd_cutoff).fit(bs_num_ifg_coords[j])
            distances, indices = nbrs.kneighbors(bs_num_ifg_coords[i])
            where = np.where(distances < rmsd_cutoff)
            for m, p in zip(where[0], where[1]):
                n = indices[m, p]
            # for m, n in indices: #np.array(np.where(distpairs < rmsd_cutoff)).T:
                g.add_nodes_from([(i, m), (j, n)])
                g.add_edge((i, m), (j, n))
        except Exception:
            traceback.print_exc()
    return g, nx.connected_components(g), bs_num_vdm_tags, bs_num_ifg_coords_noflat


# def compare_hotspots(all_hotspots, lig_cst):
#     for ifg1_hotspot in all_hotspots[0]:
#         for hotspots in all_hotspots[1:]:
#             for hotspot in hotspots:
#                 distpairs = spatial.distance.cdist(ifg1_hotspot, hotspot)
#                 if any(lig_cst_1 - 1 < distpairs < lig_cst_1 + 1)

# def compare_hotspots(ifg1_hotspots, ifg2_hotspots, ifg1_coords, ifg2_coords, lig_cst, tol, ifg1name, ifg2name):
#     ifg1_coords = np.array([ifg1_coords[spot[0]][spot[1]] for hotspot in ifg1_hotspots for spot in hotspot])
#     ifg2_coords = np.array([ifg2_coords[spot[0]][spot[1]] for hotspot in ifg2_hotspots for spot in hotspot])
#     dist = spatial.distance.cdist(ifg1_coords.reshape(-1, 3), ifg2_coords.reshape(-1, 3))
#     dist_split = np.hsplit(dist, len(ifg2_coords))
#     ifg1_num = len(ifg1_coords)
#     dist = np.array([d.reshape(ifg1_num, int(np.prod(d.shape) / ifg1_num)) for d in dist_split])
#     lig_cst = lig_cst.ravel()
#     test = (lig_cst - tol < dist) & (dist < lig_cst + tol)
#     wh = np.where(np.all(test.T, axis=0).T)
#     indices = list(zip(wh[0], wh[1]))
#     return pd.DataFrame(indices,
#                         columns=[ifg2name, ifg1name]), ifg1_coords, ifg2_coords

def compare_hotspots(ifg1_hotspots, ifg2_hotspots, ifg1_coords, ifg2_coords, lig_cst, tol, ifg1name, ifg2name):
    ifg1_coords = np.array([ifg1_coords[spot[0]][spot[1]] for hotspot in ifg1_hotspots.copy() for spot in hotspot])
    ifg2_coords = np.array([ifg2_coords[spot[0]][spot[1]] for hotspot in ifg2_hotspots.copy() for spot in hotspot])
    nbrs_high = NearestNeighbors(radius=lig_cst.reshape(-1)[0] + tol, metric='euclidean')
    nbrs_low = NearestNeighbors(radius=lig_cst.reshape(-1)[0] - tol, metric='euclidean')
    nbrs_high.fit(ifg2_coords.reshape(ifg2_coords.shape[0], np.prod(ifg2_coords.shape[1:]))[:, :3])
    nbrs_low.fit(ifg2_coords.reshape(ifg2_coords.shape[0], np.prod(ifg2_coords.shape[1:]))[:, :3])
    rng_high = nbrs_high.radius_neighbors(ifg1_coords.reshape(ifg1_coords.shape[0], np.prod(ifg1_coords.shape[1:]))[:, :3])
    rng_low = nbrs_low.radius_neighbors(ifg1_coords.reshape(ifg1_coords.shape[0], np.prod(ifg1_coords.shape[1:]))[:, :3])
    uniqueindices = [np.setdiff1d(rng_high[1][i], rng_low[1][i], assume_unique=True) for i in range(len(rng_high[1]))]
    indices = []
    for i in range(len(ifg1_coords)):
        if uniqueindices[i].size != 0:
            dist = spatial.distance.cdist(ifg1_coords[i].reshape(-1, 3), ifg2_coords[uniqueindices[i]].reshape(-1, 3))
            dist_split = np.hsplit(dist, len(ifg2_coords[uniqueindices[i]]))
            dist = np.array([d.reshape(1, int(np.prod(d.shape))) for d in dist_split])
            lig_cst = lig_cst.ravel()
            test = (lig_cst - tol < dist) & (dist < lig_cst + tol)
            wh = np.where(np.all(test.T, axis=0).T)
            indices.extend(list(zip(uniqueindices[i][wh[0]], [i]*len(wh[0]))))
    indices = np.unique(indices, axis=0)
    ifg1_coords = np.array([[hs, spot[0], spot[1], ifg1_coords[spot[0]][spot[1]]] for hs, hotspot in
                            enumerate(ifg1_hotspots) for spot in hotspot], dtype=object)
    ifg2_coords = np.array([[hs, spot[0], spot[1], ifg2_coords[spot[0]][spot[1]]] for hs, hotspot in
                            enumerate(ifg2_hotspots) for spot in hotspot], dtype=object)
    return pd.DataFrame(indices,
                        columns=[ifg2name, ifg1name]), ifg1_coords, ifg2_coords

def combine_hotspots(df_list):
    #combine pandas dfs for mutually consistent hotspots
    column_set = np.unique([item for df in df_list for item in df.columns.values])
    tomerge = [df for df in df_list if column_set[0] in df.columns.values]
    mer1 = functools.reduce(lambda left,right: pd.merge(left, right, on=column_set[0]), tomerge)
    notmerged = [df for df in df_list if column_set[0] not in df.columns.values]
    if notmerged:
        columns = notmerged[0].columns.values.tolist()
        combined_hotspots = mer1.merge(notmerged[0], on=columns)
        if len(notmerged) > 1:
            for nm in notmerged[1:]:
                columns = notmerged[0].columns.values.tolist()
                combined_hotspots = combined_hotspots.merge(nm, on=columns)
    return combined_hotspots

def superpose_ligand(lig_sel, lig_mob_coords, combined_hotspots, ifg_coords_dict, ifg_key_order, poi, clash_cutoff):
    #superpose ligands on the mc hotspots
    #remove any ligands that clash with backbone
    # needs ifg_coords_dict with keys same as in combined_hotspots
    rows = []
    ligs = []
    for row in range(len(combined_hotspots)):
        tostack = [ifg_coords_dict[key][combined_hotspots[key][row], 4] for key in ifg_key_order]
        targ_coords = functools.reduce(lambda a, b: np.vstack(a, b), tostack)
        new_lig = superpose_coords(lig_sel, lig_mob_coords, targ_coords)
        poi_coords = get_poi_coords_for_clashsc(new_lig, poi)
        if calc_clash_score(poi_coords, new_lig.select('not element H D'), clash_cutoff=clash_cutoff) < 1:
            ligs.append(new_lig)
            rows.append(row)
    return combined_hotspots.loc[rows], ligs


def hotspots_with_ligand(g, ifg_coords_dict, combined_hotspots, ligs, rmsd_cutoff=0.5):
#retrieve hotspots of ifgs of ligands.
    #find mc sets of vdms within and between each hotspot. remove protein clashes and ligand clashes.
    bs_num_ifg_coords = {}
    nbrs = {}
    for col in combined_hotspots.columns:
        bs_num_ifg_coords[col] = [coords.flatten() for coords in ifg_coords_dict[col][:, -1]]
        nbrs[col] = NearestNeighbors(radius=rmsd_cutoff, metric='euclidean')
        nbrs[col].fit(bs_num_ifg_coords[col])

    for row, lig in enumerate(ligs):
        hotspots = [(col, combined_hotspots[col][row]) for col in combined_hotspots.columns]
        for hotspot in hotspots:
            # ifg_hs = list(nx.connected_components(g[hotspot[0]]))[hotspot[1]]
            rng = nbrs[hotspot[0]].radius_neighbors(ifg_coords_dict[hotspot[0]][hotspot[1]][-1].reshape(-1))
            num_res = len(np.unique(ifg_coords_dict[hotspot[0]][rng, 1]))
            grouped_vdms = [list(g) for k,g in itertools.groupby(sorted(ifg_coords_dict[hotspot[0]][rng, [1,2]].tolist()), lambda x: x[0])]
            combinations = {}
            for i in range(1, num_res+1):
                poss_combs = itertools.product(*grouped_vdms[:i])
                for each in poss_combs:

                    spatial.distance.cdist()
                # combinations[i] =

            # remove vdm in hotspot if clashing with ligand or protein backbone
            # store remaining vdms
            # need a step where only considering vdms in a hotspot that have ifgs within a tighter rmsd of this ifg
        #from stored vdms, make all collections with different residue positions.  Can keep same residue position for
        #different ifgs if it is the same amino acid and same rotamer.

    # rows = []
    # for row in range(len(combined_hotspots)):
    #     tostack = [ifg_coords_dict[key][combined_hotspots[key][row], 4] for key in ifg_key_order]
    #     targ_coords = functools.reduce(lambda a, b: np.vstack(a, b), tostack)
    #     new_lig = superpose_coords(lig_sel, lig_mob_coords, targ_coords)
    #     poi_coords = get_poi_coords_for_clashsc(new_lig, poi)
    #     if calc_clash_score(poi_coords, new_lig.select('not element H D'), clash_cutoff=clash_cutoff) < 1:
    #         ligs.append(new_lig)
    #         rows.append(row)
    # return combined_hotspots.iloc[rows], ligs

# def #print the candidates to file

def superpose_coords(mob_obj, mob_coords, targ_coords):
    mob_sel = mob_coords
    targ_sel = targ_coords
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
    mob_obj_cen_rot = np.dot((mob_obj - mob_sel_com), R)
    mob_obj_cen_rot_targ = mob_obj_cen_rot + targ_sel_com
    mob_obj_copy.setCoords(mob_obj_cen_rot_targ)
    return mob_obj_copy

    # wh = np.where((lig_cst - tol < dist) & (dist < lig_cst + tol))
    # print(wh)
    # ifg1_ind, counts = np.unique(wh[0], return_counts=True)
    # cumsum = np.cumsum(counts)
    # ifg2_ind = [np.unique(wh[1][i:j]) for i, j in zip(np.insert(cumsum, 0, 0), cumsum)]
    # return pd.DataFrame([(i, k) for i, j in zip(ifg1_ind, ifg2_ind) for k in j], columns=[ifg1name, ifg2name]), ifg1_coords, ifg2_coords

def get_ligand_hotspots():
    pass


def print_hotspot_pdbs_bb(hotspots, rois, bs_num_vdm_tags, interactamers, inpath, outpath, top_spots):
    roi_coords = {}
    r = {}
    int_com = {}
    roi_com = {}
    interactamers = np.vstack((interactamers.values()))
    for i, roi in enumerate(rois):
        # resn = roi.getResnames()[0]
        # for phi_low, phi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
        #     for psi_low, psi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
        #         if (phi_low <= phi_psis[i][0] < phi_high) and (psi_low <= phi_psis[i][1] < psi_high):
        #             phi_psi_key = (phi_low, phi_high, psi_low, psi_high)
        # r[i] = {}
        # roi_coords[i] = {}
        # int_com[i] = {}
        # roi_com[i] = {}
        roi_coords[i] = get_roi_coords(roi, resn='BB_NH')
        r[i], int_com[i], roi_com[i] = get_rotation_translation(interactamers[0, 4], roi_coords[i])
        # roi_coords[i]['BB_CO'] = get_roi_coords(roi, resn='BB_CO')
        # r[i]['BB_CO'], int_com[i]['BB_CO'], roi_com[i]['BB_CO'] = get_rotation_translation(interactamers['BB_CO']['O_not_N'][0, 4],
        #                                                                                     roi_coords[i]['BB_CO'])
        # roi_coords[i]['BB_NH'] = get_roi_coords(roi, resn='BB_NH')
        # r[i]['BB_NH'], int_com[i]['BB_NH'], roi_com[i]['BB_NH'] = get_rotation_translation(interactamers['BB_NH']['N_not_O'][0, 4],
        #                                                                                     roi_coords[i]['BB_NH'])
    for spot_rank, spot in enumerate(sorted(hotspots, key=len, reverse=True)):
        if spot_rank+1 <= top_spots:
            pdbs = {}
            for i, roi in enumerate(rois):
                pdbs[i] = {}
                # resn = roi.getResnames()[0]
                pdbs[i] = []
                # pdbs[i]['BB_CO'] = []
                # pdbs[i]['BB_NH'] = []
            for indices in spot:
                # print(bs_num_vdm_tags[indices[0]])
                tags = [str(t) for t in bs_num_vdm_tags[indices[0]][indices[1]]]
                resn = tags[-1]
                subfolder = ''
                if 'BB' in resn:
                    resn = resn[:5]
                    subfolder = '/' + tags[-1][6:]
                # print(tags)
                pdb = [file for file in os.listdir(inpath + resn + subfolder + '/') if 'iFG_' + tags[0] + '_vdM_' + tags[1]
                       + '_vFlip_' + tags[2] + '_iFlip_' + tags[3] in file]
                # print(tags)
                # print(subfolder)
                # print(pdb)
                pdb=pdb[0]
                pdbs[indices[0]].append(pr.parsePDB(inpath + resn + subfolder + '/' + pdb))
            for bs_site in pdbs.keys():
                # for resn in pdbs[bs_site].keys():
                movepdbs = pdbs[bs_site]
                for movepdb in movepdbs:
                    movepdb.setCoords(move_interactamers(movepdb.getCoords(), r[bs_site], int_com[bs_site],
                                                         roi_com[bs_site]))
                    try:
                        os.makedirs(outpath + 'hotspots/' + str(spot_rank+1))
                    except:
                        pass
                    pr.writePDB(outpath + 'hotspots/' + str(spot_rank + 1) + '/hs' + str(spot_rank + 1) + '_'
                                            + str(bs_site) + '_' + repr(movepdb).split()[1]
                                            + '_superposed.pdb.gz', movepdb)

def print_hotspot_pdbs(hotspots, rois, bs_num_vdm_tags, interactamers, inpath, outpath, top_spots):
    roi_coords = {}
    r = {}
    int_com = {}
    roi_com = {}
    for i, roi in enumerate(rois):
        resn = roi.getResnames()[0]
        # for phi_low, phi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
        #     for psi_low, psi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
        #         if (phi_low <= phi_psis[i][0] < phi_high) and (psi_low <= phi_psis[i][1] < psi_high):
        #             phi_psi_key = (phi_low, phi_high, psi_low, psi_high)
        r[i] = {}
        roi_coords[i] = {}
        int_com[i] = {}
        roi_com[i] = {}
        roi_coords[i][resn] = get_roi_coords(roi)
        r[i][resn], int_com[i][resn], roi_com[i][resn] = get_rotation_translation(interactamers[resn][0, 4], roi_coords[i][resn])
        roi_coords[i]['BB_CO'] = get_roi_coords(roi, resn='BB_CO')
        r[i]['BB_CO'], int_com[i]['BB_CO'], roi_com[i]['BB_CO'] = get_rotation_translation(interactamers['BB_CO']['O_not_N'][0, 4],
                                                                                            roi_coords[i]['BB_CO'])
        roi_coords[i]['BB_NH'] = get_roi_coords(roi, resn='BB_NH')
        r[i]['BB_NH'], int_com[i]['BB_NH'], roi_com[i]['BB_NH'] = get_rotation_translation(interactamers['BB_NH']['N_not_O'][0, 4],
                                                                                            roi_coords[i]['BB_NH'])
    for spot_rank, spot in enumerate(sorted(hotspots, key=len, reverse=True)):
        if spot_rank+1 <= top_spots:
            pdbs = {}
            for i, roi in enumerate(rois):
                pdbs[i] = {}
                resn = roi.getResnames()[0]
                pdbs[i][resn] = []
                pdbs[i]['BB_CO'] = []
                pdbs[i]['BB_NH'] = []
            for indices in spot:
                # print(bs_num_vdm_tags[indices[0]])
                tags = [str(t) for t in bs_num_vdm_tags[indices[0]][indices[1]]]
                resn = tags[-1]
                subfolder = ''
                if 'BB' in resn:
                    resn = resn[:5]
                    subfolder = '/' + tags[-1][6:]
                # print(tags)
                pdb = [file for file in os.listdir(inpath + resn + subfolder + '/') if 'iFG_' + tags[0] + '_vdM_' + tags[1]
                       + '_vFlip_' + tags[2] + '_iFlip_' + tags[3] in file]
                # print(tags)
                # print(subfolder)
                # print(pdb)
                pdb=pdb[0]
                pdbs[indices[0]][resn].append(pr.parsePDB(inpath + resn + subfolder + '/' + pdb))
            for bs_site in pdbs.keys():
                for resn in pdbs[bs_site].keys():
                    movepdbs = pdbs[bs_site][resn]
                    for movepdb in movepdbs:
                        movepdb.setCoords(move_interactamers(movepdb.getCoords(), r[bs_site][resn], int_com[bs_site][resn],
                                                             roi_com[bs_site][resn]))
                        try:
                            os.makedirs(outpath + 'hotspots/' + str(spot_rank+1))
                        except:
                            pass
                        pr.writePDB(outpath + 'hotspots/' + str(spot_rank + 1) + '/'
                                                + str(bs_site) + '_' + repr(movepdb).split()[1]
                                                + '_superposed.pdb.gz', movepdb)
#     #parse pdb. translate it. print to file
#     if path[-1] != '/':
#         path += '/'
#     for i, spot in enumerate(sorted(hotspots, key=len, reverse=True)):
#         try:
#             os.makedirs(path + 'hotspots/' + str(i+1))
#         except:
#             pass
#         for indices in spot:
#             tags = bs_num_vdm_tags[spot[0]][spot[1]]
#             parsePDB(inpath +
#             #NEED TO KEEP TRACK OF WHAT PDB FOLDER TO GRAB PDBS FROM
#
#             pr.writePDB(path + 'hotspots/' + str(i+1) + '/'
#                         + str(indices[0]) + '_' + repr(bs_num_pdbs[indices[0]][indices[1]]).split()[1]
#                         + '_superposed.pdb.gz', bs_num_pdbs[indices[0]][indices[1]])
