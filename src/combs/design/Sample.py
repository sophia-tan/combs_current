__all__ = ['Sample']

import prody as pr
from scipy.spatial.distance import cdist
import numpy as np
import collections
from ..apps.constants import bb_type_dict
from sklearn.neighbors import NearestNeighbors
import networkx as nx
import itertools
import random
import functools
from sklearn.neighbors.kde import KernelDensity
from ..apps.fitting import score_fit, score_fit_1, rigid_body, rigid_body_highdim
from scipy.optimize import minimize, basinhopping
from .pose import Pose
from multiprocessing import Pool, Process, JoinableQueue, Manager
from scipy import stats

class Sample:
    def __init__(self):
        self.poi = None
        self.bs_residues = None
        self.rois = {}
        self.rois_bb_coords = collections.defaultdict(dict)
        self.rois_phipsi = {}
        self.ligand_conformers = None
        self.rel_bb_coords = None
        self.poi_clash_coords_NO = {}
        self.poi_clash_coords_C = {}
        self.lig_ifg_corr = collections.defaultdict(dict)
        self.lig_ifg_coords = collections.defaultdict(dict)
        self.hotspot_cliques = {}
        self.lig_names = {}
        # self.ifg_names = {}
        self.ifg_names = collections.defaultdict(dict)
        self.lig_ifg_dict = {}
        self.rel_vdms = {}
        self.poses = []
        # self._results = None
        # self._sorted_result_indices = None
        self.pose_scores = None
        self.ranked_poses_by_score = None
        self.pose_num_satisfied_iFGs = None
        self.name_pairs = {}

    def set_rois(self):
        """returns a list of prody selections (backbone only) of the residues of interest from the protein
        of interest"""
        for resnum, chid in self.bs_residues:
            self.rois[(resnum, chid)] = self.poi.select('(backbone or name H H1) and resnum ' + str(resnum)
                                                        + ' and chain ' + chid)
            if self.rois[(resnum, chid)].select('name CA').getResnames() == 'PRO':
                self.rois[(resnum, chid)] = self.poi.select('(backbone or name CD) and resnum ' + str(resnum)
                                                            + ' and chain ' + chid)


    def set_rois_phipsi(self):
        """returns a list of the phi,psi tuples of each residue of interest"""
        for resnum, chid in self.bs_residues:
            try:
                phi = pr.calcPhi(self.poi[chid, resnum])
            except:
                phi = None
            try:
                psi = pr.calcPsi(self.poi[chid, resnum])
            except:
                psi = None
            self.rois_phipsi[(resnum, chid)] = (phi, psi)

    def set_poi_clash_coords(self):
        """take prody selection roi. Returns an array of heavy atom backbone coordinates within 15A
        of the residue of interest"""
        # # for resnum_chid, roi in self.rois.items():
        # #     coords_sel = self.poi.select('(backbone and not element H D) exwithin 15 of sel', sel=roi)
        # #     self.poi_clash_coords[resnum_chid] = coords_sel.getCoords()
        # for resnum_chid, roi in self.rois.items():
        #     coords_sel = self.poi.select('(backbone and not element H D) exwithin 15 of sel', sel=roi)
        #     self.poi_clash_coords[resnum_chid] = NearestNeighbors(metric='euclidean', radius=rmsd_cutoff, n_neighbors=1)
        #     self.poi_clash_coords[resnum_chid].fit(coords_sel.getCoords())
        for resnum_chid, roi in self.rois.items():
            coords_sel_NO = self.poi.select('(name N O) exwithin 15 of sel', sel=roi)
            coords_sel_C = self.poi.select('(name C CA) exwithin 15 of sel', sel=roi)
            self.poi_clash_coords_NO[resnum_chid] = NearestNeighbors(metric='euclidean', n_neighbors=1)
            self.poi_clash_coords_NO[resnum_chid].fit(coords_sel_NO.getCoords())
            self.poi_clash_coords_C[resnum_chid] = NearestNeighbors(metric='euclidean', n_neighbors=1)
            self.poi_clash_coords_C[resnum_chid].fit(coords_sel_C.getCoords())

    def set_roi_bb_coords(self):
        """bb_atom_list order follows same convention as order in make_rel_vdm_coords"""
        for resnum_chid, roi in self.rois.items():
            for type_, bb_atom_list in bb_type_dict.items():
                try:
                    self.rois_bb_coords[resnum_chid][type_] = np.array([roi.select('name ' + name).getCoords()[0]
                                                                       for name in bb_atom_list])
                except:
                    # print(str(roi))
                    print(bb_atom_list)
                    print(resnum_chid)
                    if type_ == 'N_CA':
                        try:
                            self.rois_bb_coords[resnum_chid][type_] = np.array([roi.select('name ' + name).getCoords()[0]
                                                                            for name in ['N', 'H1', 'CA']])
                        except:
                            self.rois_bb_coords[resnum_chid][type_] = np.array(
                                [roi.select('name ' + name).getCoords()[0]
                                 for name in ['N', 'CD', 'CA']])

    # def set_ligand_ifg_correspondence(self, rel_vdm, ligand_names, ifg_names):
    #     self.lig_names[rel_vdm.name] = ligand_names.split()
    #     self.ifg_names[rel_vdm.name] = ifg_names.split()
    #     for lig, ifg in zip(ligand_names.split(), ifg_names.split()):
    #         self.lig_ifg_corr[rel_vdm.name][lig] = ifg

        # for conformer_num, conformer_sel in enumerate(self.ligand_conformers):
        #     self.lig_ifg_coords[rel_vdm.name][conformer_num] = \
        #         conformer_sel.select('name ' + ligand_name).getCoords()[0]

    # def set_ligand_ifg_coords(self, rel_vdm):
    #     for conformer_num, conformer_sel in enumerate(self.ligand_conformers):
    #         self.lig_ifg_coords[rel_vdm.name][conformer_num] = {}
    #         for lig in self.lig_ifg_corr[rel_vdm.name].keys():
    #             self.lig_ifg_coords[rel_vdm.name][conformer_num][lig] = \
    #             conformer_sel.select('name ' + lig).getCoords()[0]

    def set_ligand_ifg_coords(self, rel_vdm):
        for conformer_num, conformer_sel in enumerate(self.ligand_conformers):
            self.lig_ifg_coords[rel_vdm.name][conformer_num] = {}
            for lig in self.lig_ifg_dict[rel_vdm.name]:
                self.lig_ifg_coords[rel_vdm.name][conformer_num][lig] = \
                conformer_sel.select('name ' + lig).getCoords()[0]


    # def set_ligand_ifg_correspondence(self, rel_vdm, ligand_names, ifg_name):
    #     self.lig_ifg_corr[rel_vdm.name] = {ifg_name: ligand_names}
    #     for conformer_num, conformer_sel in enumerate(self.ligand_conformers):
    #         self.lig_ifg_coords[rel_vdm.name][conformer_num] = \
    #             np.array([conformer_sel.select('name ' + ligand_name).getCoords()[0] for ligand_name in ligand_names.split()])

    def set_ligand_ifg_distance_csts(self, rel_vdm1, rel_vdm2, lig_key1, lig_key2):
        return [cdist(self.lig_ifg_coords[rel_vdm1.name][i][lig_key1].reshape(1, -1),
                      self.lig_ifg_coords[rel_vdm2.name][i][lig_key2].reshape(1, -1))[0][0]
                for i in range(len(self.ligand_conformers))]

    @staticmethod
    def _neigh_internal(hs_reps1, hs_reps2, dist_cst, tol):
        nbrs_high = NearestNeighbors(metric='euclidean', radius=dist_cst + tol)
        nbrs_high.fit(hs_reps2)
        adj_mat_high = nbrs_high.radius_neighbors_graph(hs_reps1, mode='distance')
        lowtol = dist_cst - tol
        wh = (adj_mat_high > lowtol).nonzero()
        return zip(wh[0], wh[1])

    @staticmethod
    def _set_hs_reps(rel_vdm, index):
        return np.array([random.choice(list(nx.get_node_attributes(graph, 'ifg').values()))
                        for graph in rel_vdm.hotspot_subgraphs])[:, 3 * index:3 * index + 3]

    # @staticmethod
    # def _set_hs_reps(rel_vdm, index):
    #     hs_reps = np.array(functools.reduce(lambda x, y: np.vstack((x, y)),
    #                                      (np.array(list(nx.get_node_attributes(graph, 'ifg').values()))[:, 3 * index:3 * index + 3]
    #                                       for graph in rel_vdm.hotspot_subgraphs)))
    #     hs_reps_ind = np.array(functools.reduce(lambda x, y: np.hstack((x, y)),
    #                                      (np.array([i] * len(nx.get_node_attributes(graph, 'ifg').keys()))
    #                                       for i, graph in enumerate(rel_vdm.hotspot_subgraphs))))
    #     return hs_reps, hs_reps_ind

    # def find_ligand_cst_hotspots(self, tol=0.2):
    #     g = {}
    #     hotspot_cliques = {}
    #     max_clique_len = len(self.ifg_names) #len(list(itertools.combinations(self.rel_vdms.keys(), 2)))
    #     # max_clique_len = len(list(list(self.ifg_names.values())[0].values())[0])
    #     for y in range(max_clique_len):
    #         g[y] = {}
    #         hotspot_cliques[y] = {}
    #
    #         for i in range(len(self.ligand_conformers)):
    #             g[y][i] = nx.Graph()
    #
    #         for rel_vdm1, rel_vdm2 in itertools.combinations(self.rel_vdms.values(), 2):
    #             # names = list(zip(list(self.lig_ifg_corr[rel_vdm1.name].keys()),
    #             #                list(self.lig_ifg_corr[rel_vdm2.name].keys()),
    #             #                list(self.lig_ifg_corr[rel_vdm1.name].values()),
    #             #                list(self.lig_ifg_corr[rel_vdm2.name].values())))
    #             names = list(zip(self.lig_names[rel_vdm1.name],
    #                              self.lig_names[rel_vdm2.name],
    #                              self.ifg_names[rel_vdm1.name],
    #                              self.ifg_names[rel_vdm2.name]))
    #             lig1, lig2, ifg1, ifg2 = names[y]
    #             print(names[y])
    #             dist_csts = self.set_ligand_ifg_distance_csts(rel_vdm1, rel_vdm2, lig1, lig2)
    #             index1 = list(rel_vdm1.ifg_dict.values())[0].split().index(ifg1)
    #             index2 = list(rel_vdm2.ifg_dict.values())[0].split().index(ifg2)
    #             for k, dist_cst in enumerate(dist_csts):
    #                 hs_reps1 = self._set_hs_reps(rel_vdm1, index1)
    #                 hs_reps2 = self._set_hs_reps(rel_vdm2, index2)
    #                 hs_pair_set = set(self._neigh_internal(hs_reps1, hs_reps2, dist_cst, tol))
    #                 for j in range(1000):
    #                     print('trial ', j)
    #                     hs_reps1 = self._set_hs_reps(rel_vdm1, index1)
    #                     hs_reps2 = self._set_hs_reps(rel_vdm2, index2)
    #                     hs_pair_set_test = hs_pair_set.union(set(self._neigh_internal(hs_reps1, hs_reps2, dist_cst, tol)))
    #                     if hs_pair_set_test == hs_pair_set:
    #                         hs_pair_set = hs_pair_set_test
    #                         hs_reps1 = self._set_hs_reps(rel_vdm1, index1)
    #                         hs_reps2 = self._set_hs_reps(rel_vdm2, index2)
    #                         hs_pair_set_test = hs_pair_set.union(
    #                             set(self._neigh_internal(hs_reps1, hs_reps2, dist_cst, tol)))
    #                         if hs_pair_set_test == hs_pair_set:
    #                             print('hotspot pairs (' + rel_vdm1.name + ',' + rel_vdm2.name + ') converged after '
    #                                   + str(j+1) + ' trials')
    #                             break
    #                     hs_pair_set = hs_pair_set_test
    #                 wh = np.array(list(hs_pair_set))
    #                 # delete 3rd item below to return to normal
    #                 len_ = len(wh[:, 0])
    #                 # rel_vdm1_nodes = zip(wh[:, 0], [rel_vdm1.name] * len_, [(lig1, ifg1)] * len_)
    #                 # rel_vdm2_nodes = zip(wh[:, 1], [rel_vdm2.name] * len_, [(lig2, ifg2)] * len_)
    #                 rel_vdm1_nodes = zip(wh[:, 0], [rel_vdm1.name] * len_)
    #                 rel_vdm2_nodes = zip(wh[:, 1], [rel_vdm2.name] * len_)
    #                 g[y][k].add_nodes_from(rel_vdm1_nodes)
    #                 g[y][k].add_nodes_from(rel_vdm2_nodes)
    #                 # rel_vdm1_nodes = zip(wh[:, 0], [rel_vdm1.name] * len_, [(lig1, ifg1)] * len_)
    #                 # rel_vdm2_nodes = zip(wh[:, 1], [rel_vdm2.name] * len_, [(lig2, ifg2)] * len_)
    #                 rel_vdm1_nodes = zip(wh[:, 0], [rel_vdm1.name] * len_)
    #                 rel_vdm2_nodes = zip(wh[:, 1], [rel_vdm2.name] * len_)
    #                 edges = zip(rel_vdm1_nodes, rel_vdm2_nodes)
    #                 g[y][k].add_edges_from(edges)
    #         for conf_num, graph_ in g[y].items():
    #             # hotspot_cliques[y][conf_num] = list(filter(lambda x: len(x) == max_clique_len, nx.find_cliques(graph_)))
    #             hotspot_cliques[y][conf_num] = set(tuple(sorted(m)) for m in filter(lambda x: len(x) == max_clique_len,
    #                                                                          nx.find_cliques(graph_)))
    #     for conf_num in range(len(self.ligand_conformers)):
    #         self.hotspot_cliques[conf_num] = list(functools.reduce(lambda a, b: a.intersection(b),
    #                          (hotspot_cliques[key][conf_num] for key in hotspot_cliques.keys())))

    def set_lig_ifg_dict(self, rel_vdm, lig_names):
        self.lig_ifg_dict[rel_vdm.name] = lig_names.split()
        ifgnames = list(rel_vdm.ifg_dict.values())[0].split()
        for lname, iname in zip(self.lig_ifg_dict[rel_vdm.name], ifgnames):
            self.ifg_names[rel_vdm.name][lname] = iname


    def find_ligand_cst_hotspots(self):
        g = {}
        hotspot_cliques = {}
        max_clique_len = len(self.ifg_names) #len(list(itertools.combinations(self.rel_vdms.keys(), 2)))
        # max_clique_len = len(list(list(self.ifg_names.values())[0].values())[0])
        for y in range(max_clique_len):
            g[y] = {}
            hotspot_cliques[y] = {}

            for i in range(len(self.ligand_conformers)):
                g[y][i] = nx.Graph()

            for rel_vdm1, rel_vdm2 in itertools.combinations(self.rel_vdms.values(), 2):
                # names = list(zip(list(self.lig_ifg_corr[rel_vdm1.name].keys()),
                #                list(self.lig_ifg_corr[rel_vdm2.name].keys()),
                #                list(self.lig_ifg_corr[rel_vdm1.name].values()),
                #                list(self.lig_ifg_corr[rel_vdm2.name].values())))
                # names = list(zip(self.lig_names[rel_vdm1.name],
                #                  self.lig_names[rel_vdm2.name],
                #                  self.ifg_names[rel_vdm1.name],
                #                  self.ifg_names[rel_vdm2.name]))
                wh_sets = []
                for name_pair, tol in self.name_pairs[(rel_vdm1.name, rel_vdm2.name)]:
                    print(name_pair)
                    ifg1 = self.ifg_names[rel_vdm1.name][name_pair[0]]
                    ifg2 = self.ifg_names[rel_vdm2.name][name_pair[1]]
                    print('ifg1= ', ifg1, ' ifg2= ', ifg2)
                    dist_csts = self.set_ligand_ifg_distance_csts(rel_vdm1, rel_vdm2, name_pair[0], name_pair[1]) # needs to take the name
                    print('dist_csts= ', dist_csts)
                    index1 = list(rel_vdm1.ifg_dict.values())[0].split().index(ifg1)
                    index2 = list(rel_vdm2.ifg_dict.values())[0].split().index(ifg2)
                    for k, dist_cst in enumerate(dist_csts):
                        hs_reps1 = self._set_hs_reps(rel_vdm1, index1)
                        hs_reps2 = self._set_hs_reps(rel_vdm2, index2)
                        hs_pair_set = set(self._neigh_internal(hs_reps1, hs_reps2, dist_cst, tol))
                        for j in range(1000):
                            print('trial ', j)
                            hs_reps1 = self._set_hs_reps(rel_vdm1, index1)
                            hs_reps2 = self._set_hs_reps(rel_vdm2, index2)
                            hs_pair_set_test = hs_pair_set.union(set(self._neigh_internal(hs_reps1, hs_reps2, dist_cst, tol)))
                            if hs_pair_set_test == hs_pair_set:
                                hs_pair_set = hs_pair_set_test
                                hs_reps1 = self._set_hs_reps(rel_vdm1, index1)
                                hs_reps2 = self._set_hs_reps(rel_vdm2, index2)
                                hs_pair_set_test = hs_pair_set.union(
                                    set(self._neigh_internal(hs_reps1, hs_reps2, dist_cst, tol)))
                                if hs_pair_set_test == hs_pair_set:
                                    print('hotspot pairs (' + rel_vdm1.name + ',' + rel_vdm2.name + ') converged after '
                                          + str(j+1) + ' trials')
                                    break
                            hs_pair_set = hs_pair_set_test
                    wh_sets.append(hs_pair_set)
                hs_pair_set = set.intersection(*wh_sets)
                wh = np.array(list(hs_pair_set))
                # delete 3rd item below to return to normal
                len_ = len(wh[:, 0])
                # rel_vdm1_nodes = zip(wh[:, 0], [rel_vdm1.name] * len_, [(lig1, ifg1)] * len_)
                # rel_vdm2_nodes = zip(wh[:, 1], [rel_vdm2.name] * len_, [(lig2, ifg2)] * len_)
                rel_vdm1_nodes = zip(wh[:, 0], [rel_vdm1.name] * len_)
                rel_vdm2_nodes = zip(wh[:, 1], [rel_vdm2.name] * len_)
                g[y][k].add_nodes_from(rel_vdm1_nodes)
                g[y][k].add_nodes_from(rel_vdm2_nodes)
                # rel_vdm1_nodes = zip(wh[:, 0], [rel_vdm1.name] * len_, [(lig1, ifg1)] * len_)
                # rel_vdm2_nodes = zip(wh[:, 1], [rel_vdm2.name] * len_, [(lig2, ifg2)] * len_)
                rel_vdm1_nodes = zip(wh[:, 0], [rel_vdm1.name] * len_)
                rel_vdm2_nodes = zip(wh[:, 1], [rel_vdm2.name] * len_)
                edges = zip(rel_vdm1_nodes, rel_vdm2_nodes)
                g[y][k].add_edges_from(edges)
            for conf_num, graph_ in g[y].items():
                # hotspot_cliques[y][conf_num] = list(filter(lambda x: len(x) == max_clique_len, nx.find_cliques(graph_)))
                hotspot_cliques[y][conf_num] = set(tuple(sorted(m)) for m in filter(lambda x: len(x) == max_clique_len,
                                                                             nx.find_cliques(graph_)))
        for conf_num in range(len(self.ligand_conformers)):
            self.hotspot_cliques[conf_num] = list(functools.reduce(lambda a, b: a.intersection(b),
                             (hotspot_cliques[key][conf_num] for key in hotspot_cliques.keys())))

    # def find_ligand_cst_hotspots(self, tol=0.2):
    #     g = {}
    #     hotspot_cliques = {}
    #     max_clique_len = len(list(itertools.combinations(self.rel_vdms.keys(), 2)))
    #
    #     for y in range(len(self.ifg_names)):
    #         g[y] = {}
    #         hotspot_cliques[y] = {}
    #
    #         for i in range(len(self.ligand_conformers)):
    #             g[y][i] = nx.Graph()
    #
    #         for rel_vdm1, rel_vdm2 in itertools.combinations(self.rel_vdms.values(), 2):
    #             # names = list(zip(list(self.lig_ifg_corr[rel_vdm1.name].keys()),
    #             #                list(self.lig_ifg_corr[rel_vdm2.name].keys()),
    #             #                list(self.lig_ifg_corr[rel_vdm1.name].values()),
    #             #                list(self.lig_ifg_corr[rel_vdm2.name].values())))
    #             names = list(zip(self.lig_names[rel_vdm1.name],
    #                              self.lig_names[rel_vdm2.name],
    #                              self.ifg_names[rel_vdm1.name],
    #                              self.ifg_names[rel_vdm2.name]))
    #             lig1, lig2, ifg1, ifg2 = names[y]
    #             print(names[y])
    #             dist_csts = self.set_ligand_ifg_distance_csts(rel_vdm1, rel_vdm2, lig1, lig2)
    #             index1 = list(rel_vdm1.ifg_dict.values())[0].split().index(ifg1)
    #             index2 = list(rel_vdm2.ifg_dict.values())[0].split().index(ifg2)
    #             for k, dist_cst in enumerate(dist_csts):
    #                 hs_reps1, hs_reps1_ind = self._set_hs_reps(rel_vdm1, index1)
    #                 hs_reps2, hs_reps2_ind = self._set_hs_reps(rel_vdm2, index2)
    #                 wh = self._neigh_internal(hs_reps1, hs_reps2, dist_cst, tol)
    #
    #                 rel_vdm1_nodes = itertools.zip_longest(hs_reps1_ind[wh[0]], rel_vdm1.name, fillvalue=rel_vdm1.name)
    #                 rel_vdm2_nodes = itertools.zip_longest(hs_reps2_ind[wh[1]], rel_vdm2.name, fillvalue=rel_vdm2.name)
    #                 g[y][k].add_nodes_from(rel_vdm1_nodes)
    #                 g[y][k].add_nodes_from(rel_vdm2_nodes)
    #                 rel_vdm1_nodes = itertools.zip_longest(hs_reps1_ind[wh[0]], rel_vdm1.name, fillvalue=rel_vdm1.name)
    #                 rel_vdm2_nodes = itertools.zip_longest(hs_reps2_ind[wh[1]], rel_vdm2.name, fillvalue=rel_vdm2.name)
    #                 edges = zip(rel_vdm1_nodes, rel_vdm2_nodes)
    #                 g[y][k].add_edges_from(edges)
    #         for conf_num, graph_ in g[y].items():
    #             # hotspot_cliques[y][conf_num] = list(filter(lambda x: len(x) == max_clique_len, nx.find_cliques(graph_)))
    #             hotspot_cliques[y][conf_num] = set(tuple(sorted(m)) for m in filter(lambda x: len(x) == max_clique_len,
    #                                                                                 nx.find_cliques(graph_)))
    #     for conf_num in range(len(self.ligand_conformers)):
    #         self.hotspot_cliques[conf_num] = list(functools.reduce(lambda a, b: a.intersection(b),
    #                                                                (hotspot_cliques[key][conf_num] for key in
    #                                                                 hotspot_cliques.keys())))

    def make_densities(self):
        sorted_cliques = np.array([sorted(cli, key=lambda x: x[1]) for cliques in self.hotspot_cliques.values()
                                   for cli in cliques])  #sort by relvdm.name
        for relvdm in self.rel_vdms.values():
            hs_ind = [i for i, tup in enumerate(sorted_cliques[0]) if tup[1] == relvdm.name][0]
            hs_nums = np.array(np.unique(sorted_cliques[:, hs_ind, 0]), dtype=int)
            func = functools.partial(self.set_density, relvdm)
            list(map(func, hs_nums))

    @staticmethod
    def set_density(relvdm, hs_num):
        ifgs = np.array(list(nx.get_node_attributes(relvdm.hotspot_subgraphs[hs_num], 'ifg').values()))
        relvdm.hotspot_subgraphs[hs_num].graph['density'] = KernelDensity(kernel='gaussian', bandwidth=2.0).fit(ifgs)

    # def set_lig_ifg_dict(self, rel_vdm, lig_names):
    #     self.lig_ifg_dict[rel_vdm.name] = lig_names.split()  # e.g. 'CG CD OE1 NE2'

    # need to set lig_ifg_dict first: keys are relvdm.name, values are string with lig atom names
    # in direct correspondence to ifg names from relvdm.ifg_dict.values()
    def set_lig_coords(self):
        for relvdm_name, atom_names in self.lig_ifg_dict.items():
            for i in range(len(self.ligand_conformers)):
                self.lig_ifg_coords[i][relvdm_name] = \
                    np.array([self.ligand_conformers[i].select('name ' + name).getCoords()[0]
                              for name in atom_names]).flatten().reshape(1, -1)

    def set_protein_bb_coords(self):
        self.poi_bb_coords_NO = NearestNeighbors(metric='euclidean', n_neighbors=1)
        self.poi_bb_coords_NO.fit(self.poi.select('name N O').getCoords())
        self.poi_bb_coords_C = NearestNeighbors(metric='euclidean', n_neighbors=1)
        self.poi_bb_coords_C.fit(self.poi.select('name C CA').getCoords())

    def set_rel_vdms(self, relvdms):
        for relvdm in relvdms:
            self.rel_vdms[relvdm.name] = relvdm

    def fit_lig_to_density(self):
        # which densities
        # which coords
        # fit coords to density
        # move ligand
        # ligand clash with backbone?
        # store
        for conf_num, lig in enumerate(self.ligand_conformers):
            i=1
            for j, cliq in enumerate(self.hotspot_cliques[conf_num]):
                print('fitting clique ' + str(i))
                i += 1
                lig_coords = [self.lig_ifg_coords[conf_num][tup[1]] for tup in cliq]
                densities = [self.rel_vdms[tup[1]].hotspot_subgraphs[tup[0]].graph['density'] for tup in cliq]
                # lens = [len(self.rel_vdms[tup[1]].hotspot_subgraphs[tup[0]]) for tup in cliq]
                # partial_score_fit = functools.partial(score_fit, list(zip(lens, lig_coords, densities)))
                partial_score_fit = functools.partial(score_fit, list(zip(lig_coords, densities)))
                min_ = minimize(partial_score_fit, np.zeros(6), method='Nelder-Mead', options={'maxiter': 3000})
                # min_ = minimize(partial_score_fit, np.array([0, 0, 0, 1, 1, 1]))
                lig_copy = lig.copy()
                lig_min_coords_heavy = rigid_body(lig_copy.select('not element H D').getCoords(), min_.x)
                if self.bb_clash_test(lig_min_coords_heavy):
                    lig_min_coords_all = rigid_body(lig_copy.getCoords(), min_.x)
                    lig_copy.setCoords(lig_min_coords_all)
                    pose = Pose()
                    for tup in cliq:
                        pose.lig_ifg_coords[tup[1]] = rigid_body_highdim(self.lig_ifg_coords[conf_num][tup[1]], min_.x)
                    pose.ligand = lig_copy
                    pose.hotspot_clique = cliq
                    pose.optim_result = min_
                    pose.number = j
                    self.poses.append(pose)
                    # print(self.bb_clash_test(lig_min_coords_heavy))

    def bb_clash_test(self, test_coords_NO, test_coords_C, cutoff_NO_NO=2.5, cutoff_NO_C=2.8, cutoff_C_C=3.2):
        distances_NO_NO, indices_NO_NO = self.poi_bb_coords_NO.kneighbors(test_coords_NO)
        distances_NO_C, indices_NO_C = self.poi_bb_coords_NO.kneighbors(test_coords_C)
        distances_C_NO, indices_C_NO = self.poi_bb_coords_C.kneighbors(test_coords_NO)
        distances_C_C, indices_C_C = self.poi_bb_coords_C.kneighbors(test_coords_C)
        if (distances_NO_NO < cutoff_NO_NO).any() or (distances_NO_C < cutoff_NO_C).any() \
           or (distances_C_NO < cutoff_NO_C).any() or (distances_C_C < cutoff_C_C).any():
            return False
        else:
            return True

    # def set_results(self):
    #     self._results = [pose.optim_result.fun for pose in self.poses]

    def get_pose_scores(self):
        self.pose_scores = [pose.member_score for pose in self.poses]
        self.ranked_poses_by_score = sorted(range(len(self.pose_scores)), key=lambda k: self.pose_scores[k], reverse=True)
        self.pose_num_satisfied_iFGs = [pose.num_satisfied_iFGs for pose in self.poses]
        self.pose_num_residues = [pose.num_residues for pose in self.poses]

    # def set_sorted_results_ind(self):
    #     self._sorted_result_indices = sorted(range(len(self._results)), key=lambda k: self._results[k], reverse=True)
    #
    # def store_poses(self, poses):
    #     pass

    def pose_metrics(self, rmsd_cutoff=1.5):
        for pose in self.poses:
            pose.get_around_lig_fg(self, rmsd_cutoff=rmsd_cutoff)
            pose.get_nonclashing(clash_cutoff=2.5)
            pose.num_members = {}  # can delete this later
            for key, graph in pose.subgraphs.items():
                if graph:
                    pose.num_members[key] = len(graph)
                else:
                    pose.num_members[key] = 0
            pose.set_member_score()
            pose.total_residues = {node[0] for graph in pose.subgraphs.values() for node in graph}
            if pose.total_residues:
                pose.num_residues = len(pose.total_residues)
                pose.num_satisfied_iFGs = len([1 for gr in pose.subgraphs.values() if len(gr) > 0])
        self.get_pose_scores()


    @staticmethod
    def chunks(l, n):
        """Yield successive n-sized chunks from l."""
        for i in range(0, len(l), n):
            yield l[i:i + n]

    def get_lig_coords(self, relvdm, lig=None):
        lig = lig or self.ligand_conformers[0]
        names = self.ligand_conformers[0].select('name ' + relvdm.lig_ifg_correspondence).getNames()
        indices = [i for name in relvdm.lig_ifg_correspondence.split() for i, n in enumerate(names) if n == name]
        return lig.select('name ' + relvdm.lig_ifg_correspondence).getCoords()[indices]

    def make_lig_dmat(self, relvdm1, relvdm2):
        coords1 = self.get_lig_coords(relvdm1)
        coords2 = self.get_lig_coords(relvdm2)
        return cdist(coords1, coords2)

    # def find_cluster_pairs(self, relvdm1, relvdm2, **kwargs):
    #
    #     num_atoms_ifg1 = kwargs.get('num_atoms_ifg1', len(str(relvdm1.ifg_dict.values()).split()))
    #     num_atoms_ifg2 = kwargs.get('num_atoms_ifg2', len(str(relvdm2.ifg_dict.values()).split()))
    #     cents_chunks_1 = kwargs.get('cents_chunks_1', list(self.chunks(relvdm1._cents, 2000)))
    #     cents_chunks_2 = kwargs.get('cents_chunks_2', list(self.chunks(relvdm2._cents, 2000)))
    #     lig_dmat = kwargs.get('lig_dmat', self.make_lig_dmat(relvdm1, relvdm2))
    #     atol = kwargs.get('atol', 0.5)
    #
    #     cluster_pairs = []
    #     cpair_extend = cluster_pairs.extend
    #     isclose = np.isclose
    #     where = np.where
    #     dsplit = np.dsplit
    #     array_ = np.array
    #     vsplit = np.vsplit
    #
    #     i = 0
    #     for cents_chunk_1 in cents_chunks_1:
    #         inds_cents_chunk_1 = range(i, i + len(cents_chunk_1))
    #         dmat_coords_1 = relvdm1._all_ifgs[cents_chunk_1]
    #         sh_1 = dmat_coords_1.shape[0]
    #         dmat_coords_1 = dmat_coords_1.reshape(sh_1 * num_atoms_ifg1, 3)
    #         j = 0
    #         for cents_chunk_2 in cents_chunks_2:
    #             inds_cents_chunk_2 = range(j, j + len(cents_chunk_2))
    #             dmat_coords_2 = relvdm2._all_ifgs[cents_chunk_2]
    #             sh_2 = dmat_coords_2.shape[0]
    #             dmat_coords_2 = dmat_coords_2.reshape(sh_2 * num_atoms_ifg2, 3)
    #             dmats_array = cdist(dmat_coords_1, dmat_coords_2)
    #             dmats = array_(dsplit(array_(vsplit(dmats_array, sh_1)), sh_2))
    #             wh = where(isclose(dmats, lig_dmat, atol=atol).all(axis=2).all(axis=2))
    #             cpair_extend(((inds_cents_chunk_1[ind2], relvdm1.name), (inds_cents_chunk_2[ind1], relvdm2.name))
    #                          for ind1, ind2 in zip(wh[0], wh[1]))
    #             j += len(cents_chunk_2)
    #         i += len(cents_chunk_1)
    #     return cluster_pairs

    def find_cluster_pairs(self, relvdm1, relvdm2, **kwargs):

        print('finding cluster pairs for ' + relvdm1.name + ', ' + relvdm2.name)
        num_atoms_ifg1 = kwargs.get('num_atoms_ifg1', len(str(relvdm1.ifg_dict.values()).split()))
        num_atoms_ifg2 = kwargs.get('num_atoms_ifg2', len(str(relvdm2.ifg_dict.values()).split()))
        cents_chunks_1 = kwargs.get('cents_chunks_1', list(self.chunks(relvdm1._cents, 2000)))
        cents_chunks_2 = kwargs.get('cents_chunks_2', list(self.chunks(relvdm2._cents, 2000)))
        lig_dmat = kwargs.get('lig_dmat', self.make_lig_dmat(relvdm1, relvdm2))
        # atol = kwargs.get('atol', 0.5)
        atol1 = kwargs.get('atol1', 0.6)
        atol2 = kwargs.get('atol2', 1.5)
        pc = kwargs.get('pc', 0.5)

        cluster_pairs = []
        cpair_extend = cluster_pairs.extend
        isclose = np.isclose
        where = np.where
        dsplit = np.dsplit
        array_ = np.array
        vsplit = np.vsplit

        i = 0
        for cents_chunk_1 in cents_chunks_1:
            inds_cents_chunk_1 = range(i, i + len(cents_chunk_1))
            dmat_coords_1 = relvdm1._all_ifgs[cents_chunk_1]
            sh_1 = dmat_coords_1.shape[0]
            dmat_coords_1 = dmat_coords_1.reshape(sh_1 * num_atoms_ifg1, 3)
            j = 0
            for cents_chunk_2 in cents_chunks_2:
                inds_cents_chunk_2 = range(j, j + len(cents_chunk_2))
                dmat_coords_2 = relvdm2._all_ifgs[cents_chunk_2]
                sh_2 = dmat_coords_2.shape[0]
                dmat_coords_2 = dmat_coords_2.reshape(sh_2 * num_atoms_ifg2, 3)
                dmats_array = cdist(dmat_coords_1, dmat_coords_2)
                dmats = array_(dsplit(array_(vsplit(dmats_array, sh_1)), sh_2))
                submats = np.abs(dmats-lig_dmat)
                t1 = (submats < atol2).all(axis=2).all(axis=2)
                t2 = ((submats < atol1).sum(axis=2).sum(axis=2)/(num_atoms_ifg1*num_atoms_ifg2)) >= pc
                wh = where(t1*t2)
                # wh = where(isclose(dmats, lig_dmat, atol=atol).all(axis=2).all(axis=2))
                cpair_extend(((inds_cents_chunk_1[ind2], relvdm1.name), (inds_cents_chunk_2[ind1], relvdm2.name))
                             for ind1, ind2 in zip(wh[0], wh[1]))
                j += len(cents_chunk_2)
            i += len(cents_chunk_1)
        return cluster_pairs

    def find_ligand_cst_cliques(self, relvdms, **kwargs):
        g = nx.Graph()
        # maxlen = len(list(itertools.combinations(range(len(relvdms)), 2)))
        maxlen = kwargs.get('maxlen', len(list(itertools.combinations(range(len(relvdms)), 2))))
        for rel_vdm1, rel_vdm2 in itertools.combinations(relvdms, 2):
            cluster_pairs = self.find_cluster_pairs(rel_vdm1, rel_vdm2, **kwargs)
            g.add_edges_from(cluster_pairs)
        cliques = [cl for cl in nx.find_cliques(g) if len(cl) >= maxlen]
        self.cliques = cliques

    # def make_cluster_densities(self, relvdms):
    #     sorted_cliques = np.array([sorted(clique, key=lambda x: x[1]) for clique in self.cliques])  #sort by relvdm.name
    #     self.density = collections.defaultdict(dict)
    #     for relvdm in relvdms:
    #         hs_ind = [i for i, tup in enumerate(sorted_cliques[0]) if tup[1] == relvdm.name][0]
    #         hs_nums = np.array(np.unique(sorted_cliques[:, hs_ind, 0]), dtype=int)
    #         func = functools.partial(self.set_density_cluster, relvdm)
    #         list(map(func, hs_nums))

    def make_cluster_densities(self, relvdms):
        # sorted_cliques = np.array([sorted(clique, key=lambda x: x[1]) for clique in self.cliques])  #sort by relvdm.name
        self.density = collections.defaultdict(dict)
        for relvdm in relvdms:
            # hs_ind = [i for i, tup in enumerate(sorted_cliques[0]) if tup[1] == relvdm.name][0]
            hs_nums = [cl[0] for cliq in self.cliques for cl in cliq if cl[1] == relvdm.name]
            func = functools.partial(self.set_density_cluster, relvdm)
            list(map(func, hs_nums))

    # @staticmethod
    # def set_density_cluster(self, relvdm, hs_num):
    #     ifgs = relvdm._all_ifgs[relvdm._mems[hs_num]]
    #     self.density[relvdm.name][hs_num] = KernelDensity(kernel='gaussian', bandwidth=2.0).fit(ifgs)
    def set_density_cluster(self, relvdm, hs_num):
        ifgs = relvdm._all_ifgs[relvdm._mems[hs_num]]
        self.density[relvdm.name][hs_num] = stats.gaussian_kde(ifgs.T)

    # def fit_lig_to_cluster_density(self, relvdms, lig, lig_coords_1, q, js, cliqs):
    #     # print('len js', len(js))
    #     # print('len cliqs', len(cliqs))
    #     poses = []
    #     for j, cliq in zip(js, cliqs):
    #         print('fitting clique ' + str(j))
    #         # print(cliq)
    #         density_1 = self.density[cliq[1][1]][cliq[1][0]]
    #         partial_score_fit_1 = functools.partial(score_fit_1, [lig_coords_1, density_1])  # only fits carboxamide
    #         min1_ = minimize(partial_score_fit_1, np.zeros(6), method='Nelder-Mead', options={'maxiter': 2500})
    #         # minimizer_kwargs = {'method': 'Nelder-Mead', 'options': {'maxiter': 2500}}
    #         # min_ = basinhopping(partial_score_fit, np.zeros(6), niter=3, minimizer_kwargs=minimizer_kwargs)
    #         lig_copy = lig.copy()
    #         lig_min_coords_all = rigid_body(lig_copy.getCoords(), min1_.x)
    #         lig_copy.setCoords(lig_min_coords_all)
    #         lig_ifg_coords = {}
    #         for relvdm in relvdms:
    #             lig_ifg_coords[relvdm.name] = self.get_lig_coords(relvdm, lig=lig_copy).reshape(1, -1)
    #         lig_coords = [lig_ifg_coords[tup[1]] for tup in cliq]
    #         densities = [self.density[tup[1]][tup[0]] for tup in cliq]
    #         partial_score_fit = functools.partial(score_fit, list(zip(lig_coords, densities)))
    #         min_ = minimize(partial_score_fit, np.zeros(6), method='Nelder-Mead', options={'maxiter': 2500})
    #         # lig_min_coords_heavy = rigid_body(lig_copy.select('not element H D').getCoords(), min_.x)
    #         lig_min_coords_NO = rigid_body(lig_copy.select('element N O').getCoords(), min_.x)
    #         lig_min_coords_C = rigid_body(lig_copy.select('element C').getCoords(), min_.x)
    #         # if self.bb_clash_test(lig_min_coords_heavy):
    #         if self.bb_clash_test(lig_min_coords_NO, lig_min_coords_C):
    #             lig_min_coords_all = rigid_body(lig_copy.getCoords(), min_.x)
    #             lig_copy.setCoords(lig_min_coords_all)
    #             pose = Pose()
    #             for tup in cliq:
    #                 pose.lig_ifg_coords[tup[1]] = rigid_body_highdim(lig_ifg_coords[tup[1]], min_.x)
    #             pose.ligand = lig_copy
    #             pose.hotspot_clique = cliq
    #             pose.optim_result = min_
    #             pose.number = j
    #             poses.append(pose)
    #     # q.put(poses)
    #     # self.poses.append(poses)
    #     q.append(poses)

    def fit_lig_to_cluster_density_mult_proc(self, relvdms, lig, order, q, js, cliqs):
        # print('len js', len(js))
        # print('len cliqs', len(cliqs))
        poses = []
        lig_ifg_coords = {}
        ligcoords_start = lig.getCoords()
        for relvdm in relvdms:
            lig_ifg_coords[relvdm.name] = self.get_lig_coords(relvdm).reshape(1, -1)
        for j, cliq in zip(js, cliqs):
            print('fitting clique ' + str(j))
            # print(cliq)
            cliq_names = [c[1] for c in cliq]
            for name in order:
                if name in cliq_names:
                    lig_coords_1 = lig_ifg_coords[name]
                    ind = [i for i, c in enumerate(cliq) if c[1] == name][0]
                    density_1 = self.density[name][cliq[ind][0]]
                    break
            partial_score_fit_1 = functools.partial(score_fit_1, [lig_coords_1, density_1])  # only fits carboxamide
            min1_ = minimize(partial_score_fit_1, np.zeros(6), method='Nelder-Mead', options={'maxiter': 2500})
            lig_min_coords_all = rigid_body(ligcoords_start, min1_.x)
            lig_copy = lig.copy()
            lig_copy.setCoords(lig_min_coords_all)
            lig_ifg_coords_2 = {}
            for relvdm in relvdms:
                lig_ifg_coords_2[relvdm.name] = self.get_lig_coords(relvdm, lig=lig_copy).reshape(1, -1)
            lig_coords = [lig_ifg_coords_2[tup[1]] for tup in cliq]
            densities = [self.density[tup[1]][tup[0]] for tup in cliq]
            partial_score_fit = functools.partial(score_fit, list(zip(lig_coords, densities)))
            min_ = minimize(partial_score_fit, np.zeros(6), method='Nelder-Mead', options={'maxiter': 2500})
            # lig_min_coords_heavy = rigid_body(lig_copy.select('not element H D').getCoords(), min_.x)
            lig_min_coords_NO = rigid_body(lig_copy.select('element N O').getCoords(), min_.x)
            lig_min_coords_C = rigid_body(lig_copy.select('element C').getCoords(), min_.x)
            # if self.bb_clash_test(lig_min_coords_heavy):
            if self.bb_clash_test(lig_min_coords_NO, lig_min_coords_C):
                lig_min_coords_all = rigid_body(lig_copy.getCoords(), min_.x)
                lig_copy.setCoords(lig_min_coords_all)
                pose = Pose()
                for tup in cliq:
                    pose.lig_ifg_coords[tup[1]] = rigid_body_highdim(lig_ifg_coords_2[tup[1]], min_.x)
                pose.ligand = lig_copy
                pose.hotspot_clique = cliq
                pose.optim_result = min_
                pose.number = j
                poses.append(pose)
                # q.put(poses)
                # self.poses.append(poses)
        q.append(poses)

    def fit_lig_to_cluster_density_one_proc(self, relvdms, lig, order, js, cliqs):
        # print('len js', len(js))
        # print('len cliqs', len(cliqs))
        poses = []
        lig_ifg_coords = {}
        ligcoords_start = lig.getCoords()
        for relvdm in relvdms:
            lig_ifg_coords[relvdm.name] = self.get_lig_coords(relvdm).reshape(1, -1)
        for j, cliq in zip(js, cliqs):
            print('fitting clique ' + str(j))
            # print(cliq)
            cliq_names = [c[1] for c in cliq]
            for name in order:
                if name in cliq_names:
                    lig_coords_1 = lig_ifg_coords[name]
                    ind = [i for i, c in enumerate(cliq) if c[1] == name][0]
                    density_1 = self.density[name][cliq[ind][0]]
                    break
            partial_score_fit_1 = functools.partial(score_fit_1, [lig_coords_1, density_1])  # only fits carboxamide
            min1_ = minimize(partial_score_fit_1, np.zeros(6), method='Nelder-Mead', options={'maxiter': 2500})
            lig_min_coords_all = rigid_body(ligcoords_start, min1_.x)
            lig_copy = lig.copy()
            lig_copy.setCoords(lig_min_coords_all)
            lig_ifg_coords_2 = {}
            for relvdm in relvdms:
                lig_ifg_coords_2[relvdm.name] = self.get_lig_coords(relvdm, lig=lig_copy).reshape(1, -1)
            lig_coords = [lig_ifg_coords_2[tup[1]] for tup in cliq]
            densities = [self.density[tup[1]][tup[0]] for tup in cliq]
            partial_score_fit = functools.partial(score_fit, list(zip(lig_coords, densities)))
            min_ = minimize(partial_score_fit, np.zeros(6), method='Nelder-Mead', options={'maxiter': 2500})
            # lig_min_coords_heavy = rigid_body(lig_copy.select('not element H D').getCoords(), min_.x)
            lig_min_coords_NO = rigid_body(lig_copy.select('element N O').getCoords(), min_.x)
            lig_min_coords_C = rigid_body(lig_copy.select('element C').getCoords(), min_.x)
            # if self.bb_clash_test(lig_min_coords_heavy):
            if self.bb_clash_test(lig_min_coords_NO, lig_min_coords_C):
                lig_min_coords_all = rigid_body(lig_copy.getCoords(), min_.x)
                lig_copy.setCoords(lig_min_coords_all)
                pose = Pose()
                for tup in cliq:
                    pose.lig_ifg_coords[tup[1]] = rigid_body_highdim(lig_ifg_coords_2[tup[1]], min_.x)
                pose.ligand = lig_copy
                pose.hotspot_clique = cliq
                pose.optim_result = min_
                pose.number = j
                poses.append(pose)
        # q.put(poses)
        # self.poses.append(poses)
        # q.append(poses)
        return poses


    def fit_lig_to_cluster_densities(self, relvdms, order=None, num_proc = None):
        if num_proc:

            # num_proc = num_proc

            for lig in self.ligand_conformers:

                partial_fit_dens = functools.partial(self.fit_lig_to_cluster_density_mult_proc, relvdms, lig, order)
                ch_size = int(np.ceil(len(self.cliques) / num_proc))
                clique_ch = list(self.chunks(self.cliques, ch_size))
                ind_ch = list(self.chunks(list(range(len(self.cliques))), ch_size))

                man = Manager()
                q = man.list()

                def proc_q(q, arg1, arg2):
                    return Process(target=partial_fit_dens, args=(q, arg1, arg2))

                partial_fit_dens_q = functools.partial(proc_q, q)
                processes = [partial_fit_dens_q(i, s) for i, s in zip(ind_ch, clique_ch)]
                for p in processes:
                    p.start()
                for p in processes:
                    p.join()
                self.poses = [pose for pose_ch in list(q) for pose in pose_ch]

        else:
            for lig in self.ligand_conformers:
                js = range(len(self.cliques))
                self.poses = self.fit_lig_to_cluster_density_one_proc(relvdms, lig, order, js, self.cliques)



    # def fit_lig_to_cluster_densities(self, relvdms):
    #     lig_coords_1 = self.get_lig_coords(relvdms[1]).reshape(1, -1)
    #
    #     for lig in self.ligand_conformers:
    #         i=1
    #         sorted_cliques = [sorted(clique, key=lambda x: x[1]) for clique in self.cliques]
    #         for j, cliq in enumerate(sorted_cliques):
    #             print('fitting clique ' + str(i))
    #             i += 1
    #             density_1 = self.density[cliq[1][1]][cliq[1][0]]
    #             partial_score_fit_1 = functools.partial(score_fit_1, [lig_coords_1, density_1]) # only fits carboxamide
    #             min1_ = minimize(partial_score_fit_1, np.zeros(6), method='Nelder-Mead', options={'maxiter': 2500})
    #             # minimizer_kwargs = {'method': 'Nelder-Mead', 'options': {'maxiter': 2500}}
    #             # min_ = basinhopping(partial_score_fit, np.zeros(6), niter=3, minimizer_kwargs=minimizer_kwargs)
    #             lig_copy = lig.copy()
    #             lig_min_coords_all = rigid_body(lig_copy.getCoords(), min1_.x)
    #             lig_copy.setCoords(lig_min_coords_all)
    #             lig_ifg_coords = {}
    #             for relvdm in relvdms:
    #                 lig_ifg_coords[relvdm.name] = self.get_lig_coords(relvdm, lig=lig_copy).reshape(1, -1)
    #             lig_coords = [lig_ifg_coords[tup[1]] for tup in cliq]
    #             densities = [self.density[tup[1]][tup[0]] for tup in cliq]
    #             partial_score_fit = functools.partial(score_fit, list(zip(lig_coords, densities)))
    #             min_ = minimize(partial_score_fit, np.zeros(6), method='Nelder-Mead', options={'maxiter': 2500})
    #             # lig_min_coords_heavy = rigid_body(lig_copy.select('not element H D').getCoords(), min_.x)
    #             lig_min_coords_NO = rigid_body(lig_copy.select('element N O').getCoords(), min_.x)
    #             lig_min_coords_C = rigid_body(lig_copy.select('element C').getCoords(), min_.x)
    #             # if self.bb_clash_test(lig_min_coords_heavy):
    #             if self.bb_clash_test(lig_min_coords_NO, lig_min_coords_C):
    #                 lig_min_coords_all = rigid_body(lig_copy.getCoords(), min_.x)
    #                 lig_copy.setCoords(lig_min_coords_all)
    #                 pose = Pose()
    #                 for tup in cliq:
    #                     pose.lig_ifg_coords[tup[1]] = rigid_body_highdim(lig_ifg_coords[tup[1]], min_.x)
    #                 pose.ligand = lig_copy
    #                 pose.hotspot_clique = cliq
    #                 pose.optim_result = min_
    #                 pose.number = j
    #                 self.poses.append(pose)
    #                 # print(self.bb_clash_test(lig_min_coords_heavy))