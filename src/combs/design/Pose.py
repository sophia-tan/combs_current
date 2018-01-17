__all__ = ['Pose']
import networkx as nx
from sklearn.neighbors import NearestNeighbors
import numpy as np
import functools
import prody as pr
import os
import random
import itertools
import collections
import operator
from ..analysis.cluster import cluster_adj_mat
from scipy.spatial.distance import cdist
import traceback
import pyrosetta as py
from ..apps.functions import writePDBStream

class Pose:
    def __init__(self):
        self.ligand = None
        self.lig_ifg_coords = {}
        self.number = None
        self.hotspot_clique = None
        self.optim_result = None
        self.subgraphs = {}
        self.member_score = None
        self.num_members = {}
        # self.subgraphs_by_site = collections.defaultdict(dict)
        self.supergraph = nx.Graph()
        # self.site_nodes = {}
        self.site_graph = {}
        self.sub_poses = []
        self.total_residues = 0
        self.num_residues = 0
        self.num_satisfied_iFGs = 0

    # def get_around_lig_fg(self, sample, rmsd_cutoff=1):
    #     for tup in self.hotspot_clique:
    #         hs_num, relvdm_key = tup
    #         graph = sample.rel_vdms[relvdm_key].hotspot_subgraphs[hs_num]
    #         nbrs_node_ifgs = NearestNeighbors(metric='euclidean', radius=rmsd_cutoff)
    #         nbrs_node_ifgs.fit(list(nx.get_node_attributes(graph, 'ifg').values()))
    #         lig_ifg_coords = self.lig_ifg_coords[relvdm_key]
    #         adj_mat = nbrs_node_ifgs.radius_neighbors_graph(lig_ifg_coords)
    #         inds = adj_mat.nonzero()[1]
    #         nodes = np.array(list(nx.get_node_attributes(graph, 'ifg').keys()))[inds]
    #         nodes_list = [tuple(node) for node in nodes]
    #         self.subgraphs[relvdm_key] = graph.subgraph(nodes_list)

    def get_around_lig_fg(self, sample):
        self.mems = {}
        self.score_breakdown = {}
        for name in self.lig_ifg_coords.keys():
            self.mems[name] = sample._nbrs[name].radius_neighbors(self.lig_ifg_coords[name], return_distance=False)[0]
            self.score_breakdown[name] = len(self.mems[name])

    # def get_nonclashing(self, clash_cutoff=2.5):
    #     lig_heavy_coords = self.ligand.select('not element H D').getCoords()
    #     for relvdm_key, graph in self.subgraphs.items():
    #         nodes_set = set(graph.nodes())
    #         sc_coords_dict = nx.get_node_attributes(graph, 'sc')
    #         if sc_coords_dict:
    #             sc_coords = np.array(functools.reduce(lambda x, y: np.vstack((x, y)), sc_coords_dict.values()))
    #             sc_nodes = np.array(functools.reduce(lambda x, y: np.vstack((x, y)),
    #                                                  (coords.shape[0] * [node] for node, coords in sc_coords_dict.items())))
    #             nbrs_sc = NearestNeighbors(metric='euclidean', radius=clash_cutoff, n_neighbors=1)
    #             nbrs_sc.fit(lig_heavy_coords)
    #             dist, ind = nbrs_sc.kneighbors(sc_coords)
    #             bad_nodes = sc_nodes[(dist < nbrs_sc.radius).reshape(-1)]
    #             bad_nodes_set = {tuple(node) for node in bad_nodes}
    #             good_nodes = nodes_set - bad_nodes_set
    #             self.subgraphs[relvdm_key] = graph.subgraph(good_nodes)

    def get_nonclashing(self, sample, cutoff_NO_NO=2.5, cutoff_NO_C=2.8, cutoff_C_C=3.2):
        self.mems_prune = {}

        if self.ligand.select('element O N') is not None:
            lig_ON_coords = self.ligand.select('element O N').getCoords()
        else:
            lig_ON_coords = np.array([])
        if self.ligand.select('(not element O N) and (not element H D)') is not None:
            lig_C_coords = self.ligand.select('(not element O N) and (not element H D)').getCoords()
        else:
            lig_C_coords = np.array([])

        for name, mems in self.mems.items():
            sc_ON_coords = []
            sc_CS_coords = []
            mems_ON = []
            mems_CS = []
            bad_mems = []
            if mems.shape[0] > 0:
                for m in mems:
                    type = sample.relvdms[name]._all_type[m]
                    resn = sample.relvdms[name]._all_resn[m]
                    resnum_chid = tuple(sample.relvdms[name]._all_resnum_chid[m])
                    index = sample.relvdms[name]._all_indices[m]
                    if type == 'SC':
                        try:
                            coords = sample.relvdms[name].rel_vdm_sc_ON_coords[resnum_chid][type][resn][index]
                            sc_ON_coords.append(coords)
                            mems_ON.append([[m]] * coords.shape[0])
                        except:
                            print('exception ON')
                            traceback.print_exc()

                        try:
                            coords = sample.relvdms[name].rel_vdm_sc_CS_coords[resnum_chid][type][resn][index]
                            sc_CS_coords.append(coords)
                            mems_CS.append([[m]] * coords.shape[0])
                        except:
                            print('exception CS')
                            traceback.print_exc()

                bad_mems_ON = np.array([])
                if sc_ON_coords:
                    if lig_ON_coords.shape[0] > 0:
                        sc_ON_coords = np.array(functools.reduce(lambda x, y: np.vstack((x, y)), sc_ON_coords))
                        mems_ON = np.array(functools.reduce(lambda x, y: np.vstack((x, y)), mems_ON))
                        nbrs_sc_ON_ON = NearestNeighbors(metric='euclidean', n_neighbors=1)
                        nbrs_sc_ON_ON.fit(lig_ON_coords)
                        dist_ON_ON, ind_ON_ON = nbrs_sc_ON_ON.kneighbors(sc_ON_coords)
                        fail_ON_ON = (dist_ON_ON < cutoff_NO_NO).reshape(-1)
                    else:
                        fail_ON_ON = False

                    if lig_C_coords.shape[0] > 0:
                        nbrs_sc_C_ON = NearestNeighbors(metric='euclidean', n_neighbors=1)
                        nbrs_sc_C_ON.fit(lig_C_coords)
                        dist_C_ON, ind_ON_ON = nbrs_sc_C_ON.kneighbors(sc_ON_coords)
                        fail_C_ON = (dist_C_ON < cutoff_NO_C).reshape(-1)
                    else:
                        fail_C_ON = False

                    bad_mems_ON = mems_ON[fail_ON_ON | fail_C_ON]

                bad_mems_CS = np.array([])
                if sc_CS_coords:
                    if lig_ON_coords.shape[0] > 0:
                        sc_CS_coords = np.array(functools.reduce(lambda x, y: np.vstack((x, y)), sc_CS_coords))
                        mems_CS = np.array(functools.reduce(lambda x, y: np.vstack((x, y)), mems_CS))
                        nbrs_sc_ON_C = NearestNeighbors(metric='euclidean', n_neighbors=1)
                        nbrs_sc_ON_C.fit(lig_ON_coords)
                        dist_ON_C, ind_ON_C = nbrs_sc_ON_C.kneighbors(sc_CS_coords)
                        fail_ON_C = (dist_ON_C < cutoff_NO_C).reshape(-1)
                    else:
                        fail_ON_C = False

                    if lig_C_coords.shape[0] > 0:
                        nbrs_sc_C_C = NearestNeighbors(metric='euclidean', n_neighbors=1)
                        nbrs_sc_C_C.fit(lig_C_coords)
                        dist_C_C, ind_C_C = nbrs_sc_C_C.kneighbors(sc_CS_coords)
                        fail_C_C = (dist_C_C < cutoff_C_C).reshape(-1)
                    else:
                        fail_C_C = False

                    bad_mems_CS = mems_CS[fail_ON_C | fail_C_C]

                if bad_mems_CS.shape[0] > 0 and bad_mems_ON.shape[0] > 0:
                    bad_mems = np.union1d(bad_mems_CS, bad_mems_ON)
                elif bad_mems_CS.shape[0] > 0:
                    bad_mems = bad_mems_CS
                elif bad_mems_ON.shape[0] > 0:
                    bad_mems = bad_mems_ON

            self.mems_prune[name] = np.setdiff1d(mems, bad_mems)

    def set_member_score(self):
        self.member_score = np.prod([num_mem + 1 for num_mem in self.num_members.values()])

#compare clashes across all scs of each hotspot.
    def _print_node(self, node, sample, relvdm_key, outdir):
        resnum_chid = node[0]
        type_ = self.subgraphs[relvdm_key].node[node]['type']
        typestr = type_
        resn = self.subgraphs[relvdm_key].node[node]['resn']
        vdm_tags = self.subgraphs[relvdm_key].node[node]['vdm_tags']

        if typestr == 'PHI_PSI':
            typestr = 'PHI_PSI/' + sample.rel_vdms[relvdm_key].rel_vdm_phipsi_bin[resnum_chid]
        pdbpath = sample.rel_vdms[relvdm_key].rel_vdm_path + typestr + '/pdbs/' + resn + '/'
        filename = 'iFG_' + str(vdm_tags[0]) + '_vdM_' + str(vdm_tags[1]) + '_iFlip_' \
                   + str(vdm_tags[2]) + '_' + sample.rel_vdms[relvdm_key].name + '_' + 'oriented.pdb.gz'
        pdb = pr.parsePDB(pdbpath + filename)
        old_coords = pdb.getCoords()
        new_coords = \
            np.dot((old_coords - sample.rel_vdms[relvdm_key]._rois_rot_trans[resnum_chid][type_][resn][1]),
                   sample.rel_vdms[relvdm_key]._rois_rot_trans[resnum_chid][type_][resn][0]) \
            + sample.rel_vdms[relvdm_key]._rois_rot_trans[resnum_chid][type_][resn][2]
        pdb.setCoords(new_coords)
        newfile_path = outdir + 'hotspots/' + str(self.number) + '/'
        newfile_name = sample.rel_vdms[relvdm_key].name + '_hotspot_' + str(self.number) + '_' \
                       + ''.join(str(x) for x in resnum_chid) + '_' + type_ + '_' + filename
        try:
            os.makedirs(newfile_path)
        except:
            pass
        pr.writePDB(newfile_path + newfile_name, pdb)

    def print_graph_pdbs(self, sample, relvdm_key, outdir, randmembers=None):
        """prints pdbs of the top number of hotspots (number) to the output directory (outdir)."""
        if outdir[-1] != '/':
            outdir += '/'
        if randmembers is None:
            for node in self.subgraphs[relvdm_key]:
                self._print_node(node, sample, relvdm_key, outdir)
        elif randmembers is not None:
            for node in random.sample(list(self.subgraphs[relvdm_key]), randmembers):
                self._print_node(node, sample, relvdm_key, outdir)

    def _zip_gen(self):
        for bs_pos, graph in self.site_graph.items():
            yield zip(len(graph) * [bs_pos], graph.nodes())

    @staticmethod
    def setup_score_fn(path_to_params=None):
        if path_to_params:
            py.init('-extra_res_fa ' + path_to_params  + ' -in:Ntermini XA -in:Ctermini XA')
        else:
            py.init(' -in:Ntermini XA -in:Ctermini XA')
        fa_atr = py.rosetta.core.scoring.ScoreType.fa_atr
        fa_rep = py.rosetta.core.scoring.ScoreType.fa_rep
        fa_elec_sc_sc = py.rosetta.core.scoring.ScoreType.fa_elec_sc_sc
        fa_elec_bb_sc = py.rosetta.core.scoring.ScoreType.fa_elec_bb_sc
        hbond_sc = py.rosetta.core.scoring.ScoreType.hbond_sc
        hbond_bb_sc = py.rosetta.core.scoring.ScoreType.hbond_bb_sc
        ch_hbond_sc_sc = py.rosetta.core.scoring.ScoreType.ch_bond_sc_sc
        ch_hbond_bb_sc = py.rosetta.core.scoring.ScoreType.ch_bond_bb_sc
        fa_dun = py.rosetta.core.scoring.ScoreType.fa_dun

        scorefxn_self = py.ScoreFunction()
        scorefxn_self.set_weight(fa_atr, 0.1)
        scorefxn_self.set_weight(fa_rep, 0.44)
        scorefxn_self.set_weight(hbond_sc, 20)
        scorefxn_self.set_weight(ch_hbond_sc_sc, 1)
        scorefxn_self.set_weight(fa_elec_sc_sc, 0.1)
        scorefxn_self.set_weight(fa_dun, 0.025)
        
        scorefxn_pair = py.ScoreFunction()
        scorefxn_pair.set_weight(fa_atr, 0.1)
        scorefxn_pair.set_weight(fa_rep, 0.44)
        scorefxn_pair.set_weight(hbond_sc, 10)
        scorefxn_pair.set_weight(ch_hbond_sc_sc, 1)
        scorefxn_pair.set_weight(fa_elec_sc_sc, 0.1)

        scorefxn_pose = py.ScoreFunction()
        scorefxn_pose.set_weight(fa_atr, 0.1)
        scorefxn_pose.set_weight(fa_rep, 0.44)
        scorefxn_pose.set_weight(hbond_sc, 15)
        scorefxn_pose.set_weight(hbond_bb_sc, 15)
        scorefxn_pose.set_weight(ch_hbond_sc_sc, 1)
        scorefxn_pose.set_weight(ch_hbond_bb_sc, 1)
        scorefxn_pose.set_weight(fa_elec_sc_sc, 0.5)
        scorefxn_pose.set_weight(fa_elec_bb_sc, 0.5)
        scorefxn_pose.set_weight(fa_dun, 0.025)

        return scorefxn_self, scorefxn_pair, scorefxn_pose

    def make_sub_poses(self, sample, path_to_params=None):
        # group by res position
        # group by sc vs bb
        # group sc by res name
        # cluster scs to 0.3A rmsd and take centroids. This is the sc set
        # take rep from bb.
        # enumerate all pose sets and remove poses with clashes.
        # Do you do less than the total iFGs? For now, not interested.
        # append good poses to self.sub_poses.  Should look like: {relvdm name: mem_indices}

        scorefxn_self, scorefxn_pair, scorefxn_pose = self.setup_score_fn(path_to_params=path_to_params)

        class Vividict(dict):
            def __missing__(self, key):
                value = self[key] = type(self)()
                return value

        mem_info = Vividict()
        sc_mem_info = collections.defaultdict(dict)
        sc_mem_coords = collections.defaultdict(dict)
        for relvdm_name, mems in self.mems_prune.items():
            for mem in mems:
                resnum_chid = tuple(sample.relvdms[relvdm_name]._all_resnum_chid[mem])
                type_ = sample.relvdms[relvdm_name]._all_type[mem]
                resn = sample.relvdms[relvdm_name]._all_resn[mem]
                index = sample.relvdms[relvdm_name]._all_indices[mem]
                vdm_tags = sample.relvdms[relvdm_name]._all_vdm_tags[mem]
                mem_info[resnum_chid][type_][resn][relvdm_name][index] = vdm_tags
                if type_ == 'SC':
                    try:
                        sc_mem_info[resnum_chid][resn]
                    except:
                        sc_mem_info[resnum_chid][resn] = []
                        sc_mem_coords[resnum_chid][resn] = []
                    try:
                        coor = sample.relvdms[relvdm_name].rel_vdm_sc_coords[resnum_chid][type_][resn][index].flatten()
                        sc_mem_coords[resnum_chid][resn].append(coor)
                        sc_mem_info[resnum_chid][resn].append((index, relvdm_name, mem))
                    except:
                        print('exception SC coords')
                        traceback.print_exc()

        sc_centroids = Vividict()
        bs_sc_cens = collections.defaultdict(dict)
        # print(1)
        print('sc_mem_coords=', sc_mem_coords)
        keys = sc_mem_coords.keys()
        for resnum_chid in keys:
            # print(resnum_chid)
            # print(sc_mem_coords[resnum_chid].keys())
            if sc_mem_coords[resnum_chid].keys():
                for resn in sc_mem_coords[resnum_chid].keys():
                    # print(resn, resnum_chid)
                    if not sc_mem_coords[resnum_chid][resn]:
                        sc_mem_coords[resnum_chid].pop(resn)
            else:
                sc_mem_coords.pop(resnum_chid)

        for resnum_chid in sc_mem_coords.keys():
            for resn in sc_mem_coords[resnum_chid].keys():
                if len(sc_mem_coords[resnum_chid][resn]) > 2:
                    num_atoms = len(sc_mem_coords[resnum_chid][resn][0])
                    nbrs = NearestNeighbors(metric='euclidean', radius=np.sqrt(num_atoms) * 0.3)
                    nbrs.fit(sc_mem_coords[resnum_chid][resn])
                    adj_mat = nbrs.radius_neighbors_graph(sc_mem_coords[resnum_chid][resn])
                    m, c = cluster_adj_mat(adj_mat)
                    sc_centroids[resnum_chid][resn] = [sc_mem_info[resnum_chid][resn][cen] for cen in c]
                else:
                    sc_centroids[resnum_chid][resn] = [sc_mem_info[resnum_chid][resn][0]]
            print('sc_centroids=', sc_centroids)
            bs_sc_cens[resnum_chid]['mems'] = [c[2]
                                               for resn in sc_centroids[resnum_chid].keys()
                                               for c in sc_centroids[resnum_chid][resn]]
            bs_sc_cens[resnum_chid]['relvdm_name'] = [c[1]
                                                      for resn in sc_centroids[resnum_chid].keys()
                                                      for c in sc_centroids[resnum_chid][resn]]
            i = 1
            rnch = ''.join(str(x) for x in resnum_chid)
            for mem, relvdm_name in zip(bs_sc_cens[resnum_chid]['mems'], bs_sc_cens[resnum_chid]['relvdm_name']):
                poi = sample.poi.copy()
                sample.relvdms[relvdm_name].print_member(poi, mem, 'temp', 'temp_pdbs/temp_' + rnch, self.ligand, counter=i)
                i += 1

        Es = {}
        for resnum_chid in sc_mem_coords.keys():
            rnch = ''.join(str(x) for x in resnum_chid)
            dir_ = 'temp_pdbs/temp_' + rnch
            # for f in [file for file in os.listdir(dir_) if file[0] != '.']:
            #     os.system('/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL ' + dir_ + '/' + f + ' -c -d "save '
            #               + dir_ + '/' + f.split('.')[0] + '_py.pdb"')
            # sorted_pdbs = sorted([file for file in os.listdir(dir_) if (file[0] != '.' and file[-6:-4] == 'py')])
            sorted_pdbs = sorted([file for file in os.listdir(dir_) if file[0] != '.'])
            E_ar = np.zeros(len(sorted_pdbs))
            for i, pdb in enumerate(sorted_pdbs):
               pose = py.pose_from_pdb(dir_ + '/' + pdb)
               E_ar[i] = scorefxn_self(pose)
            Es[resnum_chid] = E_ar

        for resnum_chid1, resnum_chid2 in itertools.combinations(sc_mem_coords.keys(), 2):
            rnch1 = ''.join(str(x) for x in resnum_chid1)
            rnch2 = ''.join(str(x) for x in resnum_chid2)
            i = 0
            for mem1, relvdm_name1 in zip(bs_sc_cens[resnum_chid1]['mems'], bs_sc_cens[resnum_chid1]['relvdm_name']):
                poi = sample.poi.copy()
                bb1, pdb1, filename1 = sample.relvdms[relvdm_name1]._print_pair(poi, mem1)
                j = 0
                for mem2, relvdm_name2 in zip(bs_sc_cens[resnum_chid2]['mems'], bs_sc_cens[resnum_chid2]['relvdm_name']):
                    poi = sample.poi.copy()
                    bb2, pdb2, filename2 = sample.relvdms[relvdm_name2]._print_pair(poi, mem2)
                    if np.abs(resnum_chid1[0] - resnum_chid2[0]) > 1 and resnum_chid1[1] == resnum_chid2[1]:
                        self.print_pair(bb1, bb2, pdb1, pdb2, filename1, filename2, 'temp_pair', 'temp_pdbs_pair/temp_' + rnch1 + '_' + rnch2,
                                   counter1=i, counter2=j, chains=2)
                    else:
                        self.print_pair(bb1, bb2, pdb1, pdb2, filename1, filename2, 'temp_pair',
                                        'temp_pdbs_pair/temp_' + rnch1 + '_' + rnch2,
                                        counter1=i, counter2=j, chains=1)
                    j += 1
                i += 1

        Ep = collections.defaultdict(dict)
        for resnum_chid1, resnum_chid2 in itertools.combinations(sc_mem_coords.keys(), 2):
            rnch1 = ''.join(str(x) for x in resnum_chid1)
            rnch2 = ''.join(str(x) for x in resnum_chid2)
            dir_ = 'temp_pdbs_pair/temp_' + rnch1 + '_' + rnch2
            # for f in [file for file in os.listdir(dir_) if file[0] != '.']:
            #     os.system('/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL ' + dir_ + '/' + f + ' -c -d "save '
            #               + dir_ + '/' + f.split('.')[0] + '_py.pdb"')
            # sorted_pdbs = sorted([file for file in os.listdir(dir_) if (file[0] != '.' and file[-6:-4] == 'py')])
            sorted_pdbs = sorted([file for file in os.listdir(dir_) if file[0] != '.'])
            E_ar = np.zeros((len(sorted_pdbs),len(sorted_pdbs)))
            for pdb in sorted_pdbs:
                pose = py.pose_from_pdb(dir_ + '/' + pdb)
                i = int(pdb.split('_')[0])
                j = int(pdb.split('_')[1])
                E_ar[i][j] = scorefxn_pair(pose)
            Ep[resnum_chid1][resnum_chid2] = E_ar
            Ep[resnum_chid2][resnum_chid1] = E_ar.T

        opt_ens = []
        opt_ress = []
        del_dicts = []
        for i in reversed(range(1, len(sc_mem_coords.keys()) + 1)):
            for combo in itertools.combinations(sc_mem_coords.keys(), i):
                ep = collections.defaultdict(dict)
                es = {}
                print(combo)
                for resnum_chid in combo:
                    es[resnum_chid] = Es[resnum_chid]
                if len(combo) == 1:
                    opt_ens.append(np.min(es[resnum_chid]))
                    argmin = np.argmin(es[resnum_chid])
                    opt_ress.append([(argmin, resnum_chid)])
                    del_dicts.append({})
                else:
                    for resnum_chid1, resnum_chid2 in itertools.combinations(combo, 2):
                        ep[resnum_chid1][resnum_chid2] = Ep[resnum_chid1][resnum_chid2]
                        ep[resnum_chid2][resnum_chid1] = Ep[resnum_chid2][resnum_chid1]
                    es, ep, del_dict = self.dee(es, ep)
                    opt_en, opt_res = self.find_opt_en(es, ep)
                    opt_ens.append(opt_en)
                    opt_ress.append(opt_res)
                    del_dicts.append(del_dict)



        # Es, Ep, del_dict = self.dee(Es, Ep)
        # print("Es, Ep, del_dict=", Es, Ep, del_dict)
        #
        # opt_en, opt_res = self.find_opt_en(Es, Ep)
        # print('opt_en, opt_res=', opt_en, opt_res)
        #
        # print("bs_sc_cens=", bs_sc_cens)
        # bs_sc_cens_dee = collections.defaultdict(dict)
        self.opt_energy = np.min(opt_ens)
        ind_opt = np.argmin(opt_ens)

        for key, val in del_dicts[ind_opt].items():
            bs_sc_cens[key]['mems'] = [m for i, m in enumerate(bs_sc_cens[key]['mems']) if i not in set(val)]
            bs_sc_cens[key]['relvdm_name'] = [m for i, m in enumerate(bs_sc_cens[key]['relvdm_name']) if i not in set(val)]

        print("bs_sc_cens=", bs_sc_cens)
        # self.opt_energy = opt_en
        self.opt_residues = {}
        for tup in opt_ress[ind_opt]:
            rnchid = tup[1]
            ind = tup[0]
            self.opt_residues[rnchid] = (bs_sc_cens[rnchid]['mems'][ind], bs_sc_cens[rnchid]['relvdm_name'][ind])

    @staticmethod
    def print_pair(bb1, bb2, pdb1, pdb2, filename1, filename2, label, outdir, counter1=0, counter2=0, chains=2):
        """prints pdbs of the top number of hotspots (number) to the output directory (outdir)."""
        if outdir[-1] != '/':
            outdir += '/'

        # pdb1sel = pdb1.select('sidechain and chain X and resnum 10')
        # pdb2sel = pdb2.select('sidechain and chain X and resnum 10')

        rn1 = pdb1.select('name CA and chain X and resnum 10').getResnames()[0]
        rn2 = pdb2.select('name CA and chain X and resnum 10').getResnames()[0]

        bb1.setChids('X')
        bb1.setResnums(10)
        num_ind = len(bb1.getIndices())
        start = 1
        finish = start + num_ind
        bb1.setBetas(list(range(start, finish)))
        bb1.setResnames(rn1)

        sc1 = pdb1.select('sidechain and chain X and resnum 10')
        num_ind = len(sc1.getIndices())
        start = finish
        finish = start + num_ind
        sc1.setBetas(list(range(start, finish)))

        if chains == 2:
            bb2.setChids('W')
            bb2.setResnums(10)
        elif chains == 1:
            bb2.setChids('X')
            bb2.setResnums(11)
        num_ind = len(bb2.getIndices())
        start = finish
        finish = start + num_ind
        bb2.setBetas(list(range(start, finish)))
        bb2.setResnames(rn2)

        sc2 = pdb2.select('sidechain and chain X and resnum 10')
        num_ind = len(sc2.getIndices())
        start = finish
        finish = start + num_ind
        sc2.setBetas(list(range(start, finish)))
        if chains == 2:
            sc2.setChids('W')
        elif chains == 1:
            sc2.setResnums(11)

        newfile_path = outdir
        newfile_name = str(counter1) + '_' + str(counter2) + '_mem_' + label + '_' + filename1 \
                       + '_' + filename2
        try:
            os.makedirs(newfile_path)
        except:
            pass
        with open(newfile_path + newfile_name, 'w') as outfile:
            writePDBStream(outfile, bb1)
            writePDBStream(outfile, sc1)
            writePDBStream(outfile, bb2)
            writePDBStream(outfile, sc2)

        #make array of self energies per residue
        #make array of pair energies per residue pairs
        #perform DEE and brute force optimization

        # sites = list(range(4))
        # Es = dict((s, -np.random.rand(20, 1)) for s in sites)
        # Ep = collections.defaultdict(dict)
        # for s1, s2 in itertools.permutations(sites, 2):
        #     Ep[s1][s2] = -np.random.rand(20, 20) + 0.5 * np.random.rand(20, 20)
        # Es, Ep, del_dict = dee(Es, Ep)

    # @staticmethod
    # def dee(Es, Ep):
    #     Es = Es.copy()
    #     Ep = Ep.copy()
    #     del_dict = {}
    #     for s in Es.keys():
    #         to_del = []
    #         for i in range(len(Es[s])):
    #             for j in set(range(len(Es[s]))) - {i}:
    #                 cond = Es[s][i] - Es[s][j] + min(
    #                     Ep[s][k][i, w] - Ep[s][k][j, w] for k in set(Ep[s].keys()) - {s} for w in
    #                     range(Ep[s][k].shape[1]))
    #                 if cond > 0:
    #                     to_del.append(i)
    #                     break
    #         if to_del:
    #             Es[s] = np.delete(Es[s], to_del, 0)
    #             for key in Ep[s].keys():
    #                 Ep[s][key] = np.delete(Ep[s][key], to_del, 0)
    #             for key in set(Ep.keys()) - {s}:
    #                 Ep[key][s] = np.delete(Ep[key][s], to_del, 1)
    #             del_dict[s] = to_del
    #     return Es, Ep, del_dict

    @staticmethod
    def dee(Es, Ep):
        Es = Es.copy()
        Ep = Ep.copy()
        del_dict = {}
        for s in Es.keys():
            to_del = []
            for i in range(len(Es[s])):
                for j in set(range(len(Es[s]))) - {i}:
                    pair_ens = []
                    for k in set(Ep[s].keys()) - {s}:
                        if len(Ep[s][k].shape) > 1:
                            for w in range(Ep[s][k].shape[1]):
                                pair_ens.append(Ep[s][k][i, w] - Ep[s][k][j, w])
                        else:
                            for w in range(Ep[s][k].shape[0]):
                                pair_ens.append(Ep[s][k][i, w] - Ep[s][k][j, w])
                    cond = Es[s][i] - Es[s][j] + min(pair_ens)
                    if cond > 0:
                        to_del.append(i)
                        break
            if to_del:
                Es[s] = np.delete(Es[s], to_del, 0)
                for key in Ep[s].keys():
                    Ep[s][key] = np.delete(Ep[s][key], to_del, 0)
                for key in set(Ep.keys()) - {s}:
                    Ep[key][s] = np.delete(Ep[key][s], to_del, 1)
                del_dict[s] = to_del
        return Es, Ep, del_dict

    @staticmethod
    def find_opt_en(Es, Ep):
        l = len(Es.keys())
        opt_res = list(zip([-1] * l, Es.keys()))
        opt_en = np.inf
        ranges = map(range, map(len, Es.values()))
        for inds in itertools.product(*ranges):
            tuples = list(zip(inds, Es.keys()))
            c_en = np.sum(Es[key][ind] for ind, key in tuples)
            c_en += np.sum(Ep[t1[1]][t2[1]][t1[0], t2[0]] for t1, t2 in itertools.combinations(tuples, 2))
            if c_en < opt_en:
                opt_en = c_en
                opt_res = tuples
        return opt_en, opt_res

    def write_pose(self, sample, outdir):
        if outdir[-1] != '/':
            outdir += '/'

        with open(outdir + 'pose.pdb', 'w') as outfile:
            finish = 1
            for i, (resnum_chid, sidechain) in enumerate(sorted(self.pose_dict.items())):
                if sidechain:
                    poi = sample.poi.copy()
                    bb = poi.select('(backbone or name H) and resnum ' + str(resnum_chid[0]) + ' and chain ' + resnum_chid[1])
                    mem = sidechain[0]
                    name = sidechain[1]
                    relvdm = sample.relvdms[name]
                    pdb, resn = relvdm.parse_member(mem)
                    sc = pdb.select('chain X and resnum 10 and sidechain')

                    bb.setChids('A')
                    bb.setResnums(i + 1)
                    num_ind = len(bb.getIndices())
                    start = finish
                    finish = start + num_ind
                    bb.setBetas(list(range(start, finish)))
                    bb.setResnames(resn)

                    sc.setChids('A')
                    sc.setResnums(i + 1)
                    num_ind = len(sc.getIndices())
                    start = finish
                    finish = start + num_ind
                    sc.setBetas(list(range(start, finish)))

                    writePDBStream(outfile, bb)
                    writePDBStream(outfile, sc)
                else:
                    poi = sample.poi.copy()
                    bb = poi.select('(backbone or name H) and resnum ' + str(resnum_chid[0]) + ' chain ' + resnum_chid[1])
                    bb.setChids('A')
                    bb.setResnums(i + 1)
                    num_ind = len(bb.getIndices())
                    start = finish
                    finish = start + num_ind
                    bb.setBetas(list(range(start, finish)))
                    bb.setResnames('GLY')
                    writePDBStream(outfile, bb)

            num_ind = len(self.ligand.select('all').getIndices())
            start = finish
            finish = start + num_ind
            self.ligand.setBetas(list(range(start, finish)))
            writePDBStream(outfile, self.ligand)



    # def make_sub_poses(self, sample):
    #     #group by res position
    #     #group by sc vs bb
    #     #group sc by res name
    #     #cluster scs to 0.3A rmsd and take centroids. This is the sc set
    #     #take rep from bb.
    #     #enumerate all pose sets and remove poses with clashes.
    #     #Do you do less than the total iFGs? For now, not interested.
    #     #append good poses to self.sub_poses.  Should look like: {relvdm name: mem_indices}
    #
    #     class Vividict(dict):
    #         def __missing__(self, key):
    #             value = self[key] = type(self)()
    #             return value
    #
    #
    #     mem_info = Vividict()
    #     sc_mem_info = collections.defaultdict(dict)
    #     sc_mem_coords = collections.defaultdict(dict)
    #     for relvdm_name, mems in self.mems_prune.items():
    #         for mem in mems:
    #             resnum_chid = tuple(sample.relvdms[relvdm_name]._all_resnum_chid[mem])
    #             type_ = sample.relvdms[relvdm_name]._all_type[mem]
    #             resn = sample.relvdms[relvdm_name]._all_resn[mem]
    #             index = sample.relvdms[relvdm_name]._all_indices[mem]
    #             vdm_tags = sample.relvdms[relvdm_name]._all_vdm_tags[mem]
    #             mem_info[resnum_chid][type_][resn][relvdm_name][index] = vdm_tags
    #             if type_ == 'SC':
    #                 try:
    #                     sc_mem_info[resnum_chid][resn]
    #                 except:
    #                     sc_mem_info[resnum_chid][resn] = []
    #                     sc_mem_coords[resnum_chid][resn] = []
    #                 try:
    #                     sc_mem_coords[resnum_chid][resn].append(
    #                         sample.relvdms[relvdm_name].rel_vdm_sc_coords[resnum_chid][type_][resn][index].flatten())
    #                     sc_mem_info[resnum_chid][resn].append((index, relvdm_name))
    #                 except:
    #                     print('exception SC coords')
    #                     traceback.print_exc()
    #
    #     sc_centroids = Vividict()
    #     sc_centroids_coords = Vividict()
    #     sc_centroids_info = Vividict()
    #     sc_centroids_ON_coords_flat = {}
    #     sc_centroids_ON_info_flat = {}
    #     sc_centroids_ON_resn_flat = {}
    #     sc_centroids_CS_coords_flat = {}
    #     sc_centroids_CS_info_flat = {}
    #     sc_centroids_CS_resn_flat = {}
    #     for resnum_chid in sc_mem_coords.keys():
    #         for resn in sc_mem_coords[resnum_chid].keys():
    #             num_atoms = len(sc_mem_coords[resnum_chid][resn][0])
    #             nbrs = NearestNeighbors(metric='euclidean', radius=np.sqrt(num_atoms)*0.3)
    #             nbrs.fit(sc_mem_coords[resnum_chid][resn])
    #             adj_mat = nbrs.radius_neighbors_graph(sc_mem_coords[resnum_chid][resn])
    #             m, c = cluster_adj_mat(adj_mat)
    #             sc_centroids[resnum_chid][resn] = [sc_mem_info[resnum_chid][resn][cen] for cen in c]
    #             try:
    #                 sc_centroids_coords[resnum_chid]['CS'][resn] = \
    #                     [sample.relvdms[relvdm_name].rel_vdm_sc_CS_coords[resnum_chid]['SC'][resn][index]
    #                      for index, relvdm_name in sc_centroids[resnum_chid][resn]]
    #                 sc_centroids_info[resnum_chid]['CS'][resn] = sc_centroids[resnum_chid][resn]
    #             except:
    #                 print('exception CS coords')
    #                 traceback.print_exc()
    #             try:
    #                 sc_centroids_coords[resnum_chid]['ON'][resn] = \
    #                     [sample.relvdms[relvdm_name].rel_vdm_sc_ON_coords[resnum_chid]['SC'][resn][index]
    #                      for index, relvdm_name in sc_centroids[resnum_chid][resn]]
    #                 sc_centroids_info[resnum_chid]['ON'][resn] = sc_centroids[resnum_chid][resn]
    #             except:
    #                 print('exception ON coords')
    #                 traceback.print_exc()
    #
    #         sc_centroids_ON_coords_flat[resnum_chid] = [sc_centroids_coords[resnum_chid]['ON'][resn] for resn in
    #                                                     sc_centroids_coords[resnum_chid]['ON'].keys()]
    #         sc_centroids_ON_coords_flat[resnum_chid] = [item for list_ in sc_centroids_ON_coords_flat[resnum_chid]
    #                                                     for item in list_]
    #         sc_centroids_ON_info_flat[resnum_chid] = [sc_centroids_info[resnum_chid]['ON'][resn] for resn in
    #                                                     sc_centroids_info[resnum_chid]['ON'].keys()]
    #         sc_centroids_ON_info_flat[resnum_chid] = [item for list_ in sc_centroids_ON_info_flat[resnum_chid]
    #                                                     for item in list_]
    #         sc_centroids_ON_resn_flat[resnum_chid] = [[resn]*len(sc_centroids_info[resnum_chid]['ON'][resn]) for resn in
    #                                                   sc_centroids_info[resnum_chid]['ON'].keys()]
    #         sc_centroids_ON_resn_flat[resnum_chid] = [item for list_ in sc_centroids_ON_resn_flat[resnum_chid]
    #                                                   for item in list_]
    #
    #         sc_centroids_CS_coords_flat[resnum_chid] = [sc_centroids_coords[resnum_chid]['CS'][resn] for resn in
    #                                                     sc_centroids_coords[resnum_chid]['CS'].keys()]
    #         sc_centroids_CS_coords_flat[resnum_chid] = [item for list_ in sc_centroids_CS_coords_flat[resnum_chid]
    #                                                     for item in list_]
    #         sc_centroids_CS_info_flat[resnum_chid] = [sc_centroids_info[resnum_chid]['CS'][resn] for resn in
    #                                                   sc_centroids_info[resnum_chid]['CS'].keys()]
    #         sc_centroids_CS_info_flat[resnum_chid] = [item for list_ in sc_centroids_CS_info_flat[resnum_chid]
    #                                                   for item in list_]
    #         sc_centroids_CS_resn_flat[resnum_chid] = [[resn] * len(sc_centroids_info[resnum_chid]['CS'][resn]) for resn
    #                                                   in
    #                                                   sc_centroids_info[resnum_chid]['CS'].keys()]
    #         sc_centroids_CS_resn_flat[resnum_chid] = [item for list_ in sc_centroids_CS_resn_flat[resnum_chid]
    #                                                   for item in list_]
    #
    #         for coords1, coords2 in itertools.combinations(sc_centroids_ON_coords_flat.values(), 2):
    #
    #         for coords1, coords2 in itertools.combinations(sc_centroids_CS_coords_flat.values(), 2):
    #
    #         for coords1, coords2 in itertools.product(sc_centroids_ON_coords_flat.values(), sc_centroids_CS_coords_flat.values()):
    #
    #
    #
    #
    #
    #
    #
    #
    #         # pairwise clash calcs, add non-clashing to graph
    #         # find cliques
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #     mem_info = []
    #     for relvdm_name, mems in self.mems_prune.items():
    #         for mem in mems:
    #             resnum_chid = tuple(sample.relvdms[relvdm_name]._all_resnum_chid[mem])
    #             type_ = sample.relvdms[relvdm_name]._all_type[mem]
    #             resn =  sample.relvdms[relvdm_name]._all_resn[mem]
    #             index = sample.relvdms[relvdm_name]._all_indices[mem]
    #             mem_info.append([resnum_chid, type_, resn, index])
    #
    #     sorted_ = sorted(mem_info, key=lambda x: x[0])
    #     gr_ = itertools.groupby(sorted_, key=lambda x: x[0])
    #     gr_by_resnum_chid = []
    #     keys_resnum_chid = []
    #     for k, g in gr_:
    #         gr_by_resnum_chid.append(list(g))
    #         keys_resnum_chid.append(k)
    #
    #     gr_by_type = []
    #     keys_type = []
    #     for group in gr_by_resnum_chid:
    #         sorted_ = sorted(group, key=lambda x: x[1])
    #         gr_ = itertools.groupby(sorted_, key=lambda x: x[1])
    #         gr_int = []
    #         key_int = []
    #         for k, g in gr_:
    #             gr_int.append(list(g))
    #             key_int.append(k)
    #         gr_by_type.append(gr_int)
    #         keys_type.append(key_int)
    #
    #     gr_by_resn = []
    #     keys_resn = []
    #
    #
    #
    #
    #
    #     resnum_chids = collections.defaultdict(list)
    #     for relvdm_name, mems in self.mems_prune.items():
    #         for mem in mems:
    #             resnum_chid = tuple(sample.relvdms[relvdm_name]._all_resnum_chid[mem])
    #             resnum_chids[resnum_chid].append((relvdm_name, mem))
    #
    #     sc_dict = collections.defaultdict(list)
    #     bb_dict = collections.defaultdict(list)
    #     for resnum_chid, mems_tup in resnum_chids.items():
    #         relvdm_name = mems_tup[0]
    #         mem = mems_tup[1]
    #         if sample.relvdms[relvdm_name]._all_type[mem]
    #             sc_dict[resnum_chid]['SC'] =
    #         else:
    #             bb_dict[resnum_chid]['BB'] =
    #
    #
    #             resnum_chid = tuple(sample.relvdms[relvdm_name]._all_resnum_chid[mem])
    #             resnum_chids[resnum_chid].append((relvdm_name, mem))
    #
    #
    #
    #
    #
    #
    # def make_sub_poses(self):
    #     i = 0
    #     for name, graph in self.subgraphs.items():
    #         nx.set_node_attributes(graph, 'relvdm_name', name)
    #         bs_pos_dict = {node: node[0] for node in graph}
    #         nx.set_node_attributes(graph, 'bs_pos', bs_pos_dict)
    #         nx.relabel_nodes(graph, dict(zip(graph.nodes(), range(i, i + len(graph)))), copy=False)
    #         i += len(graph)
    #
    #     for graph in self.subgraphs.values():
    #         self.supergraph.add_nodes_from(graph.nodes(data=True))
    #
    #     site_dict = nx.get_node_attributes(self.supergraph, 'bs_pos')
    #     for bs_pos, nodes in itertools.groupby(sorted(site_dict.items(),
    #                                                   key=operator.itemgetter(1)), key=lambda x: x[1]):
    #         # self.site_graph[bs_pos] = {}
    #         nodes = [node[0] for node in nodes]
    #         site_graph_bs = self.supergraph.subgraph(nodes)
    #         type_dict = nx.get_node_attributes(site_graph_bs, 'type')
    #         type_set = set(type_dict.values())
    #         node_list = []
    #         if 'PHI_PSI' in type_set:
    #             node_list.append((key for key, val in type_dict.items() if val == 'PHI_PSI').__next__())
    #             # self.site_graph[bs_pos]['PHI_PSI'] = site_graph_bs.subgraph(node)
    #
    #         if 'C_O' in type_set:
    #             node_list.append((key for key, val in type_dict.items() if val == 'C_O').__next__())
    #             # self.site_graph[bs_pos]['C_O'] = site_graph_bs.subgraph(node)
    #
    #         if 'N_CA' in type_set:
    #             node_list.append((key for key, val in type_dict.items() if val == 'N_CA').__next__())
    #             # self.site_graph[bs_pos]['N_CA'] = site_graph_bs.subgraph(node)
    #
    #         if 'SC' in type_set:
    #             # node_list.extend(list(nx.get_node_attributes(site_graph_bs, 'sc').keys()))
    #             sc_graph = site_graph_bs.subgraph(list(nx.get_node_attributes(site_graph_bs, 'sc').keys()))
    #             sc_dict = nx.get_node_attributes(sc_graph, 'sc')
    #             for resn, resn_nodes in itertools.groupby(sorted(nx.get_node_attributes(sc_graph, 'resn').items(),
    #                                                              key=operator.itemgetter(1)), key=lambda x: x[1]):
    #                 resn_nodes = [resn_node[0] for resn_node in resn_nodes]
    #
    #                 if len(resn_nodes) > 1:
    #                     sc_coords = [sc_dict[resn_node].flatten() for resn_node in resn_nodes]
    #                     nbrs = NearestNeighbors(metric='euclidean', radius=0.6)
    #                     nbrs.fit(sc_coords)
    #                     adj_mat = nbrs.radius_neighbors_graph(sc_coords)
    #                     members, centroids = cluster_adj_mat(adj_mat)
    #                     good_nodes = [resn_nodes[cent] for cent in centroids]
    #                     node_list.extend(good_nodes)
    #                 else:
    #                     node_list.extend(resn_nodes)
    #
    #                     #     self.site_graph[bs_pos]['SC'][resn] = sc_graph.subgraph(list(resn_nodes))
    #
    #         self.site_graph[bs_pos] = self.supergraph.subgraph(node_list)
    #
    #     # gen = (itertools.zip_longest(bs_pos, graph.nodes(), fillvalue=bs_pos)
    #     #        for bs_pos, graph in self.site_graph.items())
    #
    #     for node_set in itertools.product(*self._zip_gen()):
    #         sc_coords = []
    #         test = True
    #         for each in node_set:
    #             bs_pos = each[0]
    #             node = each[1]
    #             if self.site_graph[bs_pos].node[node]['type'] == 'SC':
    #                 sc_coords.append(self.site_graph[bs_pos].node[node]['sc'])
    #         if sc_coords:
    #             if len(sc_coords) > 1:
    #                 dists = []
    #                 list(dists.extend(cdist(combo[0], combo[1]).flatten())
    #                      for combo in itertools.combinations(sc_coords, 2))
    #                 if (np.array(dists) < 2.5).any():
    #                     test = False
    #                 else:
    #                     test = True
    #         if test:
    #             subpose_graph = nx.Graph()
    #             subby_graphs = [self.site_graph[each2[0]].subgraph(each2[1]) for each2 in node_set]
    #             list(subpose_graph.add_nodes_from(gra.nodes(data=True)) for gra in subby_graphs)
    #             self.sub_poses.append(subpose_graph)


    # def make_sub_poses(self):
    #     i = 0
    #     for name, graph in self.subgraphs.items():
    #         nx.set_node_attributes(graph, 'relvdm_name', name)
    #         bs_pos_dict = {node: node[0] for node in graph}
    #         nx.set_node_attributes(graph, 'bs_pos', bs_pos_dict)
    #         nx.relabel_nodes(graph, dict(zip(graph.nodes(), range(i, i + len(graph)))), copy=False)
    #         i += len(graph)
    #
    #     for graph in self.subgraphs.values():
    #         self.supergraph.add_nodes_from(graph.nodes(data=True))
    #
    #     site_dict = nx.get_node_attributes(self.supergraph, 'bs_pos')
    #     for bs_pos, nodes in itertools.groupby(sorted(site_dict.items(),
    #                                            key=operator.itemgetter(1)), key=lambda x: x[1]):
    #         # self.site_graph[bs_pos] = {}
    #         nodes = [node[0] for node in nodes]
    #         site_graph_bs = self.supergraph.subgraph(nodes)
    #         type_dict = nx.get_node_attributes(site_graph_bs, 'type')
    #         type_set = set(type_dict.values())
    #         node_list = []
    #         if 'PHI_PSI' in type_set:
    #             node_list.append((key for key, val in type_dict.items() if val == 'PHI_PSI').__next__())
    #             # self.site_graph[bs_pos]['PHI_PSI'] = site_graph_bs.subgraph(node)
    #
    #         if 'C_O' in type_set:
    #             node_list.append((key for key, val in type_dict.items() if val == 'C_O').__next__())
    #             # self.site_graph[bs_pos]['C_O'] = site_graph_bs.subgraph(node)
    #
    #         if 'N_CA' in type_set:
    #             node_list.append((key for key, val in type_dict.items() if val == 'N_CA').__next__())
    #             # self.site_graph[bs_pos]['N_CA'] = site_graph_bs.subgraph(node)
    #
    #         if 'SC' in type_set:
    #             # node_list.extend(list(nx.get_node_attributes(site_graph_bs, 'sc').keys()))
    #             sc_graph = site_graph_bs.subgraph(list(nx.get_node_attributes(site_graph_bs, 'sc').keys()))
    #             sc_dict = nx.get_node_attributes(sc_graph, 'sc')
    #             for resn, resn_nodes in itertools.groupby(sorted(nx.get_node_attributes(sc_graph, 'resn').items(),
    #                                                       key=operator.itemgetter(1)), key=lambda x: x[1]):
    #                 resn_nodes = [resn_node[0] for resn_node in resn_nodes]
    #
    #                 if len(resn_nodes) > 1:
    #                     sc_coords = [sc_dict[resn_node].flatten() for resn_node in resn_nodes]
    #                     nbrs = NearestNeighbors(metric='euclidean', radius=0.6)
    #                     nbrs.fit(sc_coords)
    #                     adj_mat = nbrs.radius_neighbors_graph(sc_coords)
    #                     members, centroids = cluster_adj_mat(adj_mat)
    #                     good_nodes = [resn_nodes[cent] for cent in centroids]
    #                     node_list.extend(good_nodes)
    #                 else:
    #                     node_list.extend(resn_nodes)
    #
    #             #     self.site_graph[bs_pos]['SC'][resn] = sc_graph.subgraph(list(resn_nodes))
    #
    #         self.site_graph[bs_pos] = self.supergraph.subgraph(node_list)
    #
    #     # gen = (itertools.zip_longest(bs_pos, graph.nodes(), fillvalue=bs_pos)
    #     #        for bs_pos, graph in self.site_graph.items())
    #
    #     for node_set in itertools.product(*self._zip_gen()):
    #         sc_coords = []
    #         test = True
    #         for each in node_set:
    #             bs_pos = each[0]
    #             node = each[1]
    #             if self.site_graph[bs_pos].node[node]['type'] == 'SC':
    #                 sc_coords.append(self.site_graph[bs_pos].node[node]['sc'])
    #         if sc_coords:
    #             if len(sc_coords) > 1:
    #                 dists = []
    #                 list(dists.extend(cdist(combo[0], combo[1]).flatten())
    #                      for combo in itertools.combinations(sc_coords, 2))
    #                 if (np.array(dists) < 2.5).any():
    #                     test = False
    #                 else:
    #                     test = True
    #         if test:
    #             subpose_graph = nx.Graph()
    #             subby_graphs = [self.site_graph[each2[0]].subgraph(each2[1]) for each2 in node_set]
    #             list(subpose_graph.add_nodes_from(gra.nodes(data=True)) for gra in subby_graphs)
    #             self.sub_poses.append(subpose_graph)

    def print_sub_pose(self, subpose, sample, outdir):
        # for chid_num, node in enumerate(self.sub_poses[subpose]):
        if outdir[-1] != '/':
            outdir += '/'

        try:
            os.makedirs(outdir)
        except:
            pass

        with open(outdir + 'subpose_' + str(subpose) + '.pdb.gz', 'w') as outfile:
            i = 30
            for node, data in self.sub_poses[subpose].nodes(data=True):
                relvdm_key = data['relvdm_name']
                resnum_chid = data['bs_pos']
                type_ = data['type']
                typestr = type_
                resn = data['resn']
                vdm_tags = data['vdm_tags']

                if typestr == 'PHI_PSI':
                    typestr = 'PHI_PSI/' + sample.rel_vdms[relvdm_key].rel_vdm_phipsi_bin[resnum_chid]
                pdbpath = sample.rel_vdms[relvdm_key].rel_vdm_path + typestr + '/pdbs/' + resn + '/'
                filename = 'iFG_' + str(vdm_tags[0]) + '_vdM_' + str(vdm_tags[1]) + '_iFlip_' \
                           + str(vdm_tags[2]) + '_' + sample.rel_vdms[relvdm_key].name + '_' + 'oriented.pdb.gz'
                pdb = pr.parsePDB(pdbpath + filename)
                old_coords = pdb.getCoords()
                new_coords = \
                    np.dot((old_coords - sample.rel_vdms[relvdm_key]._rois_rot_trans[resnum_chid][type_][resn][1]),
                           sample.rel_vdms[relvdm_key]._rois_rot_trans[resnum_chid][type_][resn][0]) \
                    + sample.rel_vdms[relvdm_key]._rois_rot_trans[resnum_chid][type_][resn][2]
                pdb.setCoords(new_coords)
                pdb.select('chain X and resnum 10').setResnums(i)
                # pdb.select('chain Y').setChids(i+1)
                if type_ in ['N_CA', 'PHI_PSI']:
                    x_sel = pdb.select('chain X and resnum ' + str(i) + ' and name N CA')
                    # y_sel = pdb.select('chain ' + str(i+1))
                    # pdb_sel = x_sel | y_sel
                if type_ == 'C_O':
                    x_sel = pdb.select('chain X and resnum ' + str(i) + ' and name CA C O')
                    # y_sel = pdb.select('chain ' + str(i+1))
                    # pdb_sel = x_sel | y_sel
                if type_ == 'SC':
                    x_sel = pdb.select('chain X and resnum ' + str(i) + ' and (sidechain or name CA)')
                    # y_sel = pdb.select('chain ' + str(i+1))
                    # pdb_sel = x_sel | y_sel
                i += 1
                # newfile_path = outdir + 'hotspots/' + str(self.number) + '/'
                # newfile_name = sample.rel_vdms[relvdm_key].name + '_subpose_' + str(subpose) + '_' \
                #                + ''.join(str(x) for x in resnum_chid) + '_' + type_ + '_' + filename
                # pdb.select('chain X').setChid()
                pr.writePDBStream(outfile, x_sel)

    def print_sub_poses(self, sample, outdir):
        for i in range(len(self.sub_poses)):
            self.print_sub_pose(i, sample, outdir)


            # need sc_coords
    #     nbrs = NearestNeighbors(metric='euclidean', radius=rmsd_cutoff)
    #     nbrs.fit(sc_coords)
    #     adj_mat = nbrs.radius_neighbors_graph(sc_coords)
    #     members, centroids = cluster_adj_mat(adj_mat)








            # for relvdm_key, graph in self.subgraphs.items():
            #     for k, g in itertools.groupby(sorted(self.subgraphs[relvdm_key].nodes()), lambda x: x[0]):
            #         # node_groups_by_site[k][relvdm_key] = graph.subgraph(list(g))
            #         self.subgraphs_by_site[k][relvdm_key] = {}
            #         resn_dict = nx.get_node_attributes(graph.subgraph(list(g)), 'resn')
            #         for k2, g2 in itertools.groupby(sorted(resn_dict.items(), key=operator.itemgetter(1)),
            #                                         key = lambda x: x[1]):
            #             self.subgraphs_by_site[k][relvdm_key][k2] = graph.subgraph(list(g2))


                # unique_node_groups.append(k)





            # itertools.product()
            #
            # res_position --> ifgtype --> resn --> sc_coords --> rmsd
            #
            # 1. for ifgtype, get all nodes that have sc_coords, group by residue position and name.
            # (if rmsd < 0.5, it's the same'). if the sc is centroid of cluster, mark it as special: either
            # centroid of 1 ifgtype or centroid of multiple ifgtypes. next, get clashing pairs between all
            # vdms at each res position. These will be removed from consideration.



