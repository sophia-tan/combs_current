import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import prody as pr
import pickle
import numpy as np
import functools
from sklearn.neighbors import NearestNeighbors
import os
import copy
import collections


pdb_path = '/Users/npolizzi/Projects/combs/src/runs/apixaban/'
sample = combs.Sample()
poi_path = '/Users/npolizzi/Projects/combs/src/runs/apixaban/20171126/input_helices/'

# site_path = '/Users/npolizzi/Downloads/scPDB_2013/sc_gt80_le300_MWgt400_sitele35res_sites_pdbs/'

sample.poi = pr.parsePDB(poi_path + '1.pdb', chain='A')
# site = pr.parsePDB(site_path + 'helix_template.pdb')

# site_CA = site.select('name CA')

# sels = []
# for ca in site_CA:
#     sels.append(sample.poi.select('exwithin .1 of coords', coords=ca.getCoords()))
# resis = [6, 7, 9, 10, 12, 13, 14, 16, 17, 19, 20, 21]
resis = [6, 7, 9, 10, 13, 14, 17 20, 21]  # removed residues that can place apix outside bundle
sample.bs_residues = list(zip(resis, ['A']*len(resis)))
print(sample.bs_residues)
sample.set_rois()
sample.set_rois_phipsi()
sample.set_poi_clash_coords()
sample.set_roi_bb_coords()
sample.ligand_conformers = [pr.parsePDB(pdb_path + 'APX_0001.pdb')]

relvdm_backbone = combs.rel_vandermer.RelVandermer('backbone')
relvdm_backbone.rel_vdm_path = '/Users/npolizzi/Projects/combs/results/backbone/rel_vdms_hbond_carbonyl/20171022/'

no_charge = set(combs.constants.one_letter_code.keys()) - {'LYS', 'ARG', 'GLU', 'ASP', 'MSE', 'ANY',
                                                           'VAL', 'LEU', 'PHE', 'ILE', 'PRO', 'MET', 'CYS'}
relvdm_backbone.load_rel_vdms_pickle(sample, subset=no_charge)
relvdm_backbone.set_rel_vdm_bb_coords()
relvdm_backbone.set_rois_rot_trans(sample)
relvdm_backbone.set_rel_vdm_tags(sample)
print('moving vdMs')
relvdm_backbone.move_rel_vdms(sample)
print('removing clashing vdMs')
relvdm_backbone.remove_clash(sample)
relvdm_backbone.reshape_ifgs()
all_ifgs = functools.reduce(lambda a, b: np.vstack((a, b)),
                            [val for val in relvdm_backbone._ifgs.values()])
print('finding hotspots carbonyl')
nbrs = NearestNeighbors(metric='euclidean', radius=1.0, algorithm='kd_tree')
nbrs.fit(all_ifgs)
adj_mat = nbrs.radius_neighbors_graph(all_ifgs)
print('clustering...')

mems, cents = combs.analysis.cluster.greedy_cluster_pc(adj_mat, pc=0.9)

all_resnum_chid = functools.reduce(lambda a, b: np.vstack((a, b)),
                                   [np.array([tuple(key)]*len(val), dtype=object)
                                    for key, val in relvdm_backbone._vdm_tags.items() if len(val) > 0])
all_vdm_tags = functools.reduce(lambda a, b: np.vstack((a, b)),
                                [val for val in relvdm_backbone._vdm_tags.values() if len(val) > 0])
all_resn = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_backbone._resn.values() if len(val) > 0])
all_type = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_backbone._type.values() if len(val) > 0])
all_indices = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_backbone._indices.values() if len(val) > 0])
relvdm_backbone._all_vdm_tags = all_vdm_tags
relvdm_backbone._all_type = all_type
relvdm_backbone._all_resn = all_resn
relvdm_backbone._all_resnum_chid = all_resnum_chid
relvdm_backbone._all_ifgs = all_ifgs
print('carbonyl ifgs: ', len(all_ifgs))
relvdm_backbone._all_indices = all_indices



relvdm_backbone.rel_vdms_pickle = None
relvdm_backbone.rel_vdm_ifg_coords = None
relvdm_backbone.rel_vdm_bb_coords = None
relvdm_backbone.rel_vdm_sc_coords = None
relvdm_backbone.rel_vdm_tags = None
relvdm_backbone.rel_vdm_resnames = None
relvdm_backbone._ifgs = None
relvdm_backbone._resn = None
relvdm_backbone._type = None
relvdm_backbone._vdm_tags = None
relvdm_backbone._indices = None
relvdm_backbone.ifg_dict = {'ANY': 'C O'}
relvdm_backbone._mems = mems
relvdm_backbone._cents = cents
relvdm_backbone.lig_ifg_correspondence = 'C8 O3'

outdir = '/Users/npolizzi/Projects/combs/results/apixaban/20171126/output_helices/'
try:
    os.makedirs(outdir)
except:
    pass

# with open(outdir + 'relvdm_backbone.pickle', 'wb') as outfile:
#     pickle.dump(relvdm_backbone, outfile)


relvdm_backbone2 = copy.deepcopy(relvdm_backbone)
relvdm_backbone2.lig_ifg_correspondence = 'C19 O2'
relvdm_backbone2.name = 'backbone2'

relvdm_backbone3 = copy.deepcopy(relvdm_backbone)
relvdm_backbone3.lig_ifg_correspondence = 'N1 N6'
relvdm_backbone3.name = 'backbone3'

relvdm_carboxamide = combs.rel_vandermer.RelVandermer('carboxamide')
relvdm_carboxamide.rel_vdm_path = '/Users/npolizzi/Projects/combs/results/carboxamide/rel_vdms_hbond/20171023/'
relvdm_carboxamide.load_rel_vdms_pickle(sample, subset=no_charge)
relvdm_carboxamide.set_rel_vdm_bb_coords()
relvdm_carboxamide.set_rois_rot_trans(sample)
relvdm_carboxamide.set_rel_vdm_tags(sample)
print('moving vdMs')
relvdm_carboxamide.move_rel_vdms(sample)
print('removing clashing vdMs')
relvdm_carboxamide.remove_clash(sample)
relvdm_carboxamide.reshape_ifgs()
all_ifgs = functools.reduce(lambda a, b: np.vstack((a, b)),
                            [val for val in relvdm_carboxamide._ifgs.values()])
print('finding hotspots carboxamide')
nbrs = NearestNeighbors(metric='euclidean', radius=1.5, algorithm='kd_tree')
nbrs.fit(all_ifgs)
adj_mat = nbrs.radius_neighbors_graph(all_ifgs)
print('clustering...')
# mems, cents = combs.Cluster.cluster_adj_mat_gt(adj_mat, gt=20)
mems, cents = combs.analysis.cluster.greedy_cluster_pc(adj_mat, pc=0.8)
all_resnum_chid = functools.reduce(lambda a, b: np.vstack((a, b)),
                                   [np.array([tuple(key)]*len(val), dtype=object)
                                    for key, val in relvdm_carboxamide._vdm_tags.items() if len(val) > 0])
all_vdm_tags = functools.reduce(lambda a, b: np.vstack((a, b)),
                                [val for val in relvdm_carboxamide._vdm_tags.values() if len(val) > 0])
all_resn = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_carboxamide._resn.values() if len(val) > 0])
all_type = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_carboxamide._type.values() if len(val) > 0])
all_indices = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_carboxamide._indices.values() if len(val) > 0])
relvdm_carboxamide._all_vdm_tags = all_vdm_tags
relvdm_carboxamide._all_type = all_type
relvdm_carboxamide._all_resn = all_resn
relvdm_carboxamide._all_resnum_chid = all_resnum_chid
relvdm_carboxamide._all_ifgs = all_ifgs
print('carboxamide ifgs: ', len(all_ifgs))
relvdm_carboxamide._all_indices = all_indices

relvdm_carboxamide.rel_vdms_pickle = None
relvdm_carboxamide.rel_vdm_ifg_coords = None
relvdm_carboxamide.rel_vdm_bb_coords = None
relvdm_carboxamide.rel_vdm_sc_coords = None
relvdm_carboxamide.rel_vdm_tags = None
relvdm_carboxamide.rel_vdm_resnames = None
relvdm_carboxamide._ifgs = None
relvdm_carboxamide._resn = None
relvdm_carboxamide._type = None
relvdm_carboxamide._vdm_tags = None
relvdm_carboxamide._indices = None
relvdm_carboxamide.ifg_dict = {'GLN': 'CG CD OE1 NE2'}
relvdm_carboxamide._mems = mems
relvdm_carboxamide._cents = cents
relvdm_carboxamide.lig_ifg_correspondence = 'C10 C11 O1 N3'

# with open(outdir + 'relvdm_carboxamide.pickle', 'wb') as outfile:
#     pickle.dump(relvdm_carboxamide, outfile)


#
# inpath = '/Users/npolizzi/Projects/combs/results/apixaban/20171022/'
# with open(inpath + 'relvdm_backbone.pickle', 'rb') as infile:
#     relvdm_backbone = pickle.load(infile)
# with open(inpath + 'relvdm_carboxamide.pickle', 'rb') as infile:
#     relvdm_carboxamide = pickle.load(infile)
# with open(inpath + 'sample.pickle', 'rb') as infile:
#     sample = pickle.load(infile)

# carboxam_clusters_gt10 = [i for i, mem in enumerate(relvdm_carboxamide._mems) if len(mem) >= 10]
# ppc_clusters_gt = [i for i, mem in enumerate(relvdm_backbone._mems) if len(mem) >= 15]
# relvdm_carboxamide._cents = [relvdm_carboxamide._cents[i] for i in carboxam_clusters_gt10]
# relvdm_backbone._cents = [relvdm_backbone._cents[i] for i in ppc_clusters_gt]
# relvdm_backbone2._cents = relvdm_backbone._cents

# relvdm_backbone2 = copy.deepcopy(relvdm_backbone)
# relvdm_backbone2.lig_ifg_correspondence = 'C19 O2'
# relvdm_backbone2.name = 'backbone2'

relvdms = [relvdm_backbone, relvdm_carboxamide, relvdm_backbone2, relvdm_backbone3]
sample.relvdms = {}
sample.relvdms['carboxamide'] = relvdm_carboxamide
sample.relvdms['backbone'] = relvdm_backbone
sample.relvdms['backbone2'] = relvdm_backbone2
sample.relvdms['backbone3'] = relvdm_backbone3
print('finding cluster pairs')
sample.find_ligand_cst_cliques(relvdms, **{'maxlen': 2, 'atol1': 1.0, 'atol2': 1.0, 'pc': 1})

# sample.poi = None
# sample.rois = None

# with open(outdir + 'sample.pickle', 'wb') as outfile:
#     pickle.dump(sample, outfile)
# 
# sample.poi = pr.parsePDB(poi_path + '20171022/helix_template.pdb')

sample.set_protein_bb_coords()
print('making cluster densities')
sample.make_cluster_densities(relvdms)
print('fitting ligand to density')
print('number of cliques: ' + str(len(sample.cliques)))
sample.fit_lig_to_cluster_densities(relvdms, order=['carboxamide', 'backbone', 'backbone3', 'backbone2'], num_proc=8)


nbrs = {}
nbrs['carboxamide'] = NearestNeighbors(metric='euclidean', radius=1.5, algorithm='kd_tree').fit(relvdm_carboxamide._all_ifgs)
nbrs['backbone'] = NearestNeighbors(metric='euclidean', radius=1, algorithm='kd_tree').fit(relvdm_backbone._all_ifgs)
nbrs['backbone2'] = NearestNeighbors(metric='euclidean', radius=1, algorithm='kd_tree').fit(relvdm_backbone2._all_ifgs)
nbrs['backbone3'] = NearestNeighbors(metric='euclidean', radius=1, algorithm='kd_tree').fit(relvdm_backbone3._all_ifgs)

# for pose in sample.poses:
#     mems = {}
#     score_breakdown = {}
#     for name in pose.lig_ifg_coords.keys():
#         mems[name] = nbrs[name].radius_neighbors(pose.lig_ifg_coords[name], return_distance=False)[0]
#         score_breakdown[name] = len(mems[name])
#     pose.mems = mems
#     pose.score = np.prod(np.array(list(score_breakdown.values()))+np.ones(len(score_breakdown.values())))
#     pose.score_breakdown = score_breakdown
# scores = [p.score for p in sample.poses]
# scores_sorted = sorted(scores, reverse=True)
# inds_scores_sorted = sorted(range(len(scores)), reverse=True, key=lambda x: scores[x])
# sample.pose_scores = list(zip(scores_sorted, inds_scores_sorted))
#
# [pose.get_nonclashing(sample) for pose in sample.poses]
#
# for pose in sample.poses:
#     score_breakdown_prune = {}
#     for name in pose.mems_prune.keys():
#         score_breakdown_prune[name] = len(pose.mems_prune[name])
#     pose.score_prune = np.prod(np.array(list(score_breakdown_prune.values()))+np.ones(len(score_breakdown_prune.values())))
#     pose.score_breakdown_prune = score_breakdown_prune
# scores_prune = [p.score_prune for p in sample.poses]
# scores_prune_sorted = sorted(scores_prune, reverse=True)
# inds_scores_prune_sorted = sorted(range(len(scores_prune)), reverse=True, key=lambda x: scores_prune[x])
# sample.pose_scores_pruned = list(zip(scores_prune_sorted, inds_scores_prune_sorted))
#
# sample.poi = None
# sample.rois = None
# with open(outdir + 'sample_fit.pickle', 'wb') as outfile:
#     pickle.dump(sample, outfile)
#
# sample.poi = pr.parsePDB(poi_path + '20171022/helix_template.pdb')

for pose in sample.poses:
    mems = {}
    score_breakdown = {}
    for name in pose.lig_ifg_coords.keys():
        mems[name] = nbrs[name].radius_neighbors(pose.lig_ifg_coords[name], return_distance=False)[0]
        score_breakdown[name] = len(mems[name])
    pose.mems = mems
    pose.score = np.prod(np.array(list(score_breakdown.values()))+np.ones(len(score_breakdown.values())))
    pose.score_breakdown = score_breakdown
scores = [p.score for p in sample.poses]
scores_sorted = sorted(scores, reverse=True)
inds_scores_sorted = sorted(range(len(scores)), reverse=True, key=lambda x: scores[x])
sample.pose_scores = list(zip(scores_sorted, inds_scores_sorted))

[pose.get_nonclashing(sample) for pose in sample.poses]

for pose in sample.poses:
    score_breakdown_prune = {}
    for name in pose.mems_prune.keys():
        score_breakdown_prune[name] = len(pose.mems_prune[name])
    pose.score_prune = np.prod(np.array(list(score_breakdown_prune.values()))+np.ones(len(score_breakdown_prune.values())))
    pose.score_breakdown_prune = score_breakdown_prune
scores_prune = [p.score_prune for p in sample.poses]
scores_prune_sorted = sorted(scores_prune, reverse=True)
inds_scores_prune_sorted = sorted(range(len(scores_prune)), reverse=True, key=lambda x: scores_prune[x])
sample.pose_scores_pruned = list(zip(scores_prune_sorted, inds_scores_prune_sorted))

for pose in sample.poses:
    resnum_chid_all = collections.defaultdict(int)
    resnum_chid_sc = collections.defaultdict(int)
    resnum_chid_bb = collections.defaultdict(int)
    for name in pose.mems_prune.keys():
        for mem in pose.mems_prune[name]:
            rnch = tuple(sample.relvdms[name]._all_resnum_chid[mem])
            resnum_chid_all[rnch] += 1
            type_ =  sample.relvdms[name]._all_type[mem]
            if type_ == 'SC':
                resnum_chid_sc[rnch] += 1
            else:
                resnum_chid_bb[rnch] += 1
    pose.score_breakdown_resnum_chid_all = resnum_chid_all
    pose.score_breakdown_resnum_chid_sc = resnum_chid_sc
    pose.score_breakdown_resnum_chid_bb = resnum_chid_bb
    pose.score_resnum_chid_all = np.prod(np.array(list(resnum_chid_all.values())) + np.ones(len(resnum_chid_all.values())))
    pose.score_resnum_chid_sc = np.prod(np.array(list(resnum_chid_sc.values())) + np.ones(len(resnum_chid_sc.values())))
scores_resnum_chid_sc = [p.score_resnum_chid_sc for p in sample.poses]
scores_resnum_chid_sc_sorted = sorted(scores_resnum_chid_sc, reverse=True)
inds_scores_resnum_chid_sc_sorted = sorted(range(len(scores_resnum_chid_sc)), reverse=True, key=lambda x: scores_resnum_chid_sc[x])
sample.pose_scores_resnum_chid_sc = list(zip(scores_resnum_chid_sc_sorted, inds_scores_resnum_chid_sc_sorted))

scores_resnum_chid_all = [p.score_resnum_chid_all for p in sample.poses]
scores_resnum_chid_all_sorted = sorted(scores_resnum_chid_all, reverse=True)
inds_scores_resnum_chid_all_sorted = sorted(range(len(scores_resnum_chid_all)), reverse=True, key=lambda x: scores_resnum_chid_all[x])
sample.pose_scores_resnum_chid_all = list(zip(scores_resnum_chid_all_sorted, inds_scores_resnum_chid_all_sorted))

sample.poi = None
sample.rois = None
with open(outdir + 'sample_fit.pickle', 'wb') as outfile:
    pickle.dump(sample, outfile)

with open(outdir + 'designability_score.txt', 'w') as scorefile:
    for i in range(int(np.ceil(len(sample.pose_scores_resnum_chid_all) * 0.10))):
        scorefile.write(str(sample.pose_scores_resnum_chid_all[i][0]) + ', '
                        + str(sample.pose_scores_resnum_chid_all[i][1]) + '\n')

# sample.poi = pr.parsePDB(poi_path + '1.pdb', chain='A')
#
# pose_energies = []
# for i in range(int(np.ceil(len(sample.pose_scores_resnum_chid_all)*0.10))):
#     os.mkdir(str(i))
#     os.chdir(str(i))
#     pose_ind = sample.pose_scores_resnum_chid_all[i][1]
#     # sample.pose_clusters[pose_ind].make_sub_poses(sample, pdb_path + 'GNR.params')
#     sample.poses[pose_ind].make_sub_poses(sample, pdb_path + 'GNR.params')
#     # for val in sample.pose_clusters[pose_ind].opt_residues.values():
#     for val in sample.poses[pose_ind].opt_residues.values():
#         sample.relvdms[val[1]].print_cluster([val[0]], 'dee', '.')
#     # pr.writePDB('cluster/dee/lig.pdb', sample.pose_clusters[pose_ind].ligand)
#     pr.writePDB('cluster/dee/lig.pdb', sample.poses[pose_ind].ligand)
#
#     pose_dict = {}
#     # for key in sample.pose_clusters[pose_ind].score_breakdown_resnum_chid_bb.keys():
#     for key in sample.poses[pose_ind].score_breakdown_resnum_chid_bb.keys():
#         pose_dict[key] = ()
#     # for key, val in sample.pose_clusters[pose_ind].opt_residues.items():
#     for key, val in sample.poses[pose_ind].opt_residues.items():
#         pose_dict[key] = val
#     # sample.pose_clusters[pose_ind].pose_dict = pose_dict
#     # sample.pose_clusters[pose_ind].write_pose(sample,'.')
#     # scself, scpair, scpose = sample.pose_clusters[pose_ind].setup_score_fn(pdb_path + 'GNR.params')
#     sample.poses[pose_ind].pose_dict = pose_dict
#     sample.poses[pose_ind].write_pose(sample, '.')
#     scself, scpair, scpose = sample.poses[pose_ind].setup_score_fn(pdb_path + 'GNR.params')
#     pose = py.pose_from_pdb('pose.pdb')
#     pose_energies.append(scpose(pose))
#     os.chdir('..')
#
# ind_best_pose = np.argmin(pose_energies)
# print('best_pose=',ind_best_pose)
# # best_pose = sample.pose_clusters[ind_best_pose]
# best_pose = sample.poses[ind_best_pose]
# best_pose_py = py.pose_from_pdb(str(ind_best_pose) + '/pose.pdb')
# scpose.show(best_pose_py)