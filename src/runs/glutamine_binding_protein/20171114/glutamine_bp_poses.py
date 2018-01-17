import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import prody as pr
import pickle
import numpy as np
import functools
from sklearn.neighbors import NearestNeighbors
import os
import collections
import pyrosetta as py


pdb_path = '/Users/npolizzi/Projects/combs/src/runs/glutamine_binding_protein/'
sample = combs.Sample()
sample.poi = pr.parsePDB(pdb_path + '1wdn_noligandH.pdb')
sample.bs_residues = list(zip([70, 75, 119, 118, 50, 67, 10, 115, 156, 13, 68, 157, 185],
                              ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A']))
sample.set_rois()
sample.set_rois_phipsi()
sample.set_poi_clash_coords()
sample.set_roi_bb_coords()
sample.ligand_conformers = [pr.parsePDB(pdb_path + 'GNR_0001.pdb')]

relvdm_amino = combs.rel_vandermer.RelVandermer('amino')
relvdm_amino.rel_vdm_path = '/Users/npolizzi/Projects/combs/results/amino/rel_vdms_hbond/20171025/'
relvdm_amino.load_rel_vdms_pickle(sample)
relvdm_amino.set_rel_vdm_bb_coords()
relvdm_amino.set_rois_rot_trans(sample)
relvdm_amino.set_rel_vdm_tags(sample)
print('moving vdMs')
relvdm_amino.move_rel_vdms(sample)
print('removing clashing vdMs')
relvdm_amino.remove_clash(sample)
relvdm_amino.reshape_ifgs()
all_ifgs_amino = functools.reduce(lambda a, b: np.vstack((a, b)),
                            [val for val in relvdm_amino._ifgs.values()])
print('finding hotspots amino')
nbrs = NearestNeighbors(metric='euclidean', radius=1.1, algorithm='kd_tree')
nbrs.fit(all_ifgs_amino)
adj_mat = nbrs.radius_neighbors_graph(all_ifgs_amino)
print('clustering...')
mems, cents = combs.analysis.cluster.greedy_cluster_pc(adj_mat, pc=0.7)


all_resnum_chid = functools.reduce(lambda a, b: np.vstack((a, b)),
                                   [np.array([tuple(key)]*len(val), dtype=object)
                                    for key, val in relvdm_amino._vdm_tags.items()])
all_vdm_tags = functools.reduce(lambda a, b: np.vstack((a, b)),
                                [val for val in relvdm_amino._vdm_tags.values()])
all_resn = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_amino._resn.values()])
all_type = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_amino._type.values()])
all_indices = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_amino._indices.values() if len(val) > 0])
relvdm_amino._all_vdm_tags = all_vdm_tags
relvdm_amino._all_type = all_type
relvdm_amino._all_resn = all_resn
relvdm_amino._all_resnum_chid = all_resnum_chid
relvdm_amino._all_ifgs = all_ifgs_amino
relvdm_amino._all_indices = all_indices
print('amino ifgs: ', len(all_ifgs_amino))

outdir = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20171114/'
try:
    os.makedirs(outdir)
except:
    pass

relvdm_amino.rel_vdms_pickle = None
relvdm_amino.rel_vdm_ifg_coords = None
# relvdm_amino.rel_vdm_bb_coords = None
# relvdm_amino.rel_vdm_sc_coords = None
relvdm_amino.rel_vdm_tags = None
relvdm_amino.rel_vdm_resnames = None
relvdm_amino._ifgs = None
relvdm_amino._resn = None
relvdm_amino._type = None
relvdm_amino._vdm_tags = None
relvdm_amino._indices = None
relvdm_amino.ifg_dict = {'LYS': 'CE NZ'}
relvdm_amino._mems = mems
relvdm_amino._cents = cents
relvdm_amino.lig_ifg_correspondence = 'CA N'

with open(outdir + 'relvdm_amino.pickle', 'wb') as outfile:
    pickle.dump(relvdm_amino, outfile)


relvdm_carboxamide = combs.rel_vandermer.RelVandermer('carboxamide')
relvdm_carboxamide.rel_vdm_path = '/Users/npolizzi/Projects/combs/results/carboxamide/rel_vdms_hbond/20171023/'
relvdm_carboxamide.load_rel_vdms_pickle(sample)
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
mems, cents = combs.analysis.cluster.greedy_cluster_pc(adj_mat, pc=0.7)
all_resnum_chid = functools.reduce(lambda a, b: np.vstack((a, b)),
                                   [np.array([tuple(key)]*len(val), dtype=object)
                                    for key, val in relvdm_carboxamide._vdm_tags.items()])
all_vdm_tags = functools.reduce(lambda a, b: np.vstack((a, b)),
                                [val for val in relvdm_carboxamide._vdm_tags.values()])
all_resn = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_carboxamide._resn.values()])
all_type = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_carboxamide._type.values()])
all_indices = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_carboxamide._indices.values() if len(val) > 0])
relvdm_carboxamide._all_vdm_tags = all_vdm_tags
relvdm_carboxamide._all_type = all_type
relvdm_carboxamide._all_resn = all_resn
relvdm_carboxamide._all_resnum_chid = all_resnum_chid
relvdm_carboxamide._all_ifgs = all_ifgs
relvdm_carboxamide._all_indices = all_indices
print('carboxamide ifgs: ', len(all_ifgs))

outdir = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20171114/'
try:
    os.makedirs(outdir)
except:
    pass

relvdm_carboxamide.rel_vdms_pickle = None
relvdm_carboxamide.rel_vdm_ifg_coords = None
# relvdm_carboxamide.rel_vdm_bb_coords = None
# relvdm_carboxamide.rel_vdm_sc_coords = None
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
relvdm_carboxamide.lig_ifg_correspondence = 'CG CD OE1 NE2'

with open(outdir + 'relvdm_carboxamide.pickle', 'wb') as outfile:
    pickle.dump(relvdm_carboxamide, outfile)


relvdm_carboxylate = combs.rel_vandermer.RelVandermer('carboxylate')
relvdm_carboxylate.rel_vdm_path = '/Users/npolizzi/Projects/combs/results/carboxylate/rel_vdms_hbond/20171025/'
relvdm_carboxylate.load_rel_vdms_pickle(sample)
relvdm_carboxylate.set_rel_vdm_bb_coords()
relvdm_carboxylate.set_rois_rot_trans(sample)
relvdm_carboxylate.set_rel_vdm_tags(sample)
print('moving vdMs')
relvdm_carboxylate.move_rel_vdms(sample)
print('removing clashing vdMs')
relvdm_carboxylate.remove_clash(sample)
relvdm_carboxylate.reshape_ifgs()
all_ifgs = functools.reduce(lambda a, b: np.vstack((a, b)),
                            [val for val in relvdm_carboxylate._ifgs.values()])
print('finding hotspots carboxylate')
nbrs = NearestNeighbors(metric='euclidean', radius=1.5, algorithm='kd_tree')
nbrs.fit(all_ifgs)
adj_mat = nbrs.radius_neighbors_graph(all_ifgs)
print('clustering...')
mems, cents = combs.analysis.cluster.greedy_cluster_pc(adj_mat, pc=0.7)
all_resnum_chid = functools.reduce(lambda a, b: np.vstack((a, b)),
                                   [np.array([tuple(key)]*len(val), dtype=object)
                                    for key, val in relvdm_carboxylate._vdm_tags.items()])
all_vdm_tags = functools.reduce(lambda a, b: np.vstack((a, b)),
                                [val for val in relvdm_carboxylate._vdm_tags.values()])
all_resn = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_carboxylate._resn.values()])
all_type = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_carboxylate._type.values()])
all_indices = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_carboxylate._indices.values() if len(val) > 0])
relvdm_carboxylate._all_vdm_tags = all_vdm_tags
relvdm_carboxylate._all_type = all_type
relvdm_carboxylate._all_resn = all_resn
relvdm_carboxylate._all_resnum_chid = all_resnum_chid
relvdm_carboxylate._all_ifgs = all_ifgs
relvdm_carboxylate._all_indices = all_indices
print('carboxylate ifgs: ', len(all_ifgs))


outdir = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20171114/'
try:
    os.makedirs(outdir)
except:
    pass

relvdm_carboxylate.rel_vdms_pickle = None
relvdm_carboxylate.rel_vdm_ifg_coords = None
# relvdm_carboxylate.rel_vdm_bb_coords = None
# relvdm_carboxylate.rel_vdm_sc_coords = None
relvdm_carboxylate.rel_vdm_tags = None
relvdm_carboxylate.rel_vdm_resnames = None
relvdm_carboxylate._ifgs = None
relvdm_carboxylate._resn = None
relvdm_carboxylate._type = None
relvdm_carboxylate._vdm_tags = None
relvdm_carboxylate._indices = None
relvdm_carboxylate.ifg_dict = {'GLU': 'CG CD OE1 OE2'}
relvdm_carboxylate._mems = mems
relvdm_carboxylate._cents = cents
relvdm_carboxylate.lig_ifg_correspondence = 'CA C O OXT'

with open(outdir + 'relvdm_carboxylate.pickle', 'wb') as outfile:
    pickle.dump(relvdm_carboxylate, outfile)


relvdms = [relvdm_amino, relvdm_carboxamide, relvdm_carboxylate]
sample.relvdms = {}
sample.relvdms['carboxamide'] = relvdm_carboxamide
sample.relvdms['carboxylate'] = relvdm_carboxylate
sample.relvdms['amino'] = relvdm_amino
print('finding cluster pairs')
sample.find_ligand_cst_cliques(relvdms, **{'maxlen': 3, 'atol1': 1.5, 'atol2': 1.5, 'pc': 1})

sample.poi = None
sample.rois = None
outpath = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20171114/'
with open(outpath + 'sample.pickle', 'wb') as outfile:
    pickle.dump(sample, outfile)

sample.poi = pr.parsePDB(pdb_path + '1wdn_noligandH.pdb')
sample.set_protein_bb_coords()
print('making cluster densities')
sample.make_cluster_densities(relvdms)
print('fitting ligand to density')
print('number of cliques: ' + str(len(sample.cliques)))
sample.fit_lig_to_cluster_densities(relvdms, order=['carboxylate', 'amino', 'carboxamide'], num_proc=8)

sample.poi = None
outpath = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20171114/'
with open(outpath + 'sample_fit.pickle', 'wb') as outfile:
    pickle.dump(sample, outfile)
    
sample.poi = pr.parsePDB(pdb_path + '1wdn_noligandH.pdb')

nbrs = {}
nbrs['carboxamide'] = NearestNeighbors(metric='euclidean', radius=1.5, algorithm='kd_tree').fit(relvdm_carboxamide._all_ifgs)
nbrs['carboxylate'] = NearestNeighbors(metric='euclidean', radius=1.5, algorithm='kd_tree').fit(relvdm_carboxylate._all_ifgs)
nbrs['amino'] = NearestNeighbors(metric='euclidean', radius=1.1, algorithm='kd_tree').fit(relvdm_amino._all_ifgs)


# lig_coords = [pose.ligand.getCoords().reshape(-1) for pose in sample.poses]
# num_heavy_lig = len(sample.ligand_conformers[0].select('not element H'))
# nbrs_lig = NearestNeighbors(metric='euclidean', radius=0.4*np.sqrt(num_heavy_lig))
# nbrs_lig.fit(lig_coords)
# adj_mat_lig = nbrs_lig.radius_neighbors_graph(lig_coords)
# pose_clusters, pose_cents = combs.analysis.cluster.greedy_cluster(adj_mat_lig)
# #cluster poses based on mems.
# sample.pose_clusters = []
# for p_cluster, p_cent in zip(pose_clusters, pose_cents):
#     mems = collections.defaultdict(set)
#     score_breakdown = {}
#     for p in p_cluster:
#         pose = sample.poses[p]
#         for name in pose.lig_ifg_coords.keys():
#             m = nbrs[name].radius_neighbors(pose.lig_ifg_coords[name], return_distance=False)[0]
#             [mems[name].add(mem) for mem in m]
#     for name in mems.keys():
#         score_breakdown[name] = len(mems[name])
#         mems[name] = np.array(list(mems[name]))
#     super_pose = combs.design.pose.Pose()
#     super_pose.mems = mems
#     super_pose.score_breakdown = score_breakdown
#     super_pose.ligand = sample.poses[p_cent].ligand
#     sample.pose_clusters.append(super_pose)
#
# [super_pose.get_nonclashing(sample) for super_pose in sample.pose_clusters]

lig_coords = [pose.ligand.getCoords().reshape(-1) for pose in sample.poses]
num_heavy_lig = len(sample.ligand_conformers[0].select('not element H'))
nbrs_lig = NearestNeighbors(metric='euclidean', radius=0.3*np.sqrt(num_heavy_lig))
nbrs_lig.fit(lig_coords)
adj_mat_lig = nbrs_lig.radius_neighbors_graph(lig_coords)
pose_clusters, pose_cents = combs.analysis.cluster.greedy_cluster(adj_mat_lig)
#cluster poses based on mems.
sample.pose_clusters = []
for p_cluster, p_cent in zip(pose_clusters, pose_cents):
    mems = collections.defaultdict(set)
    score_breakdown = {}
    ligands = []
    for p in p_cluster:
        pose = sample.poses[p]
        ligands.append(pose.ligand)
        for name in pose.lig_ifg_coords.keys():
            m = nbrs[name].radius_neighbors(pose.lig_ifg_coords[name], return_distance=False)[0]
            [mems[name].add(mem) for mem in m]
    for name in mems.keys():
        score_breakdown[name] = len(mems[name])
        mems[name] = np.array(list(mems[name]))
    super_pose = combs.design.pose.Pose()
    super_pose.mems = mems
    super_pose.score_breakdown = score_breakdown
    scores = []
    for ligand in ligands:
        super_pose.ligand = ligand
        super_pose.get_nonclashing(sample)
        score_breakdown_prune = {}
        for name in super_pose.mems_prune.keys():
            score_breakdown_prune[name] = len(super_pose.mems_prune[name])
        score_prune = np.prod(np.array(list(score_breakdown_prune.values())) + np.ones(len(score_breakdown_prune.values())))
        scores.append(score_prune)
    best_ligand = ligands[np.argmax(scores)]
    super_pose.ligand = best_ligand
    super_pose.get_nonclashing(sample)
    sample.pose_clusters.append(super_pose)

# [super_pose.get_nonclashing(sample) for super_pose in sample.pose_clusters]



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

for pose in sample.pose_clusters:
    score_breakdown_prune = {}
    for name in pose.mems_prune.keys():
        score_breakdown_prune[name] = len(pose.mems_prune[name])
    pose.score_prune = np.prod(np.array(list(score_breakdown_prune.values()))+np.ones(len(score_breakdown_prune.values())))
    pose.score_breakdown_prune = score_breakdown_prune
scores_prune = [p.score_prune for p in sample.pose_clusters]
scores_prune_sorted = sorted(scores_prune, reverse=True)
inds_scores_prune_sorted = sorted(range(len(scores_prune)), reverse=True, key=lambda x: scores_prune[x])
sample.pose_scores_pruned = list(zip(scores_prune_sorted, inds_scores_prune_sorted))

for pose in sample.pose_clusters:
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
scores_resnum_chid_sc = [p.score_resnum_chid_sc for p in sample.pose_clusters]
scores_resnum_chid_sc_sorted = sorted(scores_resnum_chid_sc, reverse=True)
inds_scores_resnum_chid_sc_sorted = sorted(range(len(scores_resnum_chid_sc)), reverse=True, key=lambda x: scores_resnum_chid_sc[x])
sample.pose_scores_resnum_chid_sc = list(zip(scores_resnum_chid_sc_sorted, inds_scores_resnum_chid_sc_sorted))

scores_resnum_chid_all = [p.score_resnum_chid_all for p in sample.pose_clusters]
scores_resnum_chid_all_sorted = sorted(scores_resnum_chid_all, reverse=True)
inds_scores_resnum_chid_all_sorted = sorted(range(len(scores_resnum_chid_all)), reverse=True, key=lambda x: scores_resnum_chid_all[x])
sample.pose_scores_resnum_chid_all = list(zip(scores_resnum_chid_all_sorted, inds_scores_resnum_chid_all_sorted))

# lig_rmsds = [pr.calcRMSD(sample.ligand_conformers[0], pose.ligand) for pose in sample.poses]

# lig_coords = [pose.ligand.getCoords().reshape(-1) for pose in sample.poses]
# num_heavy_lig = len(sample.ligand_conformers[0].select('not element H'))
# nbrs = NearestNeighbors(metric='euclidean', radius=0.4*np.sqrt(num_heavy_lig))
# nbrs.fit(lig_coords)
# adj_mat = nbrs.radius_neighbors_graph(lig_coords)
# pose_clusters, pose_cents = combs.analysis.cluster.greedy_cluster(adj_mat)
# #cluster poses based on mems.
# sample.pose_clusters = []
# for p_cluster, p_cent in zip(pose_clusters, pose_cents):
#     clust = collections.defaultdict(set)
#     for p in p_cluster:
#         for name, mems in sample.poses[p].mems_prune.items():
#             [clust[name].add(mem) for mem in mems]
#     super_pose = combs.design.pose.Pose()
#     super_pose.mems_prune
    #sample.pose_clusters.append([clust, p_cent])
    # order of operations: 1. make super poses.
    # 2. get rid of members clashing with centroid ligand.
    # 3. calculate scores.



sample.poi = None
with open(outdir + 'sample_fit.pickle', 'wb') as outfile:
    pickle.dump(sample, outfile)

sample.poi = pr.parsePDB(pdb_path + '1wdn_noligandH.pdb')

# pose_ind = sample.pose_scores_pruned[0][1]
pose_energies = []
for i in range(int(np.ceil(len(sample.pose_scores_resnum_chid_all)*0.25))):
    os.mkdir(str(i))
    os.chdir(str(i))
    pose_ind = sample.pose_scores_resnum_chid_all[i][1]
    sample.pose_clusters[pose_ind].make_sub_poses(sample, pdb_path + 'GNR.params')
    for val in sample.pose_clusters[pose_ind].opt_residues.values():
        sample.relvdms[val[1]].print_cluster([val[0]], 'dee', '.')
    pr.writePDB('cluster/dee/lig.pdb', sample.pose_clusters[pose_ind].ligand)

    pose_dict = {}
    for key in sample.pose_clusters[pose_ind].score_breakdown_resnum_chid_bb.keys():
        pose_dict[key] = ()
    for key, val in sample.pose_clusters[pose_ind].opt_residues.items():
        pose_dict[key] = val
    sample.pose_clusters[pose_ind].pose_dict = pose_dict
    sample.pose_clusters[pose_ind].write_pose(sample,'.')
    scself, scpair, scpose = sample.pose_clusters[pose_ind].setup_score_fn(pdb_path + 'GNR.params')
    pose = py.pose_from_pdb('pose.pdb')
    pose_energies.append(scpose(pose))
    os.chdir('..')

ind_best_pose = np.argmin(pose_energies)
print('best_pose=',ind_best_pose)
best_pose = sample.pose_clusters[ind_best_pose]
best_pose_py = py.pose_from_pdb(str(ind_best_pose) + '/pose.pdb')
scpose.show(best_pose_py)


# pose_ind = sample.pose_scores_resnum_chid_sc[0][1]
# sample.pose_clusters[pose_ind].make_sub_poses(sample, pdb_path + 'GNR.params')
# for val in sample.pose_clusters[pose_ind].opt_residues.values():
#     sample.relvdms[val[1]].print_cluster([val[0]], 'dee', '.')
# pr.writePDB('cluster/dee/lig.pdb', sample.pose_clusters[pose_ind].ligand)
#
# pose_dict = {}
# for key in sample.pose_clusters[pose_ind].score_breakdown_resnum_chid_bb.keys():
#     pose_dict[key] = ()
# for key, val in sample.pose_clusters[pose_ind].opt_residues.items():
#     pose_dict[key] = val
# sample.pose_clusters[pose_ind].pose_dict = pose_dict
# sample.pose_clusters[pose_ind].write_pose(sample, '.')
# scself, scpair, scpose = sample.pose_clusters[pose_ind].setup_score_fn(pdb_path + 'GNR.params')
# pose = py.pose_from_pdb('pose.pdb')
# scpose.show(pose)