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


pdb_path = '/Users/npolizzi/Projects/combs/src/runs/apixaban/'
sample = combs.Sample()
poi_path = '/Users/npolizzi/Projects/combs/src/runs/apixaban/'

# site_path = '/Users/npolizzi/Downloads/scPDB_2013/sc_gt80_le300_MWgt400_sitele35res_sites_pdbs/'

sample.poi = pr.parsePDB(poi_path + '5tgy_3_truncH.pdb')
# site = pr.parsePDB(site_path + '5tgy_3_trunc.pdb')

# site_CA = site.select('name CA')

# sels = []
# for ca in site_CA:
#     sels.append(sample.poi.select('exwithin .1 of coords', coords=ca.getCoords()))
resis = [2,3,6,7,9,10,13,14,16,17,18,20,36,39,40,42,43,46,47,49,50,53,54,55,56,57,58,61,62,64,65,68,69,71,72,73,76]

sample.bs_residues = list(zip(resis, ['A']*len(resis)))
print(sample.bs_residues)
sample.set_rois()
sample.set_rois_phipsi()
sample.set_poi_clash_coords()
sample.set_roi_bb_coords()
sample.ligand_conformers = [pr.parsePDB(pdb_path + 'apix.pdb')]

relvdm_preprolinecarbonyl = combs.Rel_Vandermer.RelVandermer('preprolinecarbonyl')
relvdm_preprolinecarbonyl.rel_vdm_path = '/Users/npolizzi/Projects/combs/results/preproline_carbonyl/rel_vdms_hbond_carbonyl/20171022/'

no_charge = set(combs.constants.one_letter_code.keys()) - {'LYS', 'ARG', 'GLU', 'ASP', 'MSE', 'ANY',
                                                           'PRO', 'MET', 'CYS'}
relvdm_preprolinecarbonyl.load_rel_vdms_pickle(sample, subset=no_charge)
relvdm_preprolinecarbonyl.set_rel_vdm_bb_coords()
relvdm_preprolinecarbonyl.set_rois_rot_trans(sample)
relvdm_preprolinecarbonyl.set_rel_vdm_tags(sample)
print('moving vdMs')
relvdm_preprolinecarbonyl.move_rel_vdms(sample)
print('removing clashing vdMs')
relvdm_preprolinecarbonyl.remove_clash(sample)
relvdm_preprolinecarbonyl.reshape_ifgs()
all_ifgs = functools.reduce(lambda a, b: np.vstack((a, b)),
                            [val for val in relvdm_preprolinecarbonyl._ifgs.values()])
print('finding hotspots preproline carbonyl')
nbrs = NearestNeighbors(metric='euclidean', radius=1.0, algorithm='kd_tree')
nbrs.fit(all_ifgs)
adj_mat = nbrs.radius_neighbors_graph(all_ifgs)
print('clustering...')

mems, cents = combs.Cluster.greedy_cluster_pc(adj_mat, pc=0.7)

all_resnum_chid = functools.reduce(lambda a, b: np.vstack((a, b)),
                                   [np.array([tuple(key)]*len(val), dtype=object)
                                    for key, val in relvdm_preprolinecarbonyl._vdm_tags.items() if len(val) > 0])
all_vdm_tags = functools.reduce(lambda a, b: np.vstack((a, b)),
                                [val for val in relvdm_preprolinecarbonyl._vdm_tags.values() if len(val) > 0])
all_resn = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_preprolinecarbonyl._resn.values() if len(val) > 0])
all_type = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_preprolinecarbonyl._type.values() if len(val) > 0])
all_indices = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_preprolinecarbonyl._indices.values() if len(val) > 0])
relvdm_preprolinecarbonyl._all_vdm_tags = all_vdm_tags
relvdm_preprolinecarbonyl._all_type = all_type
relvdm_preprolinecarbonyl._all_resn = all_resn
relvdm_preprolinecarbonyl._all_resnum_chid = all_resnum_chid
relvdm_preprolinecarbonyl._all_ifgs = all_ifgs
print('ppc ifgs: ', len(all_ifgs))
relvdm_preprolinecarbonyl._all_indices = all_indices



relvdm_preprolinecarbonyl.rel_vdms_pickle = None
relvdm_preprolinecarbonyl.rel_vdm_ifg_coords = None
relvdm_preprolinecarbonyl.rel_vdm_bb_coords = None
relvdm_preprolinecarbonyl.rel_vdm_sc_coords = None
relvdm_preprolinecarbonyl.rel_vdm_tags = None
relvdm_preprolinecarbonyl.rel_vdm_resnames = None
relvdm_preprolinecarbonyl._ifgs = None
relvdm_preprolinecarbonyl._resn = None
relvdm_preprolinecarbonyl._type = None
relvdm_preprolinecarbonyl._vdm_tags = None
relvdm_preprolinecarbonyl._indices = None
relvdm_preprolinecarbonyl.ifg_dict = {'ANY': 'C O'}
relvdm_preprolinecarbonyl._mems = mems
relvdm_preprolinecarbonyl._cents = cents
relvdm_preprolinecarbonyl.lig_ifg_correspondence = 'C8 O3'

outdir = '/Users/npolizzi/Projects/combs/results/apixaban/20171025/5tgy_3_trunc/ppc/'
try:
    os.makedirs(outdir)
except:
    pass

# with open(outdir + 'relvdm_preprolinecarbonyl.pickle', 'wb') as outfile:
#     pickle.dump(relvdm_preprolinecarbonyl, outfile)


relvdm_preprolinecarbonyl2 = copy.deepcopy(relvdm_preprolinecarbonyl)
relvdm_preprolinecarbonyl2.lig_ifg_correspondence = 'C19 O2'
relvdm_preprolinecarbonyl2.name = 'preprolinecarbonyl2'

relvdm_preprolinecarbonyl3 = copy.deepcopy(relvdm_preprolinecarbonyl)
relvdm_preprolinecarbonyl3.lig_ifg_correspondence = 'N1 N6'
relvdm_preprolinecarbonyl3.name = 'preprolinecarbonyl3'

relvdm_carboxamide = combs.Rel_Vandermer.RelVandermer('carboxamide')
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
mems, cents = combs.Cluster.greedy_cluster_pc(adj_mat, pc=0.7)
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
# with open(inpath + 'relvdm_preprolinecarbonyl.pickle', 'rb') as infile:
#     relvdm_preprolinecarbonyl = pickle.load(infile)
# with open(inpath + 'relvdm_carboxamide.pickle', 'rb') as infile:
#     relvdm_carboxamide = pickle.load(infile)
# with open(inpath + 'sample.pickle', 'rb') as infile:
#     sample = pickle.load(infile)

# carboxam_clusters_gt10 = [i for i, mem in enumerate(relvdm_carboxamide._mems) if len(mem) >= 10]
# ppc_clusters_gt = [i for i, mem in enumerate(relvdm_preprolinecarbonyl._mems) if len(mem) >= 15]
# relvdm_carboxamide._cents = [relvdm_carboxamide._cents[i] for i in carboxam_clusters_gt10]
# relvdm_preprolinecarbonyl._cents = [relvdm_preprolinecarbonyl._cents[i] for i in ppc_clusters_gt]
# relvdm_preprolinecarbonyl2._cents = relvdm_preprolinecarbonyl._cents

# relvdm_preprolinecarbonyl2 = copy.deepcopy(relvdm_preprolinecarbonyl)
# relvdm_preprolinecarbonyl2.lig_ifg_correspondence = 'C19 O2'
# relvdm_preprolinecarbonyl2.name = 'preprolinecarbonyl2'

relvdms = [relvdm_preprolinecarbonyl, relvdm_carboxamide, relvdm_preprolinecarbonyl2, relvdm_preprolinecarbonyl3]
sample.relvdms = {}
sample.relvdms['carboxamide'] = relvdm_carboxamide
sample.relvdms['preprolinecarbonyl'] = relvdm_preprolinecarbonyl
sample.relvdms['preprolinecarbonyl2'] = relvdm_preprolinecarbonyl2
sample.relvdms['preprolinecarbonyl3'] = relvdm_preprolinecarbonyl3
print('finding cluster pairs')
sample.find_ligand_cst_cliques(relvdms, **{'maxlen': 3, 'atol1': 1.4, 'atol2': 1.4, 'pc': 1})

# sample.poi = None
# sample.rois = None
#
# with open(outdir + 'sample.pickle', 'wb') as outfile:
#     pickle.dump(sample, outfile)

# sample.poi = pr.parsePDB(poi_path + '5tgy_3_truncH.pdb')

sample.set_protein_bb_coords()
print('making cluster densities')
sample.make_cluster_densities(relvdms)
print('fitting ligand to density')
print('number of cliques: ' + str(len(sample.cliques)))
sample.fit_lig_to_cluster_densities(relvdms, order=['carboxamide', 'preprolinecarbonyl', 'preprolinecarbonyl3', 'preprolinecarbonyl2'], num_proc=8)


nbrs = {}
nbrs['carboxamide'] = NearestNeighbors(metric='euclidean', radius=1.5, algorithm='kd_tree').fit(relvdm_carboxamide._all_ifgs)
nbrs['preprolinecarbonyl'] = NearestNeighbors(metric='euclidean', radius=1, algorithm='kd_tree').fit(relvdm_preprolinecarbonyl._all_ifgs)
nbrs['preprolinecarbonyl2'] = NearestNeighbors(metric='euclidean', radius=1, algorithm='kd_tree').fit(relvdm_preprolinecarbonyl2._all_ifgs)
nbrs['preprolinecarbonyl3'] = NearestNeighbors(metric='euclidean', radius=1, algorithm='kd_tree').fit(relvdm_preprolinecarbonyl3._all_ifgs)

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

sample.poi = None
sample.rois = None
with open(outdir + 'sample_fit.pickle', 'wb') as outfile:
    pickle.dump(sample, outfile)

sample.poi = pr.parsePDB(poi_path + '5tgy_3_truncH.pdb')