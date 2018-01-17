import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import prody as pr
import pickle
import numpy as np
import functools
from sklearn.neighbors import NearestNeighbors
import os


pdb_path = '/Users/npolizzi/Projects/combs/src/runs/glutamine_binding_protein/'
sample = combs.Sample()
sample.poi = pr.parsePDB(pdb_path + '1wdn_noligandH.pdb')
sample.bs_residues = list(zip([70, 75, 119, 118, 50, 67, 10, 115, 156, 13, 68, 157, 185],
                              ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A']))
sample.set_rois()
sample.set_rois_phipsi()
sample.set_poi_clash_coords()
sample.set_roi_bb_coords()
sample.ligand_conformers = [pr.parsePDB(pdb_path + 'confs/model0.pdb')]

relvdm_amino = combs.Rel_Vandermer.RelVandermer('amino')
relvdm_amino.rel_vdm_path = '/Users/npolizzi/Projects/combs/results/amino/rel_vdms_hbond/20171015/'
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
mems, cents = combs.Cluster.cluster_adj_mat_all_gt(adj_mat, gt=4)


all_resnum_chid = functools.reduce(lambda a, b: np.vstack((a, b)),
                                   [np.array([tuple(key)]*len(val), dtype=object)
                                    for key, val in relvdm_amino._vdm_tags.items()])
all_vdm_tags = functools.reduce(lambda a, b: np.vstack((a, b)),
                                [val for val in relvdm_amino._vdm_tags.values()])
all_resn = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_amino._resn.values()])
all_type = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_amino._type.values()])
relvdm_amino._all_vdm_tags = all_vdm_tags
relvdm_amino._all_type = all_type
relvdm_amino._all_resn = all_resn
relvdm_amino._all_resnum_chid = all_resnum_chid
relvdm_amino._all_ifgs = all_ifgs_amino

outdir = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20171018/'
try:
    os.makedirs(outdir)
except:
    pass

relvdm_amino.rel_vdms_pickle = None
relvdm_amino.rel_vdm_ifg_coords = None
relvdm_amino.rel_vdm_bb_coords = None
relvdm_amino.rel_vdm_sc_coords = None
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


relvdm_carboxamide = combs.Rel_Vandermer.RelVandermer('carboxamide')
relvdm_carboxamide.rel_vdm_path = '/Users/npolizzi/Projects/combs/results/carboxamide/rel_vdms_hbond/20171015/'
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
mems, cents = combs.Cluster.cluster_adj_mat_all_gt(adj_mat, gt=4)
all_resnum_chid = functools.reduce(lambda a, b: np.vstack((a, b)),
                                   [np.array([tuple(key)]*len(val), dtype=object)
                                    for key, val in relvdm_carboxamide._vdm_tags.items()])
all_vdm_tags = functools.reduce(lambda a, b: np.vstack((a, b)),
                                [val for val in relvdm_carboxamide._vdm_tags.values()])
all_resn = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_carboxamide._resn.values()])
all_type = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_carboxamide._type.values()])
relvdm_carboxamide._all_vdm_tags = all_vdm_tags
relvdm_carboxamide._all_type = all_type
relvdm_carboxamide._all_resn = all_resn
relvdm_carboxamide._all_resnum_chid = all_resnum_chid
relvdm_carboxamide._all_ifgs = all_ifgs

outdir = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20171018/'
try:
    os.makedirs(outdir)
except:
    pass

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
relvdm_carboxamide.lig_ifg_correspondence = 'CG CD OE1 NE2'

with open(outdir + 'relvdm_carboxamide.pickle', 'wb') as outfile:
    pickle.dump(relvdm_carboxamide, outfile)


relvdm_carboxylate = combs.Rel_Vandermer.RelVandermer('carboxylate')
relvdm_carboxylate.rel_vdm_path = '/Users/npolizzi/Projects/combs/results/carboxylate/rel_vdms_hbond/20171015/'
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
mems, cents = combs.Cluster.cluster_adj_mat_all_gt(adj_mat, gt=4)
all_resnum_chid = functools.reduce(lambda a, b: np.vstack((a, b)),
                                   [np.array([tuple(key)]*len(val), dtype=object)
                                    for key, val in relvdm_carboxylate._vdm_tags.items()])
all_vdm_tags = functools.reduce(lambda a, b: np.vstack((a, b)),
                                [val for val in relvdm_carboxylate._vdm_tags.values()])
all_resn = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_carboxylate._resn.values()])
all_type = functools.reduce(lambda a, b: np.hstack((a, b)),
                            [val for val in relvdm_carboxylate._type.values()])
relvdm_carboxylate._all_vdm_tags = all_vdm_tags
relvdm_carboxylate._all_type = all_type
relvdm_carboxylate._all_resn = all_resn
relvdm_carboxylate._all_resnum_chid = all_resnum_chid
relvdm_carboxylate._all_ifgs = all_ifgs


outdir = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20171018/'
try:
    os.makedirs(outdir)
except:
    pass

relvdm_carboxylate.rel_vdms_pickle = None
relvdm_carboxylate.rel_vdm_ifg_coords = None
relvdm_carboxylate.rel_vdm_bb_coords = None
relvdm_carboxylate.rel_vdm_sc_coords = None
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

print('finding cluster pairs')
sample.find_ligand_cst_cliques(relvdms, **{'atol': 0.6})

sample.poi = None
sample.rois = None
outpath = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20171018/'
with open(outpath + 'sample.pickle', 'wb') as outfile:
    pickle.dump(sample, outfile)

sample.poi = pr.parsePDB(pdb_path + '1wdn_noligandH.pdb')
sample.set_protein_bb_coords()
print('making cluster densities')
sample.make_cluster_densities(relvdms)
print('fitting ligand to density')
sample.fit_lig_to_cluster_density(relvdms)

sample.poi = None
outpath = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20171018/'
with open(outpath + 'sample_fit.pickle', 'wb') as outfile:
    pickle.dump(sample, outfile)