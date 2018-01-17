import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import prody as pr
import pickle
import os
import numpy as np
import functools
from sklearn.neighbors import NearestNeighbors
# import time


pdb_path = '/Users/npolizzi/Projects/combs/src/runs/glutamine_binding_protein/1wdn_noligandH.pdb'
sample = combs.Sample()
sample.poi = pr.parsePDB(pdb_path)
sample.bs_residues = list(zip([70, 75, 119, 118, 50, 67, 10, 115, 156, 13, 68, 157, 185],
                              ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A']))
sample.set_rois()
sample.set_rois_phipsi()
sample.set_poi_clash_coords()
sample.set_roi_bb_coords()

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
all_ifgs = functools.reduce(lambda a, b: np.vstack((a, b)),
                            [val for val in relvdm_amino._ifgs.values()])
print('finding hotspots')
nbrs = NearestNeighbors(metric='euclidean', radius=1.1, algorithm='kd_tree')
nbrs.fit(all_ifgs)
adj_mat = nbrs.radius_neighbors_graph(all_ifgs)
print('clustering...')
mems, cents = combs.Cluster.cluster_adj_mat_all_gt_beta(adj_mat, gt=5)
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
relvdm_amino._all_ifgs = all_ifgs

outdir = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20171016/'
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
with open(outdir + 'relvdm_amino.pickle', 'wb') as outfile:
    pickle.dump(relvdm_amino, outfile)

