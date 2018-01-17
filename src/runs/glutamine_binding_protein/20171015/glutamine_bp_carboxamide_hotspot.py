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
print('finding hotspots')
nbrs = NearestNeighbors(metric='euclidean', radius=0.8, algorithm='kd_tree')
nbrs.fit(all_ifgs)
adj_mat = nbrs.radius_neighbors_graph(all_ifgs)
mems, cents = combs.Cluster.cluster_adj_mat_top_pc(adj_mat, pc=0.8)
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

# relvdm_carboxamide.find_hotspots(rmsd_cutoff=0.5)
# # print('# hotspots=', len(relvdm_carboxamide.hotspots))
# outdir = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20170912/'
# try:
#     os.makedirs(outdir)
# except:
#     pass
# relvdm_carboxamide.print_hotspots(outdir + 'carboxamide/', number=20)
# relvdm_carboxamide.rel_vdms_pickle = None
# relvdm_carboxamide.rel_vdm_ifg_coords = None
# relvdm_carboxamide.rel_vdm_bb_coords = None
# relvdm_carboxamide.rel_vdm_sc_coords = None
# relvdm_carboxamide.rel_vdm_tags = None
# relvdm_carboxamide.rel_vdm_resnames = None
# relvdm_carboxamide._ifgs = None
# relvdm_carboxamide._resn = None
# relvdm_carboxamide._type = None
# relvdm_carboxamide._vdm_tags = None
# relvdm_carboxamide._indices = None
# relvdm_carboxamide.ifg_dict = {'GLN': 'CG CD OE1 NE2'}
# # with open(outdir + 'relvdm_carboxamide.pickle', 'wb') as outfile:
# #     pickle.dump(relvdm_carboxamide, outfile)

