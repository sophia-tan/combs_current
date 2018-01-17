import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import prody as pr
import pickle
import os
# import time


pdb_path = '/Users/npolizzi/Projects/combs/src/runs/dirhodium/template.pdb'
sample = combs.Sample()
sample.poi = pr.parsePDB(pdb_path)
sample.bs_residues = list(zip([1, 2, 5, 8, 9, 39, 40, 43, 44, 47, 1, 2, 5, 8, 9, 39, 40, 43, 44, 47],
                              ['A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'B',
                               'C', 'C', 'C', 'C', 'C', 'D', 'D', 'D', 'D', 'D']))
sample.set_rois()
sample.set_rois_phipsi()
sample.set_poi_clash_coords()
sample.set_roi_bb_coords()

relvdm_carboxylate = combs.Rel_Vandermer('carboxylate')
relvdm_carboxylate.rel_vdm_path = '/Users/npolizzi/Projects/combs/results/carboxylate/rel_vdms/20170830/'
relvdm_carboxylate.load_rel_vdms_pickle(sample, subset={'ALA', 'GLY', 'SER', 'THR'})
relvdm_carboxylate.set_rel_vdm_bb_coords()
relvdm_carboxylate.set_rois_rot_trans(sample)
relvdm_carboxylate.set_rel_vdm_tags(sample)
print('moving vdMs')
relvdm_carboxylate.move_rel_vdms(sample)
print('removing clashing vdMs')
relvdm_carboxylate.remove_clash(sample)
relvdm_carboxylate.reshape_ifgs()
print('finding hotspots')
relvdm_carboxylate.find_hotspots(rmsd_cutoff=0.8)
# print('# hotspots=', len(relvdm_carboxylate.hotspots))
outdir = '/Users/npolizzi/Projects/combs/results/dirhodium/20170922/'
try:
    os.makedirs(outdir)
except:
    pass
# relvdm_carboxylate.print_hotspots(outdir + 'carboxylate/', number=20)
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
print('done')
#with open(outdir + 'relvdm_carboxylate_subset.pickle', 'wb') as outfile:
#    pickle.dump(relvdm_carboxylate, outfile)

