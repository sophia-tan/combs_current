import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import prody as pr
import pickle
import os
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

relvdm_amino = combs.Rel_Vandermer('amino')
relvdm_amino.rel_vdm_path = '/Users/npolizzi/Projects/combs/results/amino/rel_vdms/20170830/'
relvdm_amino.load_rel_vdms_pickle(sample)
relvdm_amino.set_rel_vdm_bb_coords()
relvdm_amino.set_rois_rot_trans(sample)
relvdm_amino.set_rel_vdm_tags(sample)
print('moving vdMs')
relvdm_amino.move_rel_vdms(sample)
print('removing clashing vdMs')
relvdm_amino.remove_clash(sample)
relvdm_amino.reshape_ifgs()
print('finding hotspots')
relvdm_amino.find_hotspots(rmsd_cutoff=0.6)
# print('# hotspots=', len(relvdm_amino.hotspots))
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
relvdm_amino.ifg_dict = {'LYS': 'CD CE NZ'}
outdir = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20170912/'
try:
    os.makedirs(outdir)
except:
    pass
# relvdm_amino.print_hotspots(outdir, number=20)
with open(outdir + 'relvdm_amino.pickle', 'wb') as outfile:
    pickle.dump(relvdm_amino, outfile)

