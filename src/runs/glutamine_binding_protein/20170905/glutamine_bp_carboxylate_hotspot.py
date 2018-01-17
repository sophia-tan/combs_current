import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import prody as pr
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

relvdm_carboxylate = combs.Rel_Vandermer('carboxylate')
relvdm_carboxylate.rel_vdm_path = '/Users/npolizzi/Projects/combs/results/carboxylate/rel_vdms/20170830/'
relvdm_carboxylate.load_rel_vdms_pickle(sample)
relvdm_carboxylate.set_rel_vdm_bb_coords()
relvdm_carboxylate.set_rois_rot_trans(sample)
relvdm_carboxylate.set_rel_vdm_tags(sample)
print('moving vdMs')
relvdm_carboxylate.move_rel_vdms(sample)
print('removing clashing vdMs')
relvdm_carboxylate.remove_clash(sample)
relvdm_carboxylate.reshape_ifgs()
print('finding hotspots')
relvdm_carboxylate.find_hotspots(rmsd_cutoff=0.6)
outdir = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20170905/carboxylate/'
relvdm_carboxylate.print_hotspots(outdir, number=10)

