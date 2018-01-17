import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import prody as pr
# import time


pdb_path = '/Users/npolizzi/Projects/combs/src/runs/glutamine_binding_protein/1wdn_noligandH.pdb'
sample = combs.Sample()
sample.poi = pr.parsePDB(pdb_path)
sample.bs_residues = list(zip([50, 67, 10, 115, 156, 13, 68, 70, 157, 185],
                              ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A']))
sample.set_rois()
sample.set_rois_phipsi()
sample.set_poi_clash_coords()
sample.set_roi_bb_coords()

relvdm_carboxamide = combs.Rel_Vandermer('carboxamide')
relvdm_carboxamide.rel_vdm_path = '/Users/npolizzi/Projects/combs/results/carboxamide/rel_vdms/20170830/'
relvdm_carboxamide.load_rel_vdms_pickle(sample)
relvdm_carboxamide.set_rel_vdm_bb_coords()
relvdm_carboxamide.set_rois_rot_trans(sample)
relvdm_carboxamide.set_rel_vdm_tags(sample)
print('moving vdMs')
relvdm_carboxamide.move_rel_vdms(sample)
print('removing clashing vdMs')
relvdm_carboxamide.remove_clash(sample)
relvdm_carboxamide.reshape_ifgs()
print('finding hotspots')
relvdm_carboxamide.find_hotspots(rmsd_cutoff=0.6)
outdir = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20170908/carboxamide/'
relvdm_carboxamide.print_hotspots(outdir, number=10)

