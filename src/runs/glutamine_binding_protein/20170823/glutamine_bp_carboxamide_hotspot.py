import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import prody as pr
# import time


pdb_path = '/Users/npolizzi/Projects/combs/src/runs/glutamine_binding_protein/1wdn_noligandH.pdb'
sample = combs.Sample()
sample.poi = pr.parsePDB(pdb_path)
sample.bs_residues = list(zip([50, 67, 10, 115, 156, 13], ['A', 'A', 'A', 'A', 'A', 'A']))
sample.set_rois()
sample.set_rois_phipsi()
sample.set_poi_clash_coords()
sample.set_roi_bb_coords()

relvdm_carbam = combs.Rel_Vandermer('carboxamide')
relvdm_carbam.rel_vdm_path = '/Users/npolizzi/Projects/combs/results/carboxamide/rel_vdms/20170822/'
relvdm_carbam.load_rel_vdms_pickle(sample)
relvdm_carbam.set_rel_vdm_bb_coords()
relvdm_carbam.set_rois_rot_trans(sample)
relvdm_carbam.set_rel_vdm_tags(sample)
relvdm_carbam.set_rel_vdm_resnames(sample)
relvdm_carbam.move_rel_vdms(sample)
relvdm_carbam.remove_clash(sample)

