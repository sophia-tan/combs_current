import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import prody as pr
# import time


pdb_path = '/Users/npolizzi/Projects/combs/src/runs/glutamine_binding_protein/1wdn_noligandH.pdb'
sample = combs.Sample()
sample.poi = pr.parsePDB(pdb_path)
sample.bs_residues = list(zip([68, 70, 157, 185], ['A', 'A', 'A', 'A']))
sample.set_rois()
sample.set_rois_phipsi()
sample.set_poi_clash_coords()
sample.set_roi_bb_coords()

relvdm_amino = combs.Rel_Vandermer('amino')
relvdm_amino.rel_vdm_path = '/Users/npolizzi/Projects/combs/results/amino/rel_vdms/20170829/'
relvdm_amino.load_rel_vdms_pickle(sample)
relvdm_amino.set_rel_vdm_bb_coords()
relvdm_amino.set_rois_rot_trans(sample)
relvdm_amino.set_rel_vdm_tags(sample)
# relvdm_amino.set_rel_vdm_resnames(sample)
relvdm_amino.move_rel_vdms(sample)
# i = 0
# for val in relvdm_amino.rel_vdm_tags[(50,'A')]['SC'].values():
#     i += val.shape[0]
# print(i)
relvdm_amino.remove_clash(sample)
# i = 0
# for val in relvdm_amino.rel_vdm_tags[(50,'A')]['SC'].values():
#     i += val.shape[0]
# print(i)
relvdm_amino.reshape_ifgs()
print('finding hotspots')
relvdm_amino.find_hotspots(rmsd_cutoff=0.8)
outdir = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20170830/amino/'
relvdm_amino.print_hotspots(outdir, number=25)

