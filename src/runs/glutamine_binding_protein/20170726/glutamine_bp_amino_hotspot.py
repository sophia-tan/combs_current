import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import os
import prody as pr
import time

try:
    output_dir_pdb = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20170726/amino/bb_2p4A_clash_cutoff_0p5Armsd/'
    os.makedirs(output_dir_pdb)
except:
    pass

ifg_dict = {'LYS': 'CD CE NZ'}
pdb_path = '/Users/npolizzi/Projects/combs/src/runs/glutamine_binding_protein/1wdn.pdb'
cb = combs.Comb(ifg_dict)
poi = pr.parsePDB(pdb_path)
resn_chid_pairs = zip([68, 70, 157, 185], ['A', 'A', 'A', 'A'])
rois = combs.get_rois(poi, resn_chid_pairs)
resn_chid_pairs = zip([68, 70, 157, 185], ['A', 'A', 'A', 'A'])
rois_phipsi = combs.get_rois_phi_psi(poi, resn_chid_pairs)
# print(rois_phipsi)
path_to_pickle = '/Users/npolizzi/Projects/combs/results/amino/rel_vdms/20170725/pickle/'
interactamers = combs.load_rel_vdms(path_to_pickle)
# print(interactamers['BB_CO'].keys())
poi = poi.select('protein and not resnum 227')
t0=time.time()
g, hotspots, pdb_tags = combs.find_hotspot_bb(rois, interactamers, poi, rmsd_cutoff=1,
                                           include_sc_clash=False, clash_cutoff=2.4)
print(time.time()-t0)
hotspots = list(hotspots)
print('Hotspots:')
print(sorted(hotspots, key=len, reverse=True))
print('Non-clashing pdbs used for hotspot search:')
print(list(len(val) for val in pdb_tags.values()))
# print(pdb_tags)
pdb_inpath = '/Users/npolizzi/Projects/combs/results/amino/rel_vdms/20170725/pdbs/'
# pdb_outpath = '/Users/npolizzi/Projects/combs/results/gluatamine_binding_protein/20170726/amino/'
combs.print_hotspot_pdbs_bb(hotspots, rois, pdb_tags, interactamers, pdb_inpath, output_dir_pdb, top_spots=30)
