import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import os
import prody as pr
import time

try:
    output_dir_pdb = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20170724/carboxamide/2p4A_clash_cutoff_1Armsd/'
    os.makedirs(output_dir_pdb)
except:
    pass

ifg_dict = {'ASN': 'CB CG ND2 OD1', 'GLN': 'CG CD NE2 OE1'}
pdb_path = '/Users/npolizzi/Projects/combs/src/runs/glutamine_binding_protein/1wdn.pdb'
cb = combs.Comb(ifg_dict)
poi = pr.parsePDB(pdb_path)
resn_chid_pairs = zip([50, 67, 10, 115, 156, 13], ['A', 'A', 'A', 'A', 'A', 'A'])
rois = combs.get_rois(poi, resn_chid_pairs)
resn_chid_pairs = zip([50, 67, 10, 115, 156, 13], ['A', 'A', 'A', 'A', 'A', 'A'])
rois_phipsi = combs.get_rois_phi_psi(poi, resn_chid_pairs)
# print(rois_phipsi)
path_to_pickle = '/Users/npolizzi/Projects/combs/results/carboxamide/20170724/interactamers/pickle/'
interactamers = combs.load_interactamers(path_to_pickle)
# print(interactamers['BB_CO'].keys())
poi = poi.select('protein and not resnum 227')
t0=time.time()
g, hotspots, pdb_tags = combs.find_hotspot(rois, rois_phipsi, interactamers, poi, rmsd_cutoff=1,
                                           include_sc_clash=True, clash_cutoff=2.4)
print(time.time()-t0)
hotspots = list(hotspots)
print('Hotspots:')
print(sorted(hotspots, key=len, reverse=True))
print('Non-clashing pdbs used for hotspot search:')
print(list(len(val) for val in pdb_tags.values()))
# print(pdb_tags)
pdb_inpath = '/Users/npolizzi/Projects/combs/results/carboxamide/20170724/interactamers/pdbs/'
# pdb_outpath = '/Users/npolizzi/Projects/combs/results/gluatamine_binding_protein/20170726/'
combs.print_hotspot_pdbs(hotspots, rois, rois_phipsi, pdb_tags, interactamers, pdb_inpath, output_dir_pdb, top_spots=5)
