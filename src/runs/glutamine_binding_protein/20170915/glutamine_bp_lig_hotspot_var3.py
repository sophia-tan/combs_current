import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import prody as pr
import pickle
import collections
import networkx as nx
# import time


def pickle_load(file_path):
    with open(file_path, "rb") as f:
        return pickle.load(f)

pdb_path = '/Users/npolizzi/Projects/combs/src/runs/glutamine_binding_protein/1wdn_noligandH.pdb'
sample = combs.Sample()
sample.poi = pr.parsePDB(pdb_path)
sample.bs_residues = list(zip([70, 75, 119, 118, 50, 67, 10, 115, 156, 13, 68, 157, 185],
                              ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A']))
# sample.set_rois()
# sample.set_rois_phipsi()
# sample.set_poi_clash_coords()
# sample.set_roi_bb_coords()
# sample.lig_ifg_corr = collections.defaultdict(dict)
# sample.lig_ifg_coords = collections.defaultdict(dict)
# sample.hotspot_cliques = {}

indir = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20170912/'
relvdm_carboxamide = pickle_load(indir + 'relvdm_carboxamide.pickle')
relvdm_carboxylate = pickle_load(indir + 'relvdm_carboxylate.pickle')
relvdm_amino = pickle_load(indir + 'relvdm_amino.pickle')

# relvdm_carboxamide.hotspot_subgraphs = list(
#     sorted(nx.connected_component_subgraphs(relvdm_carboxamide.hotspot_graph), key=len, reverse=True))
#
# relvdm_carboxylate.hotspot_subgraphs = list(
#     sorted(nx.connected_component_subgraphs(relvdm_carboxylate.hotspot_graph), key=len, reverse=True))
#
# relvdm_amino.hotspot_subgraphs = list(
#     sorted(nx.connected_component_subgraphs(relvdm_amino.hotspot_graph), key=len, reverse=True))

relvdm_carboxamide.hotspot_subgraphs = list(
    sorted((comp for comp in nx.connected_component_subgraphs(relvdm_carboxamide.hotspot_graph) if len(comp) > 4),
           key=len, reverse=True))

relvdm_carboxylate.hotspot_subgraphs = list(
    sorted((comp for comp in nx.connected_component_subgraphs(relvdm_carboxylate.hotspot_graph) if len(comp) > 4),
           key=len, reverse=True))

relvdm_amino.hotspot_subgraphs = list(
    sorted((comp for comp in nx.connected_component_subgraphs(relvdm_amino.hotspot_graph) if len(comp) > 4),
           key=len, reverse=True))

sample.ligand_conformers = [pr.parsePDB('/Users/npolizzi/Desktop/gln227_mm_good.pdb')]
sample.set_ligand_ifg_correspondence(relvdm_carboxamide, 'CD CD NE2', 'CD NE2 NE2')
sample.set_ligand_ifg_correspondence(relvdm_carboxylate, 'C C O', 'CD OE2 OE2')
sample.set_ligand_ifg_correspondence(relvdm_amino, 'CA CA N', 'CE NZ NZ')
sample.set_ligand_ifg_coords(relvdm_carboxamide)
sample.set_ligand_ifg_coords(relvdm_carboxylate)
sample.set_ligand_ifg_coords(relvdm_amino)
sample.set_rel_vdms([relvdm_carboxamide, relvdm_amino, relvdm_carboxylate])
sample.find_ligand_cst_hotspots(tol=0.2)

sample.set_protein_bb_coords(clash_cutoff=2.5)
sample.make_densities()
sample.set_lig_ifg_dict(relvdm_carboxamide, 'CG CD OE1 NE2')
sample.set_lig_ifg_dict(relvdm_carboxylate, 'CA C OXT O')
sample.set_lig_ifg_dict(relvdm_amino, 'CB CA N')  # NOTE: could have chosen C CA N here
sample.set_lig_coords()
sample.fit_lig_to_density()
sample.set_results()
sample.set_sorted_results_ind()


