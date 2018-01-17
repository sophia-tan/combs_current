import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import prody as pr
import pickle
import os
# import collections
import networkx as nx
# import time
import itertools


def pickle_load(file_path):
    with open(file_path, "rb") as f:
        return pickle.load(f)

indir = '/Users/npolizzi/Projects/combs/results/dirhodium/20170920/'
# relvdm_carboxamide = pickle_load(indir + 'relvdm_carboxamide.pickle')
relvdm_carboxylate1 = pickle_load(indir + 'relvdm_carboxylate.pickle')
relvdm_carboxylate1.name = 'carboxylate1'
relvdm_carboxylate2 = pickle_load(indir + 'relvdm_carboxylate.pickle')
relvdm_carboxylate2.name = 'carboxylate2'
relvdm_carboxylate3 = pickle_load(indir + 'relvdm_carboxylate.pickle')
relvdm_carboxylate3.name = 'carboxylate3'
relvdm_carboxylate4 = pickle_load(indir + 'relvdm_carboxylate.pickle')
relvdm_carboxylate4.name = 'carboxylate4'
# relvdm_amino = pickle_load(indir + 'relvdm_amino.pickle')
# sample = pickle_load(indir + 'sample.pickle')
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
# relvdm_carboxamide.hotspot_subgraphs = list(
#     sorted((comp for comp in nx.connected_component_subgraphs(relvdm_carboxamide.hotspot_graph) if len(comp) > 4),
#            key=len, reverse=True))
#
relvdm_carboxylate1.hotspot_subgraphs = list(
    sorted((comp for comp in nx.connected_component_subgraphs(relvdm_carboxylate1.hotspot_graph) if len(comp) > 10),
           key=len, reverse=True))

relvdm_carboxylate2.hotspot_subgraphs = list(
    sorted((comp for comp in nx.connected_component_subgraphs(relvdm_carboxylate2.hotspot_graph) if len(comp) > 10),
           key=len, reverse=True))

relvdm_carboxylate3.hotspot_subgraphs = list(
    sorted((comp for comp in nx.connected_component_subgraphs(relvdm_carboxylate3.hotspot_graph) if len(comp) > 10),
           key=len, reverse=True))

relvdm_carboxylate4.hotspot_subgraphs = list(
    sorted((comp for comp in nx.connected_component_subgraphs(relvdm_carboxylate4.hotspot_graph) if len(comp) > 10),
           key=len, reverse=True))

#
# relvdm_amino.hotspot_subgraphs = list(
#     sorted((comp for comp in nx.connected_component_subgraphs(relvdm_amino.hotspot_graph) if len(comp) > 4),
#            key=len, reverse=True))

confs_indir = '/Users/npolizzi/Projects/combs/src/runs/dirhodium/dirhod.pdb'
sample.ligand_conformers = [pr.parsePDB(confs_indir)]
# sample.set_ligand_ifg_correspondence(relvdm_carboxamide, 'CD CD NE2', 'CD NE2 NE2')
# sample.set_ligand_ifg_correspondence(relvdm_carboxylate1, 'C1 O1', 'CD OE2')
# sample.set_ligand_ifg_correspondence(relvdm_carboxylate2, 'C1D C1D O1D', 'CD OE2 OE2')
sample.set_ligand_ifg_correspondence(relvdm_carboxylate1, 'C1 C1 O2 O1', 'OE2 OE1 OE2 OE1')
sample.set_ligand_ifg_correspondence(relvdm_carboxylate2, 'C1D C1D O1D O2D', 'OE2 OE1 OE2 OE1')
sample.set_ligand_ifg_correspondence(relvdm_carboxylate3, 'C3 C3 O4D O3', 'OE2 OE1 OE2 OE1')
sample.set_ligand_ifg_correspondence(relvdm_carboxylate4, 'C3D C3D O3D O4', 'OE2 OE1 OE2 OE1')
# sample.set_ligand_ifg_correspondence(relvdm_amino, 'CA CA N', 'CE NZ NZ')
# sample.set_ligand_ifg_coords(relvdm_carboxamide)
sample.set_ligand_ifg_coords(relvdm_carboxylate1)
sample.set_ligand_ifg_coords(relvdm_carboxylate2)
sample.set_ligand_ifg_coords(relvdm_carboxylate3)
sample.set_ligand_ifg_coords(relvdm_carboxylate4)
# sample.set_ligand_ifg_coords(relvdm_amino)
sample.set_rel_vdms([relvdm_carboxylate1, relvdm_carboxylate2, relvdm_carboxylate3, relvdm_carboxylate4])
print('finding lig cst hotspots')
sample.find_ligand_cst_hotspots(tol=.3)
# sample.hotspot_cliques[0] = [[j] for j in itertools.zip_longest(range(len(relvdm_carboxylate.hotspot_subgraphs)),
#                                                                 ['carboxylate'], fillvalue='carboxylate')]
sample.set_protein_bb_coords(clash_cutoff=2.5)
sample.make_densities()
# sample.set_lig_ifg_dict(relvdm_carboxamide, 'CG CD OE1 NE2')
sample.set_lig_ifg_dict(relvdm_carboxylate1, 'C2 C1 O2 O1')
sample.set_lig_ifg_dict(relvdm_carboxylate2, 'C2D C1D O1D O2D')
sample.set_lig_ifg_dict(relvdm_carboxylate3, 'C4 C3 O4D O3')
sample.set_lig_ifg_dict(relvdm_carboxylate4, 'C4D C3D O3D O4')
# sample.set_lig_ifg_dict(relvdm_amino, 'C CA N')  # NOTE: could have chosen C CA N here
sample.set_lig_coords()
print('fitting lig confs to density')
sample.fit_lig_to_density()





