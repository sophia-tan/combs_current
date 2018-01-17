import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import prody as pr
import pickle


inpath = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20171016/'
with open(inpath + 'relvdm_amino.pickle', 'rb') as infile:
    relvdm_amino = pickle.load(infile)
with open(inpath + 'relvdm_carboxamide.pickle', 'rb') as infile:
    relvdm_carboxamide = pickle.load(infile)
with open(inpath + 'relvdm_carboxylate.pickle', 'rb') as infile:
    relvdm_carboxylate = pickle.load(infile)

pdb_path = '/Users/npolizzi/Projects/combs/src/runs/glutamine_binding_protein/'
sample = combs.Sample()
sample.poi = pr.parsePDB(pdb_path + '1wdn_noligandH.pdb')
sample.bs_residues = list(zip([70, 75, 119, 118, 50, 67, 10, 115, 156, 13, 68, 157, 185],
                              ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A']))
sample.set_rois()
sample.set_rois_phipsi()
sample.set_poi_clash_coords()
sample.set_roi_bb_coords()
sample.ligand_conformers = [pr.parsePDB(pdb_path + 'confs/model0.pdb')]

relvdm_amino.lig_ifg_correspondence = 'CA N'
relvdm_carboxylate.lig_ifg_correspondence = 'CA C O OXT'
relvdm_carboxamide.lig_ifg_correspondence = 'CG CD OE1 NE2'

relvdms = [relvdm_amino, relvdm_carboxylate, relvdm_carboxamide]

print('finding cluster pairs')
sample.find_ligand_cst_cliques(relvdms)

outpath = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20171017/'
with open(outpath + 'relvdm_amino.pickle', 'wb') as outfile:
    pickle.dump(relvdm_amino, outfile)
with open(outpath + 'relvdm_carboxamide.pickle', 'wb') as outfile:
    pickle.dump(relvdm_carboxamide, outfile)
with open(outpath + 'relvdm_carboxylate.pickle', 'wb') as outfile:
    pickle.dump(relvdm_carboxylate, outfile)
with open(outpath + 'sample.pickle', 'wb') as outfile:
    pickle.dump(sample, outfile)

sample.set_protein_bb_coords(clash_cutoff=2.5)
print('making cluster densities')
sample.make_cluster_densities(relvdms)
print('fitting ligand to density')
sample.fit_lig_to_cluster_density(relvdms)