import sys
sys.path.append('/netapp/home/nick.polizzi/combs/src/')
import combs
import prody as pr
import pickle
import os
import gzip
# import collections
import networkx as nx
# import time

pdb_path = '/netapp/home/nick.polizzi/combs/src/runs/glutamine_binding_protein/1wdn_noligandH.pdb'
sample = combs.Sample()
sample.poi = pr.parsePDB(pdb_path)
# sample.bs_residues = list(zip([70, 75, 119, 118, 50, 67, 10, 115, 156, 13, 68, 157, 185],
#                               ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A']))
sample.bs_residues = list(zip([10, 12, 13, 16, 17, 50, 53, 67, 68, 69, 70, 72, 75,
                               88, 115, 117, 118, 119, 120, 122, 139, 156, 157, 160,
                               183, 185],
                              ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
                               'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
                               'A', 'A']))
indir = '/netapp/home/nick.polizzi/combs/src/runs/glutamine_binding_protein/20170915/'
sample.set_rois()
sample.set_rois_phipsi()
sample.set_poi_clash_coords()
sample.set_roi_bb_coords()

relvdm_carboxylate = combs.Rel_Vandermer('carboxylate')
relvdm_carboxylate.rel_vdm_path = '/netapp/home/nick.polizzi/combs/results/carboxylate/rel_vdms/20170830/'
relvdm_carboxylate.load_rel_vdms_pickle(sample)
relvdm_carboxylate.set_rel_vdm_bb_coords()
relvdm_carboxylate.set_rois_rot_trans(sample)
relvdm_carboxylate.set_rel_vdm_tags(sample)
print('moving carboxylate vdMs')
relvdm_carboxylate.move_rel_vdms(sample)
print('removing clashing carboxylate vdMs')
relvdm_carboxylate.remove_clash(sample)
relvdm_carboxylate.reshape_ifgs()
print('finding carboxylate hotspots')
relvdm_carboxylate.find_hotspots(rmsd_cutoff=0.7)
relvdm_carboxylate.rel_vdms_pickle = None
relvdm_carboxylate.rel_vdm_ifg_coords = None
relvdm_carboxylate.rel_vdm_bb_coords = None
relvdm_carboxylate.rel_vdm_sc_coords = None
relvdm_carboxylate.rel_vdm_tags = None
relvdm_carboxylate.rel_vdm_resnames = None
relvdm_carboxylate._ifgs = None
relvdm_carboxylate._resn = None
relvdm_carboxylate._type = None
relvdm_carboxylate._vdm_tags = None
relvdm_carboxylate._indices = None
relvdm_carboxylate.ifg_dict = {'GLU': 'CG CD OE1 OE2'}

relvdm_carboxamide = combs.Rel_Vandermer('carboxamide')
relvdm_carboxamide.rel_vdm_path = '/netapp/home/nick.polizzi/combs/results/carboxamide/rel_vdms/20170830/'
relvdm_carboxamide.load_rel_vdms_pickle(sample)
relvdm_carboxamide.set_rel_vdm_bb_coords()
relvdm_carboxamide.set_rois_rot_trans(sample)
relvdm_carboxamide.set_rel_vdm_tags(sample)
print('moving carboxamide vdMs')
relvdm_carboxamide.move_rel_vdms(sample)
print('removing clashing carboxamide vdMs')
relvdm_carboxamide.remove_clash(sample)
relvdm_carboxamide.reshape_ifgs()
print('finding carboxamide hotspots')
relvdm_carboxamide.find_hotspots(rmsd_cutoff=0.7)
relvdm_carboxamide.rel_vdms_pickle = None
relvdm_carboxamide.rel_vdm_ifg_coords = None
relvdm_carboxamide.rel_vdm_bb_coords = None
relvdm_carboxamide.rel_vdm_sc_coords = None
relvdm_carboxamide.rel_vdm_tags = None
relvdm_carboxamide.rel_vdm_resnames = None
relvdm_carboxamide._ifgs = None
relvdm_carboxamide._resn = None
relvdm_carboxamide._type = None
relvdm_carboxamide._vdm_tags = None
relvdm_carboxamide._indices = None
relvdm_carboxamide.ifg_dict = {'GLN': 'CG CD OE1 NE2'}

relvdm_amino = combs.Rel_Vandermer('amino')
relvdm_amino.rel_vdm_path = '/netapp/home/nick.polizzi/combs/results/amino/rel_vdms/20170927/'
relvdm_amino.load_rel_vdms_pickle(sample)
relvdm_amino.set_rel_vdm_bb_coords()
relvdm_amino.set_rois_rot_trans(sample)
relvdm_amino.set_rel_vdm_tags(sample)
print('moving amino vdMs')
relvdm_amino.move_rel_vdms(sample)
print('removing clashing amino vdMs')
relvdm_amino.remove_clash(sample)
relvdm_amino.reshape_ifgs()
print('finding amino hotspots')
relvdm_amino.find_hotspots(rmsd_cutoff=0.7)
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
relvdm_amino.ifg_dict = {'LYS': 'CE NZ'}

relvdm_carboxamide.hotspot_subgraphs = list(
    sorted((comp for comp in nx.connected_component_subgraphs(relvdm_carboxamide.hotspot_graph) if len(comp) > 2),
           key=len, reverse=True))

relvdm_carboxylate.hotspot_subgraphs = list(
    sorted((comp for comp in nx.connected_component_subgraphs(relvdm_carboxylate.hotspot_graph) if len(comp) > 2),
           key=len, reverse=True))

relvdm_amino.hotspot_subgraphs = list(
    sorted((comp for comp in nx.connected_component_subgraphs(relvdm_amino.hotspot_graph) if len(comp) > 2),
           key=len, reverse=True))

outdir = '/netapp/home/nick.polizzi/combs/results/glutamine_binding_protein/20170927/output/'
try:
    os.makedirs(outdir)
except:
    pass

confs_indir = '/netapp/home/nick.polizzi/combs/src/runs/glutamine_binding_protein/confs/'
sample.ligand_conformers = [pr.parsePDB(confs_indir + file) for file in os.listdir(confs_indir)
                            if file in ['model0.pdb']]

sample.set_lig_ifg_dict(relvdm_carboxylate, 'CA C O OXT')
sample.set_lig_ifg_dict(relvdm_carboxamide, 'CG CD OE1 NE2')
sample.set_lig_ifg_dict(relvdm_amino, 'CA N')
sample.set_ligand_ifg_coords(relvdm_carboxylate)
sample.set_ligand_ifg_coords(relvdm_carboxamide)
sample.set_ligand_ifg_coords(relvdm_amino)
sample.set_rel_vdms([relvdm_carboxylate, relvdm_carboxamide, relvdm_amino])

def rev(name_pairs):
    return [((tup[0][1], tup[0][0]), tup[1]) for tup in name_pairs]

sample.name_pairs[('carboxylate', 'carboxamide')] = \
    [(('OXT', 'OE1'), 0.25), (('OXT', 'NE2'), 0.25), (('OXT', 'CG'), 0.3),
     (('O', 'OE1'), 0.25), (('O', 'NE2'), 0.25), (('O', 'CG'), 0.3),
     (('CA', 'OE1'), 0.3), (('CA', 'NE2'), 0.3), (('CA', 'CG'), 0.3)]
sample.name_pairs[('carboxamide', 'carboxylate')] = \
    rev(sample.name_pairs[('carboxylate', 'carboxamide')])

sample.name_pairs[('carboxylate', 'amino')] = \
    [(('OXT', 'CA'), 0.3), (('OXT', 'N'), 0.25),
     (('O', 'CA'), 0.3), (('O', 'N'), 0.25),
     (('CA', 'CA'), 0.3), (('CA', 'N'), 0.3)]
sample.name_pairs[('amino', 'carboxylate')] = \
    rev(sample.name_pairs[('carboxylate', 'amino')])

sample.name_pairs[('carboxamide', 'amino')] = \
    [(('OE1', 'CA'), 0.3), (('OE1', 'N'), 0.25),
     (('NE2', 'CA'), 0.3), (('NE2', 'N'), 0.25),
     (('CG', 'CA'), 0.3), (('CG', 'N'), 0.3)]
sample.name_pairs[('amino', 'carboxamide')] = \
    rev(sample.name_pairs[('carboxamide', 'amino')])

print('finding lig cst hotspots')
sample.find_ligand_cst_hotspots()
sample.set_protein_bb_coords(clash_cutoff=2.5)
sample.make_densities()
sample.set_lig_coords()
print('fitting lig confs to density')
sample.fit_lig_to_density()
if sample.poses:
    sample.pose_metrics(rmsd_cutoff=1.7)
    pscores = []
    ifgs = []
    num_sites = []
    poses = []
    for ind in sample.ranked_poses_by_score:
        pscore = sample.pose_scores[ind]
        ifg = sample.pose_num_satisfied_iFGs[ind]
        num_site = sample.pose_num_residues[ind]
        pose = ind
        if pscore > 1:
            pscores.append(pscore)
            ifgs.append(ifg)
            num_sites.append(num_site)
            poses.append(pose)
            for relvdm_key in sample.rel_vdms.keys():
                outdir_pdb = outdir + 'pose' + str(ind) + '/' + relvdm_key + '/'
                try:
                    os.makedirs(outdir_pdb)
                except:
                    pass
                sample.poses[ind].print_graph_pdbs(sample, relvdm_key, outdir_pdb)
            pr.writePDB(outdir + 'pose' + str(ind) + '/' + 'ligand.pdb', sample.poses[ind].ligand)
    if pscores:
        with open(outdir + 'gln_bp_pose_scores.txt', 'w') as outfile:
            outfile.write('pose scores, ' + ' '.join(str(s) for s in pscores) + '\n')
            outfile.write('number satisfied iFGs, ' + ' '.join(str(s) for s in ifgs) + '\n')
            outfile.write('number residue sites in pose, ' + ' '.join(str(s) for s in num_sites) + '\n')
            outfile.write('pose, ' + ' '.join(str(s) for s in poses) + '\n')

# if sample.pose_scores[sample.ranked_poses_by_score[0]] > 600:
sample.poi = None
sample.rois = None
with gzip.open(outdir + 'sample.pickle.gz', 'wb') as outfile:
    pickle.dump(sample, outfile)


