import sys
sys.path.append('/netapp/home/nick.polizzi/combs/src/')
import combs
import prody as pr
import pickle
import os
import networkx as nx
#import itertools
from copy import deepcopy
import gzip


def main():
    pdb_path = sys.argv[1]
    sample = combs.Sample()
    sample.poi = pr.parsePDB(pdb_path)
    sample.bs_residues = list(zip([7, 10, 11, 14, 7, 10, 11, 14, 17, 7, 10, 11, 14, 7, 10, 11, 14, 17],
                                  ['A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'B',
                                   'C', 'C', 'C', 'C', 'D', 'D', 'D', 'D', 'D']))
    sample.set_rois()
    sample.set_rois_phipsi()
    sample.set_poi_clash_coords()
    sample.set_roi_bb_coords()

    relvdm_carboxamide = combs.Rel_Vandermer('carboxamide')
    relvdm_carboxamide.rel_vdm_path = \
        '/netapp/home/nick.polizzi/combs/results/carboxamide/rel_vdms/20170830/'
    relvdm_carboxamide.load_rel_vdms_pickle(sample, subset={'ALA', 'GLY', 'SER', 'THR', 'TYR', 'TRP', 'ARG',
                                                            'ASN', 'GLN', 'HIS', 'PHE', 'LEU', 'VAL', 'ILE'})
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
    
    relvdm_preproline_carbonyl1 = combs.Rel_Vandermer('preproline_carbonyl1')
    relvdm_preproline_carbonyl1.rel_vdm_path = \
        '/netapp/home/nick.polizzi/combs/results/preproline_carbonyl/rel_vdms/20170927/'
    relvdm_preproline_carbonyl1.load_rel_vdms_pickle(sample, subset={'ALA', 'GLY', 'SER', 'THR', 'TYR', 'TRP', 'ARG',
                                                                     'ASN', 'GLN', 'HIS', 'PHE', 'LEU', 'VAL', 'ILE'})
    relvdm_preproline_carbonyl1.set_rel_vdm_bb_coords()
    relvdm_preproline_carbonyl1.set_rois_rot_trans(sample)
    relvdm_preproline_carbonyl1.set_rel_vdm_tags(sample)
    print('moving vdMs')
    relvdm_preproline_carbonyl1.move_rel_vdms(sample)
    print('removing clashing vdMs')
    relvdm_preproline_carbonyl1.remove_clash(sample)
    relvdm_preproline_carbonyl1.reshape_ifgs()
    print('finding hotspots')
    relvdm_preproline_carbonyl1.find_hotspots(rmsd_cutoff=0.7)
    
    relvdm_preproline_carbonyl1.rel_vdms_pickle = None
    relvdm_preproline_carbonyl1.rel_vdm_ifg_coords = None
    relvdm_preproline_carbonyl1.rel_vdm_bb_coords = None
    relvdm_preproline_carbonyl1.rel_vdm_sc_coords = None
    relvdm_preproline_carbonyl1.rel_vdm_tags = None
    relvdm_preproline_carbonyl1.rel_vdm_resnames = None
    relvdm_preproline_carbonyl1._ifgs = None
    relvdm_preproline_carbonyl1._resn = None
    relvdm_preproline_carbonyl1._type = None
    relvdm_preproline_carbonyl1._vdm_tags = None
    relvdm_preproline_carbonyl1._indices = None
    relvdm_preproline_carbonyl1.ifg_dict = {'ANY': 'C O CA N'}

    # with open(outdir + 'relvdm_preproline_carbonyl_subset.pickle', 'wb') as outfile:
    #     pickle.dump(relvdm_preproline_carbonyl1, outfile)

    relvdm_preproline_carbonyl2 = deepcopy(relvdm_preproline_carbonyl1)
    relvdm_preproline_carbonyl2.name = 'preproline_carbonyl2'

    relvdm_carboxamide.hotspot_subgraphs = list(
        sorted((comp for comp in nx.connected_component_subgraphs(relvdm_carboxamide.hotspot_graph) if len(comp) > 10),
               key=len, reverse=True))

    relvdm_preproline_carbonyl1.hotspot_subgraphs = list(
        sorted((comp for comp in nx.connected_component_subgraphs(relvdm_preproline_carbonyl1.hotspot_graph) if len(comp) > 10),
               key=len, reverse=True))

    relvdm_preproline_carbonyl2.hotspot_subgraphs = list(
        sorted((comp for comp in nx.connected_component_subgraphs(relvdm_preproline_carbonyl2.hotspot_graph) if len(comp) > 10),
               key=len, reverse=True))

    outdir = '/netapp/home/nick.polizzi/combs/results/apixaban/20170927/subset_TSAGYWNQHFLVIR/output/' \
             + pdb_path.split('/')[-1].split('.')[0] + '/'
    try:
        os.makedirs(outdir)
    except:
        pass
    
    confs_indir = '/netapp/home/nick.polizzi/combs/src/runs/apixaban/apix.pdb'
    sample.ligand_conformers = [pr.parsePDB(confs_indir)]
    sample.set_lig_ifg_dict(relvdm_carboxamide, 'C10 C11 O1 N3')
    sample.set_lig_ifg_dict(relvdm_preproline_carbonyl1, 'C8 O3 C13 N5')
    sample.set_lig_ifg_dict(relvdm_preproline_carbonyl2, 'C19 O2 C23 N2')
    sample.set_ligand_ifg_coords(relvdm_carboxamide)
    sample.set_ligand_ifg_coords(relvdm_preproline_carbonyl1)
    sample.set_ligand_ifg_coords(relvdm_preproline_carbonyl2)
    sample.set_rel_vdms([relvdm_carboxamide, relvdm_preproline_carbonyl1, relvdm_preproline_carbonyl2])

    def rev(name_pairs):
        return [((tup[0][1], tup[0][0]), tup[1]) for tup in name_pairs]

    sample.name_pairs[('carboxamide', 'preproline_carbonyl1')] = \
        [(('O1', 'N5'), 0.25), (('O1', 'C13'), 0.3), (('O1', 'O3'), 0.25),
         (('N3', 'N5'), 0.25), (('N3', 'C13'), 0.3), (('N3', 'O3'), 0.25),
         (('C10', 'N5'), 0.3), (('C10', 'C13'), 0.3), (('C10', 'O3'), 0.3)]
    sample.name_pairs[('preproline_carbonyl1', 'carboxamide')] = \
        rev(sample.name_pairs[('carboxamide', 'preproline_carbonyl1')])

    sample.name_pairs[('carboxamide', 'preproline_carbonyl2')] = \
        [(('O1', 'N2'), 0.25), (('O1', 'C23'), 0.3), (('O1', 'O2'), 0.25),
         (('N3', 'N2'), 0.25), (('N3', 'C23'), 0.3), (('N3', 'O2'), 0.25),
         (('C10', 'N2'), 0.3), (('C10', 'C23'), 0.3), (('C10', 'O2'), 0.3)]
    sample.name_pairs[('preproline_carbonyl2', 'carboxamide')] = \
        rev(sample.name_pairs[('carboxamide', 'preproline_carbonyl2')])

    sample.name_pairs[('preproline_carbonyl1', 'preproline_carbonyl2')] = \
        [(('N5', 'N2'), 0.25), (('N5', 'C23'), 0.25), (('N5', 'O2'), 0.3),
         (('C13', 'N2'), 0.25), (('C13', 'C23'), 0.25), (('C13', 'O2'), 0.3),
         (('O3', 'N2'), 0.3), (('O3', 'C23'), 0.3), (('O3', 'O2'), 0.3)]
    sample.name_pairs[('preproline_carbonyl2', 'preproline_carbonyl1')] = \
        rev(sample.name_pairs[('preproline_carbonyl1', 'preproline_carbonyl2')])

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
                keys = set(sample.rel_vdms.keys()) - {'preproline_carbonyl'}
                for relvdm_key in keys:
                    if relvdm_key[:3] == 'pre':
                        sample.rel_vdms['preproline_carbonyl'] = sample.rel_vdms[relvdm_key]
                        sample.rel_vdms['preproline_carbonyl'].name = 'preproline_carbonyl'
                        sample.rel_vdms['preproline_carbonyl'].rel_vdm_path = \
                            '/netapp/home/nick.polizzi/combs/results/preproline_carbonyl/rel_vdms/20170830/'
                        outdir_pdb = outdir + 'pose' + str(ind) + '/' + relvdm_key + '/'
                        try:
                            os.makedirs(outdir_pdb)
                        except:
                            pass
                        sample.poses[ind].subgraphs['preproline_carbonyl'] = sample.poses[ind].subgraphs[relvdm_key]
                        sample.poses[ind].print_graph_pdbs(sample, 'preproline_carbonyl', outdir_pdb)
                        pr.writePDB(outdir + 'pose' + str(ind) + '/' + 'ligand.pdb', sample.poses[ind].ligand)
                    else:
                        outdir_pdb = outdir + 'pose' + str(ind) + '/' + relvdm_key + '/'
                        try:
                            os.makedirs(outdir_pdb)
                        except:
                            pass
                        sample.poses[ind].print_graph_pdbs(sample, relvdm_key, outdir_pdb)
                        pr.writePDB(outdir + 'pose' + str(ind) + '/' + 'ligand.pdb', sample.poses[ind].ligand)
        if pscores:
            with open(outdir + pdb_path.split('/')[-1].split('.')[0] + '_pose_scores.txt', 'w') as outfile:
                outfile.write(sys.argv[1].split('/')[-1] + ', ' + ' '.join(str(s) for s in pscores) + '\n')
                outfile.write('number satisfied iFGs, ' + ' '.join(str(s) for s in ifgs) + '\n')
                outfile.write('number residue sites in pose, ' + ' '.join(str(s) for s in num_sites) + '\n')
                outfile.write('pose, ' + ' '.join(str(s) for s in poses) + '\n')


    if sample.pose_scores[sample.ranked_poses_by_score[0]] > 1000:
        sample.poi = None
        sample.rois = None
        with gzip.open(outdir + 'sample_TSAGYWNQHFLVI_' + pdb_path.split('/')[-1].split('.')[0] + '.pickle.gz', 'wb') as outfile:
            pickle.dump(sample, outfile)


if __name__ == "__main__":
    main()







