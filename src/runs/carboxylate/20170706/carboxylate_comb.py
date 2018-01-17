import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs

ifg_dict = {'ASP': 'CB CG OD2 OD1', 'GLU': 'CG CD OE2 OE1'}
kwargs = {'file_tag': 'carboxylate',
          'input_dir_pdb': '/Users/npolizzi/Projects/combs/database/20170702/pdb',
          'output_dir_pdb': '/Users/npolizzi/Projects/combs/results/carboxylate/20170706/vdM',
          'input_dir_dssp': '/Users/npolizzi/Projects/combs/database/20170702/dssp',
          'output_dir_csv': '/Users/npolizzi/Projects/combs/results/carboxylate/20170706/csv',
          'path_to_pdb_chain_file': '/Users/npolizzi/Projects/combs/database/20170702/clusterReps_rcsb_1p8_30seqid_r25_28jan2017_trimmed_fmt_blastclust_pt_bF_S12_subst.txt'
          }
combs.run_comb(ifg_dict, **kwargs)
