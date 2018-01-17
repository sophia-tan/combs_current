import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs

ifg_dict = {1: {'ANY': 'C O CA'}, 2: {'PRO': 'N'}}
kwargs = {'file_tag': 'preproline_carbonyl',
          'ifg_seq_str': '.P',
          'input_dir_pdb': '../20170702/pdb',
          'output_dir_pdb': '../../results/20170703/preproline_carbonyl/vdM',
          'input_dir_dssp': '../20170702/dssp',
          'output_dir_csv': '../../results/20170703/preproline_carbonyl/csv',
          'path_to_pdb_chain_file': '../20170702/clusterReps_rcsb_1p8_30seqid_r25_28jan2017_trimmed_fmt_blastclust_pt_bF_S12_subst.txt'
          }
combs.run_comb(ifg_dict, **kwargs)