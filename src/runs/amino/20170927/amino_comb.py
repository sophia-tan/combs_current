import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import os

try:
    output_dir_pdb = '/Users/npolizzi/Projects/combs/results/amino/20170927/vdM'
    os.makedirs(output_dir_pdb)
except:
    pass

try:
    output_dir_csv = '/Users/npolizzi/Projects/combs/results/amino/20170927/csv'
    os.makedirs(output_dir_csv)
except:
    pass

ifg_dict = {'LYS': 'CE NZ'}
kwargs = {'file_tag': 'amino',
          'input_dir_pdb': '/Users/npolizzi/Projects/combs/database/20170719/pdb/reduce',
          'output_dir_pdb': output_dir_pdb,
          'input_dir_dssp': '/Users/npolizzi/Projects/combs/database/20170719/dssp',
          'output_dir_csv': output_dir_csv,
          'path_to_pdb_chain_file': '/Users/npolizzi/Projects/combs/database/20170719/combs_pdb_list_vast_50seqid_blastclust.txt',
	  'path_to_reduce': '/Users/npolizzi/Applications/hhhhreduce',
	  'reduce': 'reduce.3.23.130521.macosx'
          }
combs.run_comb(ifg_dict, **kwargs)
