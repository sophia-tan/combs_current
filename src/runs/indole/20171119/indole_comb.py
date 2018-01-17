import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import os

try:
    output_dir_pdb = '/Users/npolizzi/Projects/combs/results/indole/20171119/vdM'
    os.makedirs(output_dir_pdb)
except:
    pass

try:
    output_dir_csv = '/Users/npolizzi/Projects/combs/results/indole/20171119/csv'
    os.makedirs(output_dir_csv)
except:
    pass

ifg_dict = {'TRP': 'CG CD1 NE1 CE2 CZ2 CH2 CZ3 CE3 CD2'}
kwargs = {'file_tag': 'indole',
          'input_dir_pdb': '/Users/npolizzi/Projects/combs/database/20170719/pdb/reduce',
          'output_dir_pdb': output_dir_pdb,
          'input_dir_dssp': '/Users/npolizzi/Projects/combs/database/20170719/dssp',
          'output_dir_csv': output_dir_csv,
          'path_to_pdb_chain_file': '/Users/npolizzi/Projects/combs/database/20170719/combs_pdb_list_vast_50seqid_blastclust.txt',
	  'path_to_reduce': '/Users/npolizzi/Applications/reduce',
	  'reduce': 'reduce.3.23.130521.macosx'
          }
combs.run_comb(ifg_dict, **kwargs)
