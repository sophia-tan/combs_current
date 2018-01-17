import sys                                                                   
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs

ifg_dict = {'LYS': 'CE NZ'}
csv_path = '/Users/npolizzi/Projects/combs/results/amino/20170725/csv/'
outpath='/Users/npolizzi/Projects/combs/results/amino/interactamers_hbond/20171027/'
cb = combs.Comb(ifg_dict)
combs.Interactamer.make_all_interactamers_hbond_partial_ifg(csv_path, outpath, cb, seq_distance=7)
