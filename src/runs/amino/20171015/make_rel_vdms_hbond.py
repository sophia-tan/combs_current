import sys                                                                   
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs

ifg_dict = {'LYS': 'CE NZ'}
csv_path = '/Users/npolizzi/Projects/combs/results/amino/20170725/csv/'
outpath='/Users/npolizzi/Projects/combs/results/amino/rel_vdms_hbond/20171015/'
cb = combs.Comb(ifg_dict)
combs.make_all_rel_vdms_hbond(csv_path, outpath, cb, seq_distance=7)
