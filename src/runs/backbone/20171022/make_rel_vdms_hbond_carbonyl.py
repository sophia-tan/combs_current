import sys                                                                   
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs

ifg_dict = {'ANY': 'C O'}
csv_path = '/Users/npolizzi/Projects/combs/results/backbone/20171014/csv/'
outpath='/Users/npolizzi/Projects/combs/results/backbone/rel_vdms_hbond_carbonyl/20171022/'
cb = combs.Comb(ifg_dict)
combs.make_all_rel_vdms_hbond_partial_ifg(csv_path, outpath, cb, ifg_name='O', seq_distance=7)
