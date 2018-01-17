import sys                                                                   
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs

ifg_dict = {'ASP': 'CB CG OD1 OD2', 'GLU': 'CG CD OE1 OE2'}
csv_path = '/Users/npolizzi/Projects/combs/results/carboxylate/20170725/csv/'
outpath='/Users/npolizzi/Projects/combs/results/carboxylate/rel_vdms/20170830/'
cb = combs.Comb(ifg_dict)
combs.make_all_rel_vdms(csv_path, outpath, cb)
