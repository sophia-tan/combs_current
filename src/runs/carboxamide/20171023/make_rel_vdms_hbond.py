import sys                                                                   
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs

ifg_dict = {'ASN': 'CB CG OD1 ND2', 'GLN': 'CG CD OE1 NE2'}
csv_path = '/Users/npolizzi/Projects/combs/results/carboxamide/20170725/csv/'
outpath='/Users/npolizzi/Projects/combs/results/carboxamide/rel_vdms_hbond/20171023/'
cb = combs.Comb(ifg_dict)
combs.make_all_rel_vdms_hbond(csv_path, outpath, cb, seq_distance=7)
