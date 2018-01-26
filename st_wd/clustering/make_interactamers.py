# will output SC, N_CA, and C_O directories
import os
import sys                              
from IFG_dicts import *
sys.path.append('/home/gpu/Sophia/combs/src/')
import combs
from sys import argv

script, ifg = argv

outpath = '/home/gpu/Sophia/STcombs/20171118/%s/clusters/' %ifg
ifg_dir = '/home/gpu/Sophia/STcombs/20171118/%s/' %ifg

try:
    os.makedirs(outpath)
except:
    pass

ifg_sele = combs.apps.constants.ifg_sele_dict[ifg]
cb = combs.Comb(ifg_sele)

if ifg=='lonepair_imidazole':
    lonepair_imidazole=True
else:
    lonepair_imidazole=False

combs.make_bb_sc_rel_vdms(ifg_dir, cb,lonepair_imidazole=lonepair_imidazole)
