from os import listdir , makedirs
from os.path import isfile, join
from sys import argv
import pickle as pkl
import math
from sklearn.neighbors import NearestNeighbors

script, ifg = argv

num_atoms = {}
num_atoms['backboneCO'] = 2
num_atoms['serCOH'] = 2
num_atoms['thrCOH'] = 2
num_atoms['tyrCOH'] = 2
num_atoms['indole'] = 3
num_atoms['imidazole'] = 5
num_atoms['lonepair_imidazole'] = 5
num_atoms['amino'] = 3
num_atoms['carboxamide'] = 4
num_atoms['carboxylate'] = 4
num_atoms['guanidino'] = 4
num_atoms['phenyl'] = 6

outpath = '../Lookups/fitted_NN/%s/' %ifg
try:
    makedirs(outpath)
except:
    pass

for typ in ['SC', 'N_CA', 'C_O']:
    path = '/home/gpu/Sophia/STcombs/20171118/%s/clusters/%s/pickle/'%(ifg,typ)
    onlyfiles = [f for f in listdir(path) if isfile(join(path,f))]
    for pklfile in onlyfiles:
        pklf = pkl.load(open(path+pklfile, 'rb'))
        ifgcoords = pklf[:,6]
        flat = [x.flatten() for x in ifgcoords]
        dist = 1.4/math.sqrt(1/num_atoms[ifg])
        nbrs = NearestNeighbors(n_neighbors=1,metric='euclidean',radius=dist)
        vdmname = pklfile.split('_')[0]
        x = nbrs.fit(flat)

        pkl.dump(x, open(outpath+'NNfit_%s_with_%s_%s.pkl' %(ifg,vdmname,typ), 'wb'))

