from os import listdir , makedirs
from os.path import isfile, join
from sys import argv
import pickle as pkl, numpy as np
import math
from sklearn.neighbors import NearestNeighbors
import sys

script, ifg = argv

num_atoms = {}

lookup_dir = '/home/gpu/Sophia/combs/st_wd/Lookups/'
outpath = '../Lookups/expected_NN/' 
try:
    makedirs(outpath)
except:
    pass

def neighb_exp(ifg,typ,resn,lookup_dir):
    '''helper function to calculate expected # of NNs
    Takes all the vdms of a certain AA type that interact with the iFG, calculates
    their NNs, and then takes avg'''
    NNs = []
    path = '/home/gpu/Sophia/STcombs/20171118/%s/clusters/%s/pickle/%s_rel_vdms.pickle'%(ifg,typ,resn)
    ifgcoords = pkl.load(open(path,'rb'))[:,6]
    for coord in ifgcoords:
        coord = coord.reshape(1,-1)
        nnfit = pkl.load(open(lookup_dir+'fitted_NN/%s/NNfit_%s_with_%s_%s.pkl'%(ifg,ifg,resn,typ),'rb'))
        number_nn = len(nnfit.radius_neighbors(coord,return_distance=False)[0])
        NNs.append(number_nn)
    return np.mean(NNs)

expdict = {} # key = (typ,resn), val = exp_nn
for typ in ['SC', 'N_CA', 'C_O']:
    for resn in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR','VAL']:
        try:
            exp = neighb_exp(ifg,typ,resn,lookup_dir)
            expdict[(typ,resn)] = exp
        except:
            pass
pkl.dump(expdict, open(outpath+'expNN_%s_ifg.pkl' % ifg, 'wb'))

