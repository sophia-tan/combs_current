# generate lookup table for burial (min dist to convex hull) in nrPDB database
import sys, os, subprocess
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis
import numpy as np, pandas as pd, pickle as pkl
from numpy import concatenate, array_equal, dot, cross, linalg, abs
from scipy.spatial import ConvexHull, Delaunay
from prody import *
import traceback

# directory of PDBs in database that have been put through reduce
pdb_dir = '/home/gpu/Sophia/STcombs/20171118/database/reduce/' 

####### functions for convex hull ######
def pnt_in_cvex_hull(hull, pnt):
    '''
    Checks if `pnt` is inside the convex hull.
    `hull` -- a QHull ConvexHull object
    `pnt` -- point array of shape (3,)
    '''
    new_hull = ConvexHull(concatenate((hull.points, [pnt])))
    if array_equal(new_hull.vertices, hull.vertices):
        return 1
    return -1

def normal(ps):
    v1 = ps[1] - ps[0]
    v2 = ps[2] - ps[0]
    crossprod = cross(v1, v2)
    return crossprod / linalg.norm(crossprod)

def distance(p, ps):
    return abs(dot(normal(ps), p - ps[0]))
####### functions for convex hull ######

# compile lists of pdbs and chains
pdbnames = []
chains = []

with open('combs_vast_50seqid_blastclust_07252017.txt') as fo:
    for line in fo:
        pdbnames.append(line[:4])
        chains.append(line[4])
zipped = zip(pdbnames, chains)

# create table 
col = ['pdbcode', 'chain', 'residuenum', 'residuename', 'burialdist']
list_rows = []
for pdb, chain in zipped:
    try:
        pdbfile = pdb_dir+pdb+'H.pdb'
        orig_pdb = parsePDB(pdbfile, altloc='A', model=1) # to preserve resnum and resname
        orig_pdb = orig_pdb.select('chain %s and name CA'%chain) # get CA instead of CB bc some glycines
        
        # use Rosetta's fixbb to mutate to all-ala to get C-beta of every residue
        parsed = parsePDB(pdb_dir+pdb+'H_0001.pdb', altloc='A', model=1) 
        
        # select c-betas from all ala and get min dist to convex hull
        parsed = parsed.select('chain %s and name CB' % chain)
        cb_coords = parsed.getCoords()
        resnums = orig_pdb.getResnums()
        resnames = orig_pdb.getResnames()
        tri = Delaunay(cb_coords)
        min_distances = []
        for p in tri.points:
            distances = [distance(p, tri.points[tri.convex_hull[i]])
                         for i in range(len(tri.convex_hull))]
            min_distances.append(min(distances))
        resnums_distances = list(zip(resnums, resnames, min_distances))
        
        for num, name, dist in resnums_distances:
            row = [pdb, chain, num, name, np.round(dist,2)]
            list_rows.append(row)
        print(pdb)
    
    except Exception:
        print(traceback)
        print(pdb, 'pdb could not be parsed')
    

df = pd.DataFrame(list_rows, columns=col)
pkl.dump(df, open('../database_mindistconvexhull.pkl','wb'))

#########################################################


#total_AAi = sum([val for key, val in raw_counts_lookup.items()])
#print('total AA in database: ', total_AAi)
#
#AAi_freq = {} # lookup for frequencies
#for key, val in raw_counts_lookup.items():
#    AAi_freq[key] = val / total_AAi
#
#total_freq = sum([val for key, val in AAi_freq.items()])
#print('total freq: ',total_freq)
#
#x = [raw_counts_lookup, AAi_freq]
#
#
######
#
#
#
#AAs = ['ALA', 'ARG', 'VAL', 'MET', 'GLU', 'GLN', 'GLY', \
#        'ASP', 'ASN', 'ILE', 'LEU', 'LYS', 'HIS', 'TRP', 'TYR', 'PHE', 'PRO', 'THR', 'SER', 'CYS']
#for aa in AAs:
#    sasa_dict[aa] = []
#
## fill sasa_dict with df values
#for ix, row in ifg_sasa.iterrows():
#    if row['resname_vdm'] in sasa_dict.keys():
#        resname = row['resname_vdm']
#        sasa_dict[resname].append(row['vdM_sasa_CB_3A_probe'])
#
## col = aa, row = sasa cutoffs
#index = [x*10 for x in range(20)][1:]
#
#all_sasa = [] # for every aa
#for aa, x in sasa_dict.items():
#    x = [float(z) for z in x ]
#    for i in x:
#        all_sasa.append(i)
#
#lookup = pkl.load(open('../Lookups/AAi_freq/AAi_database_lookups.pkl','rb')) # lookup table for freq(AAi) in db!
#
#df = pd.DataFrame(index = index, columns = AAs)
#for _bin in index:
#    bin_sasa = [z<_bin and z>(_bin-20) for z in all_sasa]
#    total_aa_bin = sum(bin_sasa) # counts of all aa's in this bin
#    for aa, x in sasa_dict.items():
#        x = [float(z) for z in x ]
#        x = [z<_bin and z>(_bin-20) for z in x]
#        aai_counts = sum(x) # at whatever cutoff bin
#        frac_aai_bin = aai_counts / total_aa_bin # f(AAi in binj / all aa in binj)
#        freq_aai_db = lookup[1][aa]
#        score = frac_aai_bin / freq_aai_db
#        df.ix[_bin, aa] = round(score,2)
#
#pkl.dump(df, open('../Lookups/sasa_scores/scores_largeprobe_sasa_%s_lookup.pkl' % ifg,'wb'))
#
