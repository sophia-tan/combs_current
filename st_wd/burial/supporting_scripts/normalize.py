import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis
from sys import argv
import numpy as np, pandas as pd, pickle as pkl
from numpy import concatenate, array_equal, dot, cross, linalg, abs
from scipy.spatial import ConvexHull, Delaunay

script, ifg = argv

database = pkl.load(open('../database_mindistconvexhull.pkl','rb'))
# get the max depth of each pdb
depths = {}
grouped  = database.groupby('pdbcode')
for pdb, dframe in grouped:
    depths[pdb] = max(dframe['burialdist'])
# only return residue if it's within 10% most buried 
# relative to its protein
def norm(row):
    if row['burialdist'] > 0.8*depths[row['pdbcode']]:
        return True
    else:
        return False

df=database[database.apply(norm, axis=1)]


AAs = ['ALA', 'ARG', 'VAL', 'MET', 'GLU', 'GLN', 'GLY', \
        'ASP', 'ASN', 'ILE', 'LEU', 'LYS', 'HIS', 'TRP', 'TYR', 'PHE', 'PRO', 'THR', 'SER', 'CYS']
db_aai_freq = pkl.load(open('../../Lookups/AAi_freq/AAi_database_lookups.pkl','rb')) # lookup table for freq(AAi) in db!
# col = aa, row = burial dist cutoffs
index = [0]
db_df = pd.DataFrame(index = index, columns = AAs)

total = len(df)
for aa in AAs:
    aai = df[df['residuename']==aa]['burialdist']
    aai_counts = len(aai)
    frac_aai = aai_counts/total
    freq_aai_db = db_aai_freq[1][aa]
    score = frac_aai / freq_aai_db
    db_df.ix[0, aa] = round(score,2)

pkl.dump(db_df, open('/home/gpu/Sophia/combs/st_wd/Lookups/burial_scores/percentage_%s_burial_scores.pkl'%ifg,'wb'))
