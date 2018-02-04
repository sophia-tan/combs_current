import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis
from sys import argv
import numpy as np, pandas as pd, pickle as pkl
from numpy import concatenate, array_equal, dot, cross, linalg, abs
from scipy.spatial import ConvexHull, Delaunay


database = pkl.load(open('../database_mindistconvexhull.pkl','rb'))


AAs = ['ALA', 'ARG', 'VAL', 'MET', 'GLU', 'GLN', 'GLY', \
        'ASP', 'ASN', 'ILE', 'LEU', 'LYS', 'HIS', 'TRP', 'TYR', 'PHE', 'PRO', 'THR', 'SER', 'CYS']
db_aai_freq = pkl.load(open('../../Lookups/AAi_freq/AAi_database_lookups.pkl','rb')) # lookup table for freq(AAi) in db!
# col = aa, row = burial dist cutoffs
index = [x for x in range(11)][1:]
db_df = pd.DataFrame(index = index, columns = AAs)

df = database
for _bin in index: #_bin is the upperbound (doesn't include)
    total_bin = sum([z<_bin and z>=(_bin-1) for z in database['burialdist']])
    for aa in AAs:
        aai = df[df['residuename']==aa]['burialdist']
        aai_bin_counts = sum([z<_bin and z>=(_bin-1) for z in aai])
        frac_aai_bin = aai_bin_counts / total_bin # f(AAi in binj / all aa in binj)
        freq_aai_db = db_aai_freq[1][aa]
        score = -np.log10(frac_aai_bin / freq_aai_db)
        db_df.ix[_bin, aa] = round(score,2)

pkl.dump(db_df, open('/home/gpu/Sophia/combs/st_wd/Lookups/burial_scores/db_burial_scores.pkl','wb'))
