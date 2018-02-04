import pandas as pd
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
from scipy import stats
from sys import argv

script, no_arg = argv

_bin = 0

db_df = pkl.load(open('/home/gpu/Sophia/combs/st_wd/Lookups/burial_scores/percentage_db_burial_scores.pkl','rb'))

Eisen = [0.620  ,-2.530  ,-0.780  ,-0.900  , 0.290  ,-0.850  ,-0.740  , 0.480  ,-0.400  , 1.380  , 1.060  ,-1.500  , 0.640  , 1.190  , 0.120  ,-0.180  ,-0.050  , 0.810  , 0.260  , 1.080 ]
AAs = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
df = pd.DataFrame(index=AAs, columns=['db', 'Eisenberg'])

for index, aa in enumerate(AAs):
    df.ix[aa,'Eisenberg'] = Eisen[index]
x = []
y = []
list_df = [db_df]
for aa in AAs:
    if no_arg == 'False':
        for index, each in enumerate(['db']):
            first_bin = list_df[index].ix[_bin, aa]
            df.ix[aa,each] = first_bin
    if no_arg == 'True':
        if aa != 'ARG':
            for index, each in enumerate(['db']):
                first_bin = list_df[index].ix[_bin, aa]
                df.ix[aa,each] = first_bin

df = df.dropna(axis=0, how='any')
corr_db = stats.pearsonr(df['db'], df['Eisenberg']) 

print(corr_db[0])

if no_arg=='True':
    AAs = ['ALA','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']

for label, x, y in zip(AAs, df['db'], df['Eisenberg']):
    plt.annotate(
        label,
        xy=(x, y), xytext=(10, 10),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.1', fc='yellow', alpha=0.2),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

plt.ylabel('Eisenberg score')
plt.xlabel('min dist to convex hull score')

plt.scatter(df['db'],df['Eisenberg'])
plt.title('10% depth vs. hydrophobicity scale')

plt.show()
