import pandas as pd
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
from scipy import stats

rel_sasa = pkl.load(open('rel_largeprobe_sasa_database_lookups.pkl','rb'))


Eisen = [0.620  ,-2.530  ,-0.780  ,-0.900  , 0.290  ,-0.850  ,-0.740  , 0.480  ,-0.400  , 1.380  , 1.060  ,-1.500  , 0.640  , 1.190  , 0.120  ,-0.180  ,-0.050  , 0.810  , 0.260  , 1.080 ]

AAs = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
df = pd.DataFrame(index=AAs, columns=['1.4A', '1.8A', '3A', '4A', 'Eisenberg'])

for index, aa in enumerate(AAs):
    df.ix[aa,'Eisenberg'] = Eisen[index]


for ix, row in rel_sasa[0].iterrows():
    aa = row.index
    for index, each in enumerate(['1.4A', '1.8A', '3A', '4A']):
        first_bin = rel_sasa[index].ix[10, aa]
        second_bin = rel_sasa[index].ix[20, aa]
        avgs = [np.mean([first_bin[ix], second_bin[ix]]) for ix in first_bin.index.values]
        df.ix[aa, each] = avgs

corr_14 = stats.pearsonr(df['1.4A'], df['Eisenberg']) 
corr_18 = stats.pearsonr(df['1.8A'], df['Eisenberg']) 
corr_3 = stats.pearsonr(df['3A'], df['Eisenberg']) 
corr_4 = stats.pearsonr(df['4A'], df['Eisenberg']) 

print('pearson coefficients for 1.4A, 1.8A, 3A, and 4A: ')
print(corr_14[0], corr_18[0], corr_3[0], corr_4[0])

plt.scatter(df['1.4A'],df['Eisenberg'])
plt.scatter(df['1.8A'],df['Eisenberg'])
plt.scatter(df['3A'],df['Eisenberg'])
plt.scatter(df['4A'],df['Eisenberg'])
plt.title('cb large probe sasa\'s vs. hydrophobicity scale; b=1.4A, o=1.8A, g=3A, r=4A')

plt.show()


