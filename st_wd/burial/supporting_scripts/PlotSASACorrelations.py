# commandline argument is fs, dssp, 3A, 4A, or 5A. need 2
import numpy as np
import pickle as pkl
from pprint import pprint
from sys import argv
import matplotlib.pyplot as plt

script, m, n = argv

nrpdb_dict=pkl.load(open('nrPDB_sasa.pkl','rb'))

megadict = {'3A': [], '4A': [], '5A': [], 'dssp':[], 'fs':[]} # for plotting

for key, val in nrpdb_dict.items():
    for k, v in val.items():
        megadict[k] = megadict[k] + v


x = megadict[m]
print(len(x))
x = x[:9000]
x = [float(z) for z in x ]#if (np.isnan(z) == False)]
#x = [z for z in x if not np.isnan(z)]
#x = [z for z in x if z < 200 and z > 1] ### WTF???
#print(max(x), min(x))

y = megadict[n]
print(len(y))
y = y[:9000]
y = [float(z) for z in y ]#if (np.isnan(z) == False)]
#y = [z for z in y if not np.isnan(z)]



# count how many nans at some point
plt.scatter(x,y,s=0.5)
plt.xlabel(m)
plt.ylabel(n)
plt.show()
#AAs = ['ALA', 'ARG', 'VAL', 'MET', 'GLU', 'GLN', 'GLY', \
#        'ASP', 'ASN', 'ILE', 'LEU', 'LYS', 'HIS', 'TRP', 'TYR', 'PHE', 'PRO', 'THR', 'SER', 'CYS']
