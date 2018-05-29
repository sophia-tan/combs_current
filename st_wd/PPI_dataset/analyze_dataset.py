''' get statistics on the dataset '''
import matplotlib.pyplot as plt
import pandas as pd

zemu = pd.read_csv('kortemme_flexddg_dataset.csv')
print(len(set(zemu['PDBFileID'])), 'unique PDBs')

ddgs = zemu['Experimental ddG']
plt.hist(ddgs,bins=16)
plt.xlabel('Experimental ddG')
plt.ylabel('counts')
plt.show()
