from pprint import pprint
from prody import *
import pickle as pkl
from sys import argv
import traceback
# generate lookup table for counts of AAi in nrPDB database
# need to provide textfile in command line!

# textfile is file that contains pdb and chains of all the pdbs in non-redundant PDB we're using
script, textfile = argv 

# compile lists of pdbs and chains
pdbnames = []
chains = []

with open(textfile) as fo:
    for line in fo:
        pdbnames.append(line[:4])
        chains.append(line[4])

zipped = zip(pdbnames, chains)

# create dictionary for lookup table
raw_counts_lookup = {} # key = AA, value = counts

for pdb, chain in zipped:
    try:
        print(pdb)
        parsed = parsePDB(pdb)
        parsed = parsed.select('chain %s and calpha' % chain)
        parsed = parsed.getResnames()
        for res in parsed:
            try:
                raw_counts_lookup[res] += 1
            except:
                raw_counts_lookup[res] = 1
    except Exception:
        print(traceback)
        print(pdb, 'pdb could not be parsed')

total_AAi = sum([val for key, val in raw_counts_lookup.items()])
print('total AA in database: ', total_AAi)

AAi_freq = {} # lookup for frequencies
for key, val in raw_counts_lookup.items():
    AAi_freq[key] = val / total_AAi

total_freq = sum([val for key, val in AAi_freq.items()])
print('total freq: ',total_freq)

x = [raw_counts_lookup, AAi_freq]
#pkl.dump(x, open('AAi_database_lookups.pkl', 'wb'))
