# in command line, provide lookup table pkl file for freq(aai) in db
# also provide aa_source (db or carboxamdide)

import numpy as np
import pickle as pkl
from pprint import pprint
from sys import argv
import matplotlib.pyplot as plt
import pandas as pd

script, freq_aai_db, aa_source = argv

if aa_source == 'db':
    aa_dict=pkl.load(open('sasa_database_from_07252017.pkl','rb'))
elif aa_source == 'carboxylate': 
    aa_dict=pkl.load(open('carboxylate_largeprobesasa_rawcountsdict.pkl','rb'))

MaxASA = {'ALA': 129, 'ARG': 274, 'ASN': 195, 'ASP':193, 'CYS':167, 'GLU':223, 'GLN':225, 'GLY':104,
          'HIS':224, 'ILE':197, 'LEU':201, 'LYS':236, 'MET':224, 'PHE':240, 'PRO':159, 
          'SER':155, 'THR':172, 'TRP':285, 'TYR':263, 'VAL':174} 

# col = aa, row = sasa cutoffs
index = [x*10 for x in range(20)][1:]

all_sasa = [] # for every aa
for aa, resdict in aa_dict.items():
    x = resdict['3A'] # all sasa's for this aa w/ probe size 3A
    x = [float(z) for z in x ]
    for i in x:
        all_sasa.append(i)

lookup = pkl.load(open(freq_aai_db,'rb')) # lookup table for freq(AAi) in db!

df = pd.DataFrame(index = index, columns = MaxASA.keys())
for _bin in index:
    bin_sasa = [z<_bin and z>(_bin-20) for z in all_sasa]
    total_aa_bin = sum(bin_sasa) # counts of all aa's in this bin
    for aa, resdict in aa_dict.items():
        x = resdict['3A']
        x = [float(z) for z in x ]
        x = [z<_bin and z>(_bin-20) for z in x]
        aai_counts = sum(x) # at whatever cutoff bin
        frac_aai_bin = aai_counts / total_aa_bin # f(AAi in binj / all aa in binj)
        freq_aai_db = lookup[1][aa]
        if aa == 'ASP' or aa=='GLU':
            print(aa, aai_counts, total_aa_bin, frac_aai_bin, freq_aai_db)
            #print(aa, _bin, frac_aai_bin, freq_aai_db)
        score = frac_aai_bin / freq_aai_db
        df.ix[_bin, aa] = round(score,2)
        #print(aa, _bin, score)

#pkl.dump(df, open('scores_largeprobe_sasa_%s_lookup.pkl' % aa_source,'wb'))
