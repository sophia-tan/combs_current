# commandline argument is fs, dssp, 3A, 4A, or 5A
import numpy as np
import pickle as pkl
from pprint import pprint
from sys import argv
import matplotlib.pyplot as plt
import pandas as pd

script, m = argv

nrpdb_dict=pkl.load(open('nrPDB_sasa.pkl','rb'))

MaxASA = {'ALA': 129, 'ARG': 274, 'ASN': 195, 'ASP':193, 'CYS':167, 'GLU':223, 'GLN':225, 'GLY':104,
          'HIS':224, 'ILE':197, 'LEU':201, 'LYS':236, 'MET':224, 'PHE':240, 'PRO':159, 
          'SER':155, 'THR':172, 'TRP':285, 'TYR':263, 'VAL':174} 

# col = aa, row = RSA cutoffs
index = [x/20 for x in range(20)][1:]
df = pd.DataFrame(index = index, columns = MaxASA.keys())

for perc in index:
    
    all_sasa = []
    for aa, resdict in nrpdb_dict.items():
    
        x = resdict[m]
        x = [float(z) for z in x ]
        #n = [z for z in x if np.isnan(z)]
        #print(aa, len(n))
        #x = [z for z in x if not np.isnan(z)]
        x = np.hstack(x)
        x = x/MaxASA[aa]
        for i in x:
            all_sasa.append(i)
    
    all_sasa = sorted(all_sasa)
    all_sasa = [z<perc and z>(perc-0.05) for z in all_sasa]
    total_counts = sum(all_sasa) # for whatever % cutoff for all aa sasa
    for aa, resdict in nrpdb_dict.items():
        x = resdict[m]
        total_aa = len(x)
        x = [float(z) for z in x ]
        x = np.hstack(x)
        x = x/MaxASA[aa]
        x = sorted(x)
        x = [z<perc and z>(perc-0.05) for z in x]
        aa_counts = sum(x) # at whatever % cutoff
        frac_aa_bin = aa_counts / total_aa # f(aa in argv bin for that aa)
        frac_allsasa_bin = total_counts / len(all_sasa) # f(aas in argv bin for all aas)
        score = frac_aa_bin / frac_allsasa_bin
        df.ix[perc, aa] = round(score,2)

pkl.dump(df, open('rsa_dssp_database_lookups.pkl', 'wb'))

#    f, axarr = plt.subplots(4,5,figsize=(12,6.5))
#    fig_row = 0
#    fig_col = 0
#    
#        axarr[fig_row, fig_col].hist(x, bins=10)
#        axarr[fig_row, fig_col].set_title(aa)
#    
#        if fig_col < 4:
#            fig_col += 1
#        else:
#            fig_col = 0
#            fig_row += 1
#
#plt.tight_layout()
#plt.suptitle('relative sasa')
#plt.subplots_adjust(top=0.85) # so suptitle is compatible with tight layout
#plt.show()
####
