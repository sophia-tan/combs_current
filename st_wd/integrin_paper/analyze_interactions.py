#for each targetresi's int_res: log(sum all counts / sum all interactions)

import sys, os
from scipy import stats
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
from residues_integrin import *
import pickle as pkl, numpy as np, pandas as pd

freq = pkl.load(open('/home/gpu/Sophia/combs/st_wd/Lookups/AAi_freq/AAi_database_lookups.pkl','rb'))

def rmsd_filter(method,cutoff,score):
    print('--------------------------')
    print('method, cutoff, score =',method, cutoff, score)
    
    def triage(array,cutoff):
        new_ls = []
        for i in array:
            ix,rmsd = i
            if rmsd < cutoff:
                new_ls.append(True)
            else:
                new_ls.append(False)
        return new_ls

    def atomtriage(array,cutoff,num_atoms):
        new_ls = []
        for i in array:
            ix,rmsd = i
            newcut = num_atoms*cutoff
            if rmsd < newcut:
                new_ls.append(True)
            else:
                new_ls.append(False)
        return new_ls

    targetresis = []
    clus_for_all_targets = []
    avgclus_for_all_targets = []
    
    for targetres, resn in integrin_res.items():
    
        #print('----------------')
        #print('Target Res: ', targetres)
        clus_size = []
        avgclus = []
        for pklf in os.listdir('./output_data/'):
            if pklf.startswith(targetres[:4]) and pklf.endswith('rmsds_'+method+'.pkl'):
                ifg = pklf.split('_')[1][3:]
                vdm = pklf.split('_')[2][3:]
                lookup = pkl.load(open('./output_data/'+pklf,'rb')) # last element has atoms info
                num_atoms = lookup[-1][0]
                footnote_info = lookup[-1]
                integrin_clus = footnote_info[-2]
                avg_size_clus = footnote_info[-1]
                num_atoms = footnote_info[0]
                lookup = lookup[:-1]
                #triaged = atomtriage(lookup,cutoff=cutoff,num_atoms=num_atoms)
                triaged = triage(lookup,cutoff=cutoff)
                clus = integrin_clus
                clus_size.append(clus)
                avgclus.append(avg_size_clus)

                #print('Interacting residue: ', pklf.split('_')[2])
                #print(rawcount,integrin_clus,avg_size_clus)
        if clus_size == []:
            clus_size = [0.001]
            avgclus = [1]
        clus_size = np.array(clus_size)
        avgclus = np.array(avgclus)

        if score == 'a':
            #print('method: for each interacting res, take log(num/denom). then, sum.')
            clus_for_all_targets.append(sum(clus_size))
            avgclus_for_all_targets.append(sum((clus_size/avgclus)))
            #avgclus_for_all_targets.append(sum(np.log10(clus_size/avgclus)))
        
        if score == 'b':
            #print('method: for each interacting res, take log(num/denom). then, average.')
            clus_for_all_targets.append(np.mean(clus_size))
            avgclus_for_all_targets.append(np.mean((clus_size/avgclus)))
            #avgclus_for_all_targets.append(np.mean(np.log10(clus_size/avgclus)))

        if score == 'c':
            #print('method: sum across all num and denom for each interacting res. then, take the log')
            clus_for_all_targets.append(sum(clus_size))
            avgclus_for_all_targets.append((sum(clus_size)/sum(avgclus)))
            #avgclus_for_all_targets.append(np.log10(sum(clus_size)/sum(avgclus)))

        if score == 'd':
            #print('method: avg across all num and denom for each interacting res. then, take the log')
            clus_for_all_targets.append(np.mean(clus_size))
            avgclus_for_all_targets.append((np.mean(clus_size)/np.mean(avgclus)))
            #avgclus_for_all_targets.append(np.log10(np.mean(clus_size)/np.mean(avgclus)))

    return [clus_for_all_targets, avgclus_for_all_targets]

def get_correlation(megalist):
    AAs=['R671', 'I673', 'N753','F755','S758','V760', 'E785','H787','R900','L959']
    #     'D552', 'Y594', 'T603','H626','T656','K658', 'V664']
    activation = np.array([0.54,0.23,0.37,0.4,0.64,0.83,0.47,0.17,0.31,0.05])
    for calc in megalist:
        corr = stats.pearsonr(activation,calc)
        print(corr[0],corr[1])
        #print(calc,'calc')
        #corr = stats.spearmanr(activation,calc)
        #print(corr[0],corr[1])

for method in ['planar_group']:
    for score in ['a','b','c','d']:
        #for cut in [.5]:
        #for cut in [.005,.01,.03,.05]:
        for cut in [.3,.4,.5,.6]:
            megalist = rmsd_filter(method,cut,score)
            get_correlation(megalist)

#
    #print(sorted(set(targetresis)))
    #print(rawcountslist)

