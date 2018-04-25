import sys, os
from scipy import stats
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
import pickle as pkl, numpy as np


fgintrxns = [35160,6336,10710,24775,18809,7348,8275]

for cutoff in [.5]:
    print('cutoff=:',cutoff)
    
    def triage(row,cutoff):
        rmsd,num_atoms = row['rmsds']
        if rmsd < cutoff:
            return True
        else:
            return False
    
    freq = pkl.load(open('/home/gpu/Sophia/combs/st_wd/Lookups/AAi_freq/AAi_database_lookups.pkl','rb'))
    
    targetresis = []
    for pklf in sorted(os.listdir('./output_data/')):
        if pklf.endswith('.pkl'):
            targetresis.append(pklf.split('_')[0]+pklf.split('_')[4])
    
    rawcountslist = []
    num_ifg_vdm_nums_list = []
    fg_intrxns_list = []

    for index,targetres in enumerate(sorted(set(targetresis))):
        #print('----------------')
        #print('Target Res: ', targetres)
        rawcounts = 0
        all_poss_intrxn= 0
        freq_method = 0
        for pklf in os.listdir('./output_data/'):
            if pklf.startswith(targetres[:4]):
                ifg = constants.AAname_rev[targetres[4:]]
                vdm = constants.AAname_rev[pklf.split('_')[6]]
                lookup = pkl.load(open('./output_data/'+pklf,'rb'))
                no_nan = lookup.dropna(axis=0)
                nan_count = len(lookup)-len(no_nan)
                triaged = no_nan.apply(triage,axis=1,cutoff=cutoff)
                rawcounts += sum(triaged)
                all_poss_intrxn += len(triaged)
                #print('Interacting residue: ', pklf.split('_')[6])
                #print('nan',nan_count)
        rawcountslist.append(rawcounts)
        num_ifg_vdm_nums_list.append(rawcounts/all_poss_intrxn)
        fg_intrxns_list.append(rawcounts/fgintrxns[index])
    print(rawcountslist)
    
    AAs=['R671', 'I673', 'N753','F755','V760', 'E785','R900']
    activation = np.array([0.61,0.23,0.54,0.42,0.83,0.47,0.31])
    
    # rosetta
    calc = [1.53,0.12,1.32,1.76,0.71,1.20,0.54]
    
    for calc in [rawcountslist,num_ifg_vdm_nums_list,fg_intrxns_list]:
        corr = stats.pearsonr(activation,calc)
        print(corr[0])
        corr = stats.spearmanr(activation,calc)
        print(corr[0])
#

#
#print(sorted(set(targetresis)))
#print(rawcountslist)
#print(num_ifg_vdm_nums_list)
#print(freq_method_list)
#
