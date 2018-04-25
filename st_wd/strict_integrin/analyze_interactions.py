import sys, os
from scipy import stats
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
import pickle as pkl, numpy as np, pandas as pd


freq = pkl.load(open('/home/gpu/Sophia/combs/st_wd/Lookups/AAi_freq/AAi_database_lookups.pkl','rb'))

def rmsd_filter(method,cutoff):
    print('--------------------------')
    print('method, cutoff =',method, cutoff)
    
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
            if rmsd/num_atoms < cutoff:
                new_ls.append(True)
            else:
                new_ls.append(False)
        return new_ls
    


    targetresis = []
    for pklf in sorted(os.listdir('./output_data/')):
        if pklf.endswith('_'+method+'.pkl'):
            #x=pkl.load(open('./output_data/'+pklf,'rb'))
            targetresis.append(pklf.split('_')[0]+pklf.split('_')[1][3:])
    

    rawcountslist = []
    num_ifg_vdm_nums_list = []
    freq_method_list = []
    interacting_method = []
    
    for targetres in sorted(set(targetresis)):
        #print('----------------')
        #print('Target Res: ', targetres)
        rawcounts = 0
        all_poss_intrxn= 0
        freq_method = 0
        interacting = 0
        for pklf in os.listdir('./output_data/'):
            if pklf.startswith(targetres[:4]) and pklf.endswith('_'+method+'.pkl'):
                ifg = targetres[4:]
                vdm = pklf.split('_')[2][3:]
                lookup = pkl.load(open('./output_data/'+pklf,'rb')) # index 0 has atoms info
                num_atoms = lookup[0][0]
                num_interacting_vdms = lookup[0][-1]
                lookup = lookup[1:]
                #triaged = atomtriage(lookup,cutoff=cutoff,num_atoms=num_atoms)
                triaged = triage(lookup,cutoff=cutoff)
                rawcounts += sum(triaged)
                all_poss_intrxn += len(triaged) 
                freq_method += freq[0][ifg]*freq[1][vdm] 
                #interacting += 1
                interacting += num_interacting_vdms

                    #print(dist, interacting)
                #print('Interacting residue: ', pklf.split('_')[2])
                #print(lookup.iloc[0]['rmsds'][1])
                ##''''''''print('rawcounts,', '# of vdMs,','# of iFGs * freq(vdm AA)')
                ##''''''''print(rawcounts, all_poss_intrxn, freq_method,interacting)
                #print(all_poss_intrxn)
                #print(freq_method)
                ##print('nan',nan_count)
        rawcountslist.append(rawcounts)
        num_ifg_vdm_nums_list.append(rawcounts/all_poss_intrxn)
        freq_method_list.append(rawcounts/freq_method)
        interacting_method.append(rawcounts/interacting)
        #print(rawcounts,freq_method)
    
    
    AAs=['R671', 'I673', 'N753','F755','S758','V760', 'E785','H787','R900','L959']
    activation = np.array([0.54,0.23,0.37,0.4,0.64,0.83,0.47,0.17,0.31,0.05])
    #calc = [1.53,0.12,1.32,1.76,0.71,1.20,0.54]
    #for calc in [freq_method_list,interacting_method]:
    for calc in [rawcountslist,num_ifg_vdm_nums_list,freq_method_list,interacting_method]:
        corr = stats.pearsonr(activation,calc)
        print(corr[0],corr[1], 'correlation')
        #print(calc,'calc')
        #corr = stats.spearmanr(activation,calc)
        #print(corr[0],corr[1])

for method in ['planar_group']:
#for method in ['whole_res', 'BBorSC', 'FGorBBorSC', 'direct_intrxn']:
    for cutoff in [.2,.3,.5]:
    #for cutoff in [.005,.01,.05,.1]:
        rmsd_filter(method,cutoff)



#
    #print(sorted(set(targetresis)))
    #print(rawcountslist)

