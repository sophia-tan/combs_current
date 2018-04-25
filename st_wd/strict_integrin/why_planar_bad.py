#for each targetresi's int_res: log(sum all counts / sum all interactions)

import sys, os
from scipy import stats
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
import pickle as pkl, numpy as np, pandas as pd


freq = pkl.load(open('/home/gpu/Sophia/combs/st_wd/Lookups/AAi_freq/AAi_database_lookups.pkl','rb'))

def rmsd_filter(method,cutoff,score):
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
            targetresis.append(pklf.split('_')[0]+pklf.split('_')[1][3:])
    
    rawcounts_for_all_targets = []
    num_vdms_method_for_all_targets = []
    freq_method_for_all_targets = []
    num_interacting_vdms_method_for_all_targets = []
    
    for targetres in sorted(set(targetresis)):
        print('----------------')
        print('Target Res: ', targetres)
        rawcounts = []
        num_vdms = []
        exp_based_on_freq = []
        num_interacting_vdms = []
        for pklf in os.listdir('./output_data/'):
            if pklf.startswith(targetres[:4]) and pklf.endswith('_'+method+'.pkl'):
                ifg = targetres[4:]
                vdm = pklf.split('_')[2][3:]
                lookup = pkl.load(open('./output_data/'+pklf,'rb')) # index 0 has atoms info
                num_atoms = lookup[0][0]
                num_vdms_interacting = lookup[0][-1]
                lookup = lookup[1:]
                #triaged = atomtriage(lookup,cutoff=cutoff,num_atoms=num_atoms)
                triaged = triage(lookup,cutoff=cutoff)
                rawcount = sum(triaged)
                if rawcount==0:
                    rawcount=0.001
                rawcounts.append(rawcount)
                num_vdms.append(len(triaged))
                exp_based_on_freq.append(freq[0][ifg]*freq[1][vdm])
                if num_vdms_interacting == 0:
                    num_vdms_interacting = 1
                num_interacting_vdms.append(num_vdms_interacting)

                print('Interacting residue: ', pklf.split('_')[2])
                print(rawcount)
                #print(lookup.iloc[0]['rmsds'][1])
        rawcounts = np.array(rawcounts)
        num_vdms = np.array(num_vdms)
        exp_based_on_freq = np.array(exp_based_on_freq)
        num_interacting_vdms = np.array(num_interacting_vdms)
        if score == 'a':
            #print('method: for each interacting res, take log(num/denom). then, sum.')
            #print('score method "a"')
            rawcounts_for_all_targets.append(sum(rawcounts))
            num_vdms_method_for_all_targets.append(sum(np.log10(rawcounts/num_vdms)))
            freq_method_for_all_targets.append(sum(np.log10(rawcounts/exp_based_on_freq)))
            num_interacting_vdms_method_for_all_targets.append(sum(np.log10(rawcounts/num_interacting_vdms)))
        
        if score == 'b':
            #print('method: for each interacting res, take log(num/denom). then, average.')
            #print('score method "b"')
            rawcounts_for_all_targets.append(np.mean(rawcounts))
            num_vdms_method_for_all_targets.append(np.mean(np.log10(rawcounts/num_vdms)))
            freq_method_for_all_targets.append(np.mean(np.log10(rawcounts/exp_based_on_freq)))
            num_interacting_vdms_method_for_all_targets.append(np.mean(np.log10(rawcounts/num_interacting_vdms)))

        if score == 'c':
            #print('method: sum across all num and denom for each interacting res. then, take the log')
            #print('score method "c"')
            rawcounts_for_all_targets.append(sum(rawcounts))
            num_vdms_method_for_all_targets.append(np.log10(sum(rawcounts)/sum(num_vdms)))
            freq_method_for_all_targets.append(np.log10(sum(rawcounts)/sum(exp_based_on_freq)))
            num_interacting_vdms_method_for_all_targets.append(np.log10(sum(rawcounts)/sum(num_interacting_vdms)))

        if score == 'd':
            #print('method: avg across all num and denom for each interacting res. then, take the log')
            #print('score method "d"')
            rawcounts_for_all_targets.append(np.mean(rawcounts))
            num_vdms_method_for_all_targets.append(np.log10(np.mean(rawcounts)/np.mean(num_vdms)))
            freq_method_for_all_targets.append(np.log10(np.mean(rawcounts)/np.mean(exp_based_on_freq)))
            num_interacting_vdms_method_for_all_targets.append(np.log10(np.mean(rawcounts)/np.mean(num_interacting_vdms)))
    
    print([round(x) for x in rawcounts_for_all_targets])
    megalist = [rawcounts_for_all_targets, num_vdms_method_for_all_targets,freq_method_for_all_targets,\
                num_interacting_vdms_method_for_all_targets]
    return megalist

def get_correlation(megalist):
    AAs=['R671', 'I673', 'N753','F755','S758','V760', 'E785','H787','R900','L959']
    print(AAs)
    activation = np.array([0.54,0.23,0.37,0.4,0.64,0.83,0.47,0.17,0.31,0.05])
    for calc in megalist:
        corr = stats.pearsonr(activation,calc)
        print(corr[0],corr[1])
        #print(calc,'calc')
        #corr = stats.spearmanr(activation,calc)
        #print(corr[0],corr[1])


for method in ['planar_group']:
    for cutoff in [.4]:
    #for cutoff in [.5]:
        for score in ['a']:
            print('+++++++++++++++++++++++++++++++')
            print('score method',score)
            megalist = rmsd_filter(method,cutoff,score)
            get_correlation(megalist)



#
    #print(sorted(set(targetresis)))
    #print(rawcountslist)

