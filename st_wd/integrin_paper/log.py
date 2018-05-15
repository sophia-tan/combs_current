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
            if newcut > .6:
                newcut = .6
            if newcut < .4:
                print(newcut)
            if rmsd < newcut:
                new_ls.append(True)
            else:
                new_ls.append(False)
        return new_ls

    rawcounts_sum = []
    numclus_sum = []
    rawcounts_avg = []
    numclus_avg = []
    rawcount_norm_by_avgclus = []
    rawcount_norm_by_total_clus = []
    numclus_norm_by_avgclus = []
    numclus_norm_by_total_clus = []
    
    for targetres, resn in integrin_res.items():
    
        #print('----------------')
        #print('Target Res: ', targetres)
        clus_size = [] # num obs
        rawcounts = [] # num obs
        avgclus = [] # option to normalize by
        num_clustered = [] # option to normalize by
        for pklf in os.listdir('./output_data/'):
            if pklf.startswith(targetres[:4]) and pklf.endswith('rmsds_'+method+'.pkl'):
                ifg = pklf.split('_')[1][3:]
                vdm = pklf.split('_')[2][3:]
                lookup = pkl.load(open('./output_data/'+pklf,'rb')) # last element has atoms info
                footnote_info = lookup[-1]
                integrin_clus = footnote_info[-3]
                num_atoms = footnote_info[0]
                lookup = lookup[:-1]
                #triaged = atomtriage(lookup,cutoff=cutoff,num_atoms=num_atoms)
                triaged = triage(lookup,cutoff=cutoff)
                rawcount = sum(triaged)
                rawcounts.append(rawcount)
                num_clustered.append(footnote_info[-2])
                clus_size.append(integrin_clus)
                avgclus.append(footnote_info[-1])

                #print('Interacting residue: ', pklf.split('_')[2])
                #print(rawcount,integrin_clus,avg_size_clus)
        if rawcounts == []:
            rawcounts = [0.1]
            clus_size = [0.1]
            #rawcounts = [0]
            #clus_size = [0]
            avgclus = [1]
            num_clustered = [1]
        for r in rawcounts:
            if r==0:
                new = []
                for w in rawcounts:
                    if w!=0:
                        new.append(w)
                    else:
                        new.append(0.1)
                rawcounts=new

        clus_size = np.array(clus_size)
        rawcounts = np.array(rawcounts)
        avgclus = np.array(avgclus)
        num_clustered = np.array(num_clustered)

        rawcounts_sum.append(sum(rawcounts))
        numclus_sum.append(sum(clus_size))
        rawcounts_avg.append(np.mean(rawcounts))
        numclus_avg.append(np.mean(clus_size))
        
        if score == 'a':
            def calc(a,b): # for num 'a' and denom 'b'
                return sum(np.log10(a/b))
            
        if score == 'b':
            def calc(a,b): # for num 'a' and denom 'b'
                return np.log10(sum(a)/sum(b))
            
        if score == 'c':
            def calc(a,b): # for num 'a' and denom 'b'
                return np.mean(np.log10(a/b))

        if score == 'd':
            def calc(a,b): # for num 'a' and denom 'b'
                return np.log10(np.mean(a)/np.mean(b))

        rawcount_norm_by_avgclus.append(calc(rawcounts,avgclus))
        rawcount_norm_by_total_clus.append(calc(rawcounts,num_clustered))
        numclus_norm_by_avgclus.append(calc(clus_size,avgclus))
        numclus_norm_by_total_clus.append(calc(clus_size,num_clustered))

    return [rawcount_norm_by_avgclus, \
        rawcount_norm_by_total_clus, numclus_norm_by_avgclus, numclus_norm_by_total_clus] 
    #return [rawcounts_sum, numclus_sum, rawcounts_avg, numclus_avg, rawcount_norm_by_avgclus, \
    #    rawcount_norm_by_total_clus, numclus_norm_by_avgclus, numclus_norm_by_total_clus] 
    
def get_correlation(megalist):
    AAs=['R671', 'I673', 'N753','F755','S758','V760', 'E785','H787','R900','L959', \
         'D552', 'Y594', 'T603','H626','K658', 'V664']
    activation = np.array([0.54,0.23,0.37,0.4,0.64,0.83,0.47,0.17,0.31,0.05, \
        .12, .44, .399, .10, .20, .42])

    #activation = np.array([0.54,0.23,0.37,0.4,0.64,0.83,0.47,0.17,0.31,0.05])
    ''' ughhhhhh taking out the funky serine'''
    activation = [activation[x] for x in range(len(activation)) if x != 4] # delete 
    for origcalc in megalist:
        ''' ughhhhhh taking out the funky serine'''
        calc = [origcalc[x] for x in range(len(origcalc)) if x != 4] # delete

        spearcorr = stats.spearmanr(activation,calc)
        pearscorr = stats.pearsonr(activation,calc)
        #if spearcorr[0] > .5 and pearscorr[0] > .75:
        print(spearcorr[0],pearscorr[0],pearscorr[1])
        #print(origcalc,'calc')

for method in ['planar_group']:
    #for score in ['a']:
    for score in ['b']:
        for cut in [.5]:
            megalist = rmsd_filter(method,cut,score)
            get_correlation(megalist)

