import sys, os
from scipy import stats
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
from residues_integrin import *
import pickle as pkl, numpy as np, pandas as pd

db = pkl.load(open('/home/gpu/Sophia/combs/st_wd/Lookups/AAi_freq/AAi_database_lookups.pkl','rb'))

def rmsd_filter(method,cutoff,score):
    print('--------------------------')
    print('method, cutoff, score =',method, cutoff, score)
    
    norm_by_num_vdms = []
    norm_by_direct_vdms = []
    norm_by_avg_num_NNs = []
    norm_by_avg_num_NNs_nosing = []
    norm_by_med_num_NNs = []
    norm_by_med_num_NNs_nosing = []

    for targetres, resn in integrin_res.items():
    
        print('----------------')
        print('Target Res: ', targetres)
        rawcounts = [] # num obs
        num_vdms = [] # option to normalize by
        num_clustered = [] # option to normalize by
        avgclus = [] # option to normalize by
        avgclus_wo_sing = [] # option to normalize by
        medclus = [] # option to normalize by
        medclus_wo_sing = [] # option to normalize by
        for pklf in os.listdir('./output_data/'):
            if pklf.startswith(targetres[:4]) and pklf.endswith(
                'matches_{}_{}.pkl'.format(method,cutoff)):
                matches = pkl.load(open('./output_data/'+pklf,'rb')) 
                [rmsd, resi, ifgresn, ifgresi, vdmresn, vdmresi,
                        num_nn, norm_metrics] = matches
                num_all_vdms, num_direct, [avgsize, avgsize_no_sing,
                    medsize, medsize_no_sing] = norm_metrics

                print(num_all_vdms, num_direct)
                rawcounts.append(num_nn+2) # it's actually effectively +1, but I'm stupid 
                # and returned num_nn - 1 in get_NN() when I shouldn't have
                num_vdms.append(num_all_vdms + 1)
                print(num_direct, db[1][vdmresn], db[1][ifgresn])
                num_direct = num_direct*db[1][vdmresn]* db[1][ifgresn]

                num_clustered.append(num_direct + 1) 
                avgclus.append(avgsize + 1)
                avgclus_wo_sing.append(avgsize_no_sing + 1)
                medclus.append(medsize + 1)
                medclus_wo_sing.append(medsize_no_sing + 1)

                print('Interacting residue: ', pklf.split('_')[2])
                print(num_nn, num_direct)
        if rawcounts == []:
            rawcounts = [1]
            num_vdms = [10000]
            num_clustered = [10000]
            for ls in [avgclus, avgclus_wo_sing, medclus, medclus_wo_sing]:
                ls.append(10) 

        ######################################################################
        # Diff methods of calculating score for list of observed counts 'a' 
        # and list of expected counts 'b'
        ######################################################################
        def list_to_ndarray(a, b):
            return [np.array(a), np.array(b)]
        
        if score == 'a':
            def calc(a,b): 
                a, b = list_to_ndarray(a, b)
                return np.log10(sum(a/b))
            
        if score == 'b':
            def calc(a,b): 
                a, b = list_to_ndarray(a, b)
                return sum(np.log10(a/b))
            
        if score == 'c':
            def calc(a,b): 
                a, b = list_to_ndarray(a, b)
                return np.log10(sum(a)/sum(b))
            
        if score == 'd':
            def calc(a,b): 
                a, b = list_to_ndarray(a, b)
                return np.log10(np.mean(a/b))

        if score == 'e':
            def calc(a,b): 
                a, b = list_to_ndarray(a, b)
                return np.mean(np.log10(a/b))

        if score == 'f':
            def calc(a,b): 
                a, b = list_to_ndarray(a, b)
                return np.log10(np.mean(a)/np.mean(b))

        #print(rawcounts)
        #print(avgclus)
        
        norm_by_direct_vdms.append(calc(rawcounts, num_clustered))

    #return [norm_by_num_vdms, norm_by_direct_vdms, norm_by_avg_num_NNs,
    #    norm_by_avg_num_NNs_nosing, norm_by_med_num_NNs, norm_by_med_num_NNs_nosing]
    return [norm_by_direct_vdms] 

def get_correlation(megalist):
    AAs=['R671', 'I673', 'N753','F755','S758','V760', 'E785','H787','R900','L959', \
         'D552', 'Y594', 'T603','H626','K658', 'V664', 'E534']
    activation = np.array([0.54,0.23,0.37,0.4,0.64,0.83,0.47,0.17,0.31,0.05, \
        .12, .44, .399, .10, .20, .42, 0.76])
    #AAs=['R671', 'I673', 'N753','F755','S758','V760', 'E785','H787','R900','L959', \
    #     'D552', 'Y594', 'T603','H626','K658', 'V664']
    #activation = np.array([0.54,0.23,0.37,0.4,0.64,0.83,0.47,0.17,0.31,0.05, \
    #    .12, .44, .399, .10, .20, .42])

    ''' ughhhhhh taking out the funky serine'''
    activation = [activation[x] for x in range(len(activation)) if x != 4] # delete 
    for origcalc in megalist:
        ''' ughhhhhh taking out the funky serine'''
        calc = [origcalc[x] for x in range(len(origcalc)) if x != 4] # delete
        spearcorr = stats.spearmanr(activation,calc)
        pearscorr = stats.pearsonr(activation,calc)
        print('corr', pearscorr)
        print(origcalc,'calc')

for method in ['planar_group_no_bb']:
    for score in ['a','b','c','d','e','f']:
        for cut in [.5]:
        #for cut in [.3,.4,.5]:
            megalist = rmsd_filter(method,cut,score)
            get_correlation(megalist)
