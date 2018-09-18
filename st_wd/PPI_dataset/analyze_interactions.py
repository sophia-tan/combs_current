import sys, os
from scipy import stats
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
import pickle as pkl, numpy as np, pandas as pd
import matplotlib.pyplot as plt
sys.path.append('/home/gpu/Sophia/combs/st_wd/Scoring_Function')
from CalculateScores import *

##### is this for PPI_dataset or integrin? ####
study = 'PPI_dataset'

def add_one(row):
    return row+1

### ANALYZING DETAILS ###
def rmsd_filter(method, score_method, nonmembrane, buried, rmsd):
    print('--------------------------')
    print('method', method)
    print('score method', score_method)
    print(buried, nonmembrane, rmsd)
    
    norm_by_num_vdms = []
    norm_by_direct_vdms = []
    norm_by_avg_num_NNs = []
    norm_by_avg_num_NNs_nosing = []
    norm_by_med_num_NNs = []
    norm_by_med_num_NNs_nosing = []
    colors = []
    labels = []
    ddgs_of_no_intrxn = []
    actual_ddgs = []

    f, axarr = plt.subplots(3,2)
    for pdb_index in range(54):
        pklf = 'pdb_ix_{}_matches_{}.pkl'.format(pdb_index, method)
        # make variable 'matches' that contains all pdb_index pklfs
        if pdb_index == 0:
            matches = pkl.load(open('./output_data/'+pklf,'rb'))
        else:
            new = pkl.load(open('./output_data/'+pklf,'rb'))
            matches = pd.concat([matches,new])
        
        mutation_data = pkl.load(open('mutation_data.pkl','rb'))
        pdb_datadf = mutation_data[1][int(pdb_index)]

        filtered_matches = matches[matches['rmsd']==rmsd]
        filtered_matches = filtered_matches[filtered_matches['nonmembrane']==nonmembrane]
        filtered_matches = filtered_matches[filtered_matches['buried']==buried]

        '''For each experimental sample, make computational prediction'''
        
        for row_ix, row in pdb_datadf.iterrows():
            ddg = row['ddG']
            combo_mutation = row['Mutation']
            predictions = filtered_matches[filtered_matches[
                'mut string'].astype(str) == str(combo_mutation)]
            
            ## ignore where num NN = 0
            #predictions = predictions[predictions['num NN'] > 0]
            
            num_nn = predictions['num NN'].apply(add_one)
            num_vdms = predictions['num vdms'].apply(add_one)
            direct_vdms = predictions[
                'num directly interacting vdms'].apply(add_one)
            avg_NNs = predictions['avg num NNs'].apply(add_one)
            avg_NNs_nosing = predictions[
                'avg num NNs w/o singles'].apply(add_one)
            med_NNs = predictions['median num NNs'].apply(add_one)
            med_NNs_nosing = predictions[
                'median num NNs w/o singles'].apply(add_one)
                
            if len(predictions) == 0:
                ddgs_of_no_intrxn.append(ddg)

            elif len(predictions) > 0:
                actual_ddgs.append(ddg)
                norm_by_num_vdms.append(calc_score(num_nn, num_vdms, score_method))
                norm_by_direct_vdms.append(calc_score(num_nn, direct_vdms, score_method))
                norm_by_avg_num_NNs.append(calc_score(num_nn, avg_NNs, score_method))
                norm_by_avg_num_NNs_nosing.append(calc_score(num_nn, avg_NNs_nosing, score_method))
                norm_by_med_num_NNs.append(calc_score(num_nn, med_NNs, score_method))
                norm_by_med_num_NNs_nosing.append(calc_score(num_nn, med_NNs_nosing, score_method))
                wt = [x[1][0] for x in combo_mutation]
                labels.append(wt)
            
            if len(predictions)==1:
                colors.append('b')
            elif len(predictions)==2:
                colors.append('green')
            else:
                colors.append('black')
        
    r,c=0, 0
    for pred, typ in zip([norm_by_num_vdms, norm_by_direct_vdms,
        norm_by_avg_num_NNs, norm_by_avg_num_NNs_nosing,
        norm_by_med_num_NNs, norm_by_med_num_NNs_nosing], 
        ['num vdms', 'direct vdms', 'avg # NN', 'avg # NN no single',
        'med # NN', 'med # NN no single']):
        pearson_r = stats.pearsonr(actual_ddgs, pred)
        print(pearson_r, 'pearson')
        spearman_r = stats.spearmanr(actual_ddgs, pred)
        axarr[r,c].scatter(pred,actual_ddgs,marker='o',
                lw=1,color=colors,s=10)
        #no_intrxn = [0 for i in ddgs_of_no_intrxn]
        #axarr[r,c].scatter(no_intrxn,ddgs_of_no_intrxn,marker='o',
        #        lw=1,color='red',s=10)
        axarr[r,c].set_title(typ)

        #for label, x, y in zip(labels, pred, all_ddgs):
        #    text(label,x,y,axarr[r,c])
    
        if c==1:
            c=0
            r+=1
        else:
            c+=1
    plt.suptitle('{}_{}_{}_{}'.format(rmsd, score_method, nonmembrane, buried))
    #plt.show()
    plt.savefig('./output_data/{}_{}_{}_{}_{}.png'.format(rmsd, score_method, nonmembrane, buried, method))


for method in ['sc_only_ifg']:
#for method in ['planar_group_no_bb']:
    for score_method in ['d', 'e','f','g','h','i','j','k','l','m','n','o','p','q']:
        for rmsd in [.2,.3,.4,.5,.6,.7,.8,.9]:
            for nonmembrane in [True, False]:
                for buried in [True, False]:
                    rmsd_filter(method, score_method, nonmembrane, buried, rmsd)
