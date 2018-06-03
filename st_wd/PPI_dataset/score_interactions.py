import sys, os
import matplotlib.pyplot as plt
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
import prody as pr
import numpy as np, pickle as pkl, pandas as pd
from itertools import *
from Scoring import *
from Functions import *
from scipy import stats

def add_one(row):
    return row+1

for rmsd in [.5]:
    #for calc_method in ['a']:
    for calc_method in ['a','b','c','d']:
        f, axarr = plt.subplots(3,2)
        all_ddgs = []
        all_rawscores = []
        all_sums = []
        all_means = []
        labels = []
        norm_by_num_vdms = []
        norm_by_direct_vdms = []
        norm_by_avg_num_NNs = []
        norm_by_avg_num_NNs_nosing = []
        norm_by_med_num_NNs = []
        norm_by_med_num_NNs_nosing = []

        for pdb_index in range(54):
            pklf='pdb_ix_{}_matches.pkl'.format(pdb_index)
            if pklf in os.listdir('./output_data'):
                if pdb_index==0:
                    matches = pkl.load(open('./output_data/'+pklf,'rb'))
                else:
                    new = pkl.load(open('./output_data/'+pklf,'rb'))
                    matches = pd.concat([matches,new])
        
                mutation_data = pkl.load(open('mutation_data.pkl','rb'))
                pdb_datadf = mutation_data[1][int(pdb_index)]

                filtered_matches = matches[matches['rmsd']==rmsd]
                '''For each experimental sample, make computational prediction'''
                
                for row_ix, row in pdb_datadf.iterrows():
                    ddg = row['ddG']
                    combo_mutation = str(row['Mutation'])
                    predictions = filtered_matches[filtered_matches[
                        'mut string'].astype(str) == combo_mutation]
                    # later, do the opposite to get it where it's 0
                    predictions = predictions[predictions['num NN'] > 0]
                    labels.append(combo_mutation)
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

                    if len(predictions) > 0:
                        all_ddgs.append(ddg)
                        if calc_method == 'a':
                            def calc(a,b): # for num 'a' and denom 'b'
                                return np.log10(sum(a/b))
                            
                        if calc_method == 'b':
                            def calc(a,b): # for num 'a' and denom 'b'
                                return sum(np.log10(a/b))
                            
                        if calc_method == 'c':
                            def calc(a,b): # for num 'a' and denom 'b'
                                return np.log10(np.mean(a/b))
                
                        if calc_method == 'd':
                            def calc(a,b): # for num 'a' and denom 'b'
                                return np.mean(np.log10(a/b))
                
                        norm_by_num_vdms.append(calc(num_nn, num_vdms))
                        norm_by_direct_vdms.append(calc(num_nn, direct_vdms))
                        norm_by_avg_num_NNs.append(calc(num_nn, avg_NNs))
                        norm_by_avg_num_NNs_nosing.append(calc(num_nn, avg_NNs_nosing))
                        norm_by_med_num_NNs.append(calc(num_nn, med_NNs))
                        norm_by_med_num_NNs_nosing.append(calc(num_nn, med_NNs_nosing))
                
        r,c=0, 0
        for pred, typ in zip([norm_by_num_vdms, norm_by_direct_vdms,
            norm_by_avg_num_NNs, norm_by_avg_num_NNs_nosing,
            norm_by_med_num_NNs, norm_by_med_num_NNs_nosing], 
            ['num vdms', 'direct vdms', 'avg # NN', 'avg # NN no single',
            'med # NN', 'med # NN no single']):
            print(norm_by_num_vdms)
            pearson_r = stats.pearsonr(all_ddgs, pred)
            print(pearson_r, 'pearson')
            spearman_r = stats.spearmanr(all_ddgs, pred)
            print(spearman_r, 'spearman')
            axarr[r,c].scatter(pred,all_ddgs,marker='o',edgecolors='gray',
                    lw=1)
            axarr[r,c].set_title(typ)
            if c==1:
                c=0
                r+=1
            else:
                c+=1
        plt.suptitle(calc_method)
        plt.show()
