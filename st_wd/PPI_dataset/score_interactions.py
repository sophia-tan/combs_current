import sys
import matplotlib.pyplot as plt
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
import prody as pr
import numpy as np, pickle as pkl, pandas as pd
from itertools import *
from Scoring import *
from Functions import *
from scipy import stats

for pdb_index in [0]:
#for pdb_index in range(54):
    matches = pkl.load(open('./output_data/pdb_ix_{}_matches.pkl'.
        format(pdb_index),'rb'))
    mutation_data = pkl.load(open('mutation_data.pkl','rb'))
    pdb_datadf = mutation_data[1][int(pdb_index)]

    for rmsd in [.5]:
        for calc_method in ['d']:
            filtered_matches = matches[matches['rmsd']==rmsd]
            '''For each experimental sample, make computational prediction'''
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
            
            for row_ix, row in pdb_datadf.iterrows():
                ddg = row['ddG']
                combo_mutation = str(row['Mutation'])
                predictions = filtered_matches[filtered_matches['mut string'].
                    astype(str) == combo_mutation]
                labels.append(combo_mutation)
                all_ddgs.append(ddg)
                num_nn = predictions['num NN']
                num_vdms = predictions['num vdms']
                direct_vdms = predictions['num directly interacting vdms']
                avg_NNs = predictions['avg num NNs']
                avg_NNs_nosing = predictions['avg num NNs w/o singles']
                med_NNs = predictions['median num NNs']
                med_NNs_nosing = predictions['median num NNs w/o singles']

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
            
            print(rmsd, 'rmsd')

            for pred, title in zip([norm_by_num_vdms, norm_by_direct_vdms], ['num vdms', 'direct vdms']):
                print(title, 'type of normalizing')
                print(calc_method, 'type of scoring')
                pearson_r = stats.pearsonr(all_ddgs, pred)
                print(pearson_r, 'pearson')
                spearman_r = stats.spearmanr(all_ddgs, pred)
                print(spearman_r, 'spearman')
                plt.scatter(pred,all_ddgs,marker='o',edgecolors='gray',
                        lw=1)
                plt.title(title)
                plt.show()
