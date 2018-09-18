import sys, os
from scipy import stats
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
import pickle as pkl, numpy as np, pandas as pd
import matplotlib.pyplot as plt
sys.path.append('/home/gpu/Sophia/combs/st_wd/Scoring_Function')
from CalculateScores import *
from AiibB3_residues import *

##### is this for PPI_dataset or integrin? ####
study = 'AiibB3'
resdict = integrin_res # from AiibB3_residues.py
# get list of pklfiles that contain interaction geometry
# information from calc_rmsd.py
matches_list = get_intrxn_geom_info(study, 'planar_group_ala_mutants', resdict=resdict)

no_ser = False
#############################################################

### ANALYZING DETAILS ###
def rmsd_filter(method, score_method, nonmembrane, buried, rmsd, matches_list):
    norm_by_num_vdms = []
    norm_by_direct_vdms = []
    norm_by_avg_num_NNs = []
    norm_by_avg_num_NNs_nosing = []
    norm_by_med_num_NNs = []
    norm_by_med_num_NNs_nosing = []
    colors = []
    labels = []
    ddgs_of_no_intrxn = [] # for PPI_dataset
    experimental_values = []

    for matches_index, matches_df in enumerate(matches_list):
    # matches_df could be whole integrin if study=AiibB3, or each pdb if study=PPI_dataset
    
        experimental = get_experimental_df(study)
        filtered_matches = get_filtered_matches(matches_df, rmsd, nonmembrane, buried)

        '''For each experimental mutant, make computational prediction'''
        for row_ix, row in experimental.iterrows():
            expt_val = row['expt val']
            predictions, combo_mutation = get_predictions_for_mutant(row, study, filtered_matches)

            ## ignore where num NN = 0
            #predictions = predictions[predictions['num NN'] > 0]
            predictions = predictions[predictions['switched ifg/vdm?']=='not_switched']
            num_nn = predictions['num NN']
            num_vdms = predictions['num vdms']
            direct_vdms = predictions[
                'num directly interacting vdms']
            avg_NNs = predictions['avg num NNs']
            avg_NNs_nosing = predictions[
                'avg num NNs w/o singles']
            med_NNs = predictions['median num NNs']
            med_NNs_nosing = predictions[
                'median num NNs w/o singles']
                
            #num_nn = predictions['num NN'].apply(add_one)
            #num_vdms = predictions['num vdms'].apply(add_one)
            #direct_vdms = predictions[
            #    'num directly interacting vdms'].apply(add_one)
            #avg_NNs = predictions['avg num NNs'].apply(add_one)
            #avg_NNs_nosing = predictions[
            #    'avg num NNs w/o singles'].apply(add_one)
            #med_NNs = predictions['median num NNs'].apply(add_one)
            #med_NNs_nosing = predictions[
            #    'median num NNs w/o singles'].apply(add_one)
            
            if len(predictions) == 0: 
                ddgs_of_no_intrxn.append(expt_val) 
                num_nn = [0]
                num_vdms = [10000]
                direct_vdms = [10000]
                avg_NNs = [10]
                avg_NNs_nosing = [10]
                med_NNs = [10]
                med_NNs_nosing = [10]
            
            lbl = get_label(combo_mutation, study)
            print(lbl)
            print(predictions[['num NN', 'vdm resi', 'vdm resn', 'num vdms']])
            colors.append('b')
            
            if no_ser == True and lbl == 'A758':
                pass
            else:
                experimental_values.append(expt_val)
                norm_by_num_vdms.append(calc_score(num_nn, num_vdms, score_method))
                norm_by_direct_vdms.append(calc_score(num_nn, direct_vdms, score_method))
                norm_by_avg_num_NNs.append(calc_score(num_nn, avg_NNs, score_method))
                norm_by_avg_num_NNs_nosing.append(calc_score(num_nn, avg_NNs_nosing, score_method))
                norm_by_med_num_NNs.append(calc_score(num_nn, med_NNs, score_method))
                norm_by_med_num_NNs_nosing.append(calc_score(num_nn, med_NNs_nosing, score_method))
            
    r,c=0, 0

    f, axarr = plt.subplots(3,2)
    for pred, typ in zip([norm_by_num_vdms, norm_by_direct_vdms,
        norm_by_avg_num_NNs, norm_by_avg_num_NNs_nosing,
        norm_by_med_num_NNs, norm_by_med_num_NNs_nosing], 
        ['num vdms', 'direct vdms', 'avg # NN', 'avg # NN no single',
        'med # NN', 'med # NN no single']):
        pearson_r = stats.pearsonr(experimental_values, pred)
        if pearson_r[0]**2 > 0.4:
            print('--------------------------')
            print('substructure method', method, typ)
            print('score method', score_method)
            print(buried, nonmembrane, rmsd)
            print(pred)
            print(typ)
            print(experimental_values)
    
            #print(pearson_r)
            print(pearson_r[0]**2, 'pearson')
##      
        #spearman_r = stats.spearmanr(experimental_vals, pred)
        axarr[r,c].scatter(pred,experimental_values,marker='o',
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
    #plt.savefig('./output_data/{}_{}_{}_{}_{}.png'.format(rmsd, score_method, nonmembrane, buried, method))
##
for method in ['planar_group_ala_mutants']:
#for method in ['planar_group_ala_mutants', 'sc']:
    #for score_method in ['d','e','f','g','h','i','j','k','l','m','n','o','p','q']:
    ####for score_method in ['f','g','h','i','j','k']:
    for score_method in ['e']:
        for rmsd in [.5]:
        #for rmsd in [.4,.5,.6,.7]:
            for nonmembrane in [True]:
            #for nonmembrane in [True, False]:
                for buried in [True]:
                #for buried in [True, False]:
                    rmsd_filter(method, score_method, nonmembrane, buried, rmsd, matches_list)
