import numpy as np, pickle as pkl, pandas as pd
import traceback, os

def get_intrxn_geom_info(study, method, resdict=None):
    '''Make list of DF of all interaction geom info.
    Args: study = 'integrin' or 'PPI_dataset' 
          method = which ifg/vdm substructure to analyze geom of'''
    pklf_list = []
    if study == 'AiibB3':
        pkl_dir = '/home/gpu/Sophia/combs/st_wd/AiibB3/output_data/'
        #pkl_dir = '/home/gpu/Sophia/combs/st_wd/K658/output_data/'
        for project in [resdict]:
            for targetres, resn in resdict.items():
                pklf = '{}_matches_{}.pkl'.format(targetres, method)
                if pklf in os.listdir(pkl_dir):
                    try: 
                        new = pkl.load(open(pkl_dir + pklf, 'rb'))
                        matches = pd.concat([matches, new])
                    except:
                        matches = pkl.load(open(pkl_dir + pklf, 'rb'))
                else: 
                    pass
        pklf_list = [matches]

    elif study == 'PPI_dataset':
        pkl_dir = '/home/gpu/Sophia/combs/st_wd/PPI_dataset/output_data/'
        for pdb_index in range(54):
            pklf = 'pdb_ix_{}_matches_{}.pkl'.format(pdb_index, method)
            if pdb_index == 0:
                matches = pkl.load(open(pkl_dir + pklf, 'rb'))
            else:
                new = pkl.load(open(pkl_dir + pklf, 'rb'))
                matches = pd.concat([matches, new])
            pklf_list.append(matches)

    return pklf_list

def get_experimental_df(study):
    if study == 'PPI_dataset': 
        pklfi = '/home/gpu/Sophia/combs/st_wd/PPI_dataset/experimental_data_PPI_dataset.pkl'
        experimental = pkl.load(open(pklfi, 'rb'))
        experimental = mutation_data[1][int(matches_index)]
    
    elif study == 'AiibB3':
        pklfi = '/home/gpu/Sophia/combs/st_wd/AiibB3/experimental_data_AiibB3.pkl'
        experimental = pkl.load(open(pklfi, 'rb'))
    
    return experimental 
    
def get_predictions_for_mutant(row, study, filtered_matches):
    combo_mutation = row['Mutation']
    if study == 'PPI_dataset': 
        mut_str = 'mut string'
    elif study == 'AiibB3':
        mut_str = 'hotspot res'
    
    predictions = filtered_matches[filtered_matches[
       mut_str].astype(str) == str(combo_mutation)]

    return predictions, combo_mutation
    
def get_filtered_matches(matches_df, rmsd, nonmembrane, buried):
    filtered_matches = matches_df[matches_df['rmsd']==rmsd]
    filtered_matches = filtered_matches[filtered_matches['nonmembrane']==nonmembrane]
    filtered_matches = filtered_matches[filtered_matches['buried']==buried]
    return filtered_matches 

def add_one(row):
    return row+0.1
            
def get_label(combo_mutation, study):
    if study == 'PPI_dataset': 
        lbl = [x[1][0] for x in combo_mutation]
    elif study == 'AiibB3':
        lbl = combo_mutation
    return lbl

######################################################################
# Diff methods of calculating score for list of observed counts 'a' 
# and list of expected counts 'b'
######################################################################
def list_to_ndarray(a, b):
    return [np.array(a), np.array(b)]

def calc_score(a,b,score): 
    a, b = list_to_ndarray(a, b)
    # all variations of just rawcounts, no normalization
    if score == 'a':
        return sum(a)
    if score == 'b':
        return sum(np.log10(a))
    if score == 'c':
        return np.log10(sum(a))

    # all variations of arrays 'a' and 'b' treated equally
    if score == 'd':
        return sum(a)/sum(b)
    if score == 'e':
        return sum(a/b)
    if score == 'f':
        return np.log10(sum(a/b))
    if score == 'g':
        return np.log10(sum(a)/sum(b))
    if score == 'h':
        return sum(np.log10(a/b))
    if score == 'i':
        return sum(np.log10(a)/np.log10(b))
    if score == 'j':
        return np.log10(sum(a))/np.log10(sum(b))
    if score == 'k':
        return sum(np.log10(a))/sum(np.log10(b))

    # all variations with log of only array 'a' and not 'b'
    if score == 'l':
        return sum(np.log10(a)/b)
    if score == 'm':
        return np.log10(sum(a))/sum(b)
    if score == 'n':
        return sum(np.log10(a))/sum(b)
    
    # all variations with log of only array 'b' and not 'a'
    if score == 'o':
        return sum(a/np.log10(b))
    if score == 'p':
        return sum(a)/np.log10(sum(b))
    if score == 'q':
        return sum(a)/sum(np.log10(b))


def best_fit_slope(xs,ys):
    '''for plotting'''
    xs,ys=np.array(xs),np.array(ys)
    m = (((np.mean(xs)*np.mean(ys)) - np.mean(xs*ys)) / \
             ((np.mean(xs)*np.mean(xs)) - np.mean(xs*xs)))
    b = np.mean(ys) - m*np.mean(xs)
    return m,b
