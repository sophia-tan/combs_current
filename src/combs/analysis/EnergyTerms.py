__all__ = ['make_freqaai_df']

from .. import cluster, analysis, apps
import os
import pickle as pkl, pandas as pd, numpy as np
import itertools

def make_freqaai_df(df):
    '''
    Inputs: df of skip10, no repeats vdms
    Returns: dict of dictionaries. (1) counts for bb vdms. 
    (2) its freq. (3) counts for sc vdms. (4) its freq. 
    Note: some vdms will be counted as both bb and sc.'''
    
    bb_atoms = ['N', 'CA', 'C', 'O', 'OXT']
    df = df[['resname_vdm', 'atom_names_vdm', 'dist_info']]

    # get bb and sc vdms that are within 3.5A
    bb = df[df.apply(cluster.Interactamer.has_bb_or_sc, bb_or_sc='bb', threepfive=True,axis=1)]
    sc = df[df.apply(cluster.Interactamer.has_bb_or_sc, bb_or_sc='sc', threepfive=True,axis=1)]
    num_within_threefive = len(set(list(pd.concat([bb, sc]).index.values)))
    print((num_within_threefive), 'num vdMs within 3.5A')
    print(len(bb), 'num bb interactamers within 3.5A')
    print(len(sc), 'num sc interactamers within 3.5A')
    bbcounts = pd.DataFrame(bb['resname_vdm'].value_counts())
    sccounts = pd.DataFrame(sc['resname_vdm'].value_counts())
    aai_df = pd.merge(sccounts,bbcounts, how='outer',right_index=True, left_index=True, suffixes=('_sc', '_bb'))
    
    # drop rare AAs 
    for ix, row in aai_df.iterrows():
        rare = ['MSE', 'SEP', 'TPO', 'CSO']
        if ix in rare:
            aai_df.drop(ix, axis=0,inplace=True)
    
    # get frequencies from the counts
    for col in aai_df.columns.values:
        interaction_type = col.split('_')[2]
        total = aai_df[col].sum()
        freq_col = pd.Series(aai_df[col]/total, name='sdf')
        aai_df = pd.merge(aai_df, pd.DataFrame(freq_col), right_index=True, left_index=True)
        aai_df = aai_df.rename(index=str, columns={'sdf':'vdm_freq_'+interaction_type})
    return aai_df

def AAi_db_lookup(lookup_dir):
    aa_dict = pkl.load(open(lookup_dir+'AAi_freq/AAi_database_lookups.pkl','rb'))[1]
    return aa_dict

def correlation(df):
    ''' Returns 2 dictionaries: (1) scores, and (2) raw values for num_obs and num_exp'''
    scores = {}
    obs_exp = {}
    df = df[df.apply(get_resnames, axis=1)]
    vdm_gr = df.groupby('iFG_count')
    vdm_gr_agg_pairs = vdm_gr['resname_vdm'].agg([num_pairs, to_pairs])
    vdm_gr_agg_pairs_gte1 = vdm_gr_agg_pairs[vdm_gr_agg_pairs['num_pairs'] > 0]
    n_pairs_total = vdm_gr_agg_pairs_gte1['num_pairs'].sum()
    resns = set(df.groupby('resname_vdm').groups)
    print(n_pairs_total,'n pairs total')
    
    ### make a dictionary that gives the freq of finding a given AA in a pair ###
    ### this will be used to calculate the expectation num                    ###
    ### need to also account for if this AA is in a pair with itself!         ###
        
    def count_obs_pairs(row, a, b):
        # order doesn't matter
        if a != b:
            return row['to_pairs'].count((a, b)) + row['to_pairs'].count((b,a))
        else: 
            return row['to_pairs'].count((a, b))
            
    for res1, res2 in itertools.combinations_with_replacement(sorted(list(resns)), 2):
        n_obs_res1_res2_pairs = vdm_gr_agg_pairs_gte1.apply(count_obs_pairs,axis=1,args=(res1,res2)).sum()

        if res1 != res2:
            # freq_dict is a dict where values is a list: second element is the freq
            # need to multiply by 2 because need to count prob of (AB) AND (BA)
            n_exp_res1_res2_pairs = n_pairs_total * 0.05**2 *2
        elif res1 == res2:
            n_exp_res1_res2_pairs = n_pairs_total * 0.05**2

        scores[(res1, res2)] = -np.log10(n_obs_res1_res2_pairs / n_exp_res1_res2_pairs)
        obs_exp[(res1, res2)] = [n_obs_res1_res2_pairs, n_exp_res1_res2_pairs]
        #print(res1, res2, obs_exp[(res1, res2)])
    return scores, obs_exp

def num_pairs(x):
    # return len(list(itertools.combinations(sorted(list(x)), 2)))
    # use permutations because
    # using pair[0] and pair[1] in lines 33,36
    return len(list(itertools.combinations(sorted(list(x)), 2)))

def to_pairs(x):
    return list(itertools.combinations(sorted(list(x)), 2))

def get_resnames(row):
    if row['resname_vdm'] in set(apps.resname_dict.keys()):
        return True
    else:
        return False

