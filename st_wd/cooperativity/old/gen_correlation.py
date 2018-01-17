import itertools
import numpy as np

resname_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                   'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                   'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                   'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
                   }

def correlation(df):
    ''' Returns 2 dictionaries: (1) scores, and (2) raw values for num_obs and num_exp'''
    scores = {}
    obs_exp = {}
    df = df[df.apply(get_resnames, axis=1)]
    vdm_gr = df.groupby('iFG_count')
    vdm_gr_agg_pairs = vdm_gr['resname'].agg([num_pairs, to_pairs])
    vdm_gr_agg_pairs_gte1 = vdm_gr_agg_pairs[vdm_gr_agg_pairs['num_pairs'] > 0]
    n_pairs_total = vdm_gr_agg_pairs_gte1['num_pairs'].sum()
    resns = set(df.groupby('resname').groups)
    print(n_pairs_total,'n pairs total')
    
    ### make a dictionary that gives the freq of finding a given AA in a pair ###
    ### this will be used to calculate the expectation num                    ###
    ### need to also account for if this AA is in a pair with itself!         ###
    freq_dict = {}
    def count_res(row):
        c = np.sum(1 for pair in row['to_pairs'] if res in pair)
        dbl = np.sum(1 for pair in row['to_pairs'] if res in pair[0] and res in pair[1])
        return c+dbl
    for res in resns:
        count = vdm_gr_agg_pairs_gte1.apply(count_res, axis=1).sum()
        freq_dict[res] = count
    totalres = sum(list(freq_dict.values()))
    # convert counts to frequencies. values is now a list where first element
    # is the raw counts, and second element is the frequency
    for key, value in freq_dict.items():
        freq_dict[key] = [value, value/totalres]

        
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
            n_exp_res1_res2_pairs = n_pairs_total * freq_dict[res1][1] * freq_dict[res2][1] * 2
        elif res1 == res2:
            n_exp_res1_res2_pairs = n_pairs_total * freq_dict[res1][1] * freq_dict[res2][1]

        scores[(res1, res2)] = -np.log10(n_obs_res1_res2_pairs / n_exp_res1_res2_pairs)
        obs_exp[(res1, res2)] = [n_obs_res1_res2_pairs, n_exp_res1_res2_pairs, freq_dict[res1][1], freq_dict[res2][1]]
        print(res1, res2, obs_exp[(res1, res2)])

    return scores, obs_exp


def num_pairs(x):
    # return len(list(itertools.combinations(sorted(list(x)), 2)))
    # use permutations because
    # using pair[0] and pair[1] in lines 33,36
    return len(list(itertools.combinations(sorted(list(x)), 2)))

def to_pairs(x):
    return list(itertools.combinations(sorted(list(x)), 2))

def get_resnames(row):
    if row['resname'] in set(resname_dict.keys()):
        return True
    else:
        return False

######## START CODE HERE ###########
from sys import argv 
import pandas as pd, pickle as pkl

script, vdm_pdb_info_csv = argv
df = pd.read_csv(vdm_pdb_info_csv)

ifg = vdm_pdb_info_csv.split('_')[0]
corr = correlation(df)
pkl.dump(corr, open('ifg_correlation.pkl'%ifg,'wb'))

