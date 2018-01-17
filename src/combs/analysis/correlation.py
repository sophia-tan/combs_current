import itertools
import numpy as np

resname_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                   'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                   'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                   'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
                   }

# resname_dict = {'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
#                 'THR': 'T', 'ASN': 'N', 'HIS': 'H', 'ARG': 'R', 'TRP': 'W',
#                 'GLU': 'E', 'TYR': 'Y'}


def correlation(df):
    scores = {}
    obs_exp = {}
    df = df[df.apply(get_resnames, axis=1)]
    vdm_gr = df.groupby('iFG_count')
    vdm_gr_agg_pairs = vdm_gr['resname'].agg([num_pairs, to_pairs])
    vdm_gr_agg_pairs_gte1 = vdm_gr_agg_pairs[vdm_gr_agg_pairs['num_pairs'] > 0]
    n_pairs_total = vdm_gr_agg_pairs_gte1['num_pairs'].sum()
    resns = set(df.groupby('resname').groups)
    print(resns)
    print(n_pairs_total)

    for res1, res2 in itertools.combinations_with_replacement(sorted(list(resns)), 2):

        def count_pairs(row):
            return row['to_pairs'].count((res1, res2))

        def count_res1_in_pairs(row):
            return np.sum(1 for pair in row['to_pairs'] if res1 in pair[0])

        def count_res2_in_pairs(row):
            return np.sum(1 for pair in row['to_pairs'] if res2 in pair[1])

        n_obs_res1_res2_pairs = vdm_gr_agg_pairs_gte1.apply(count_pairs, axis=1).sum()

        n_res1_in_pairs = vdm_gr_agg_pairs_gte1.apply(count_res1_in_pairs, axis=1).sum()

        n_res2_in_pairs = vdm_gr_agg_pairs_gte1.apply(count_res2_in_pairs, axis=1).sum()

        n_exp_res1_res2_pairs = n_pairs_total * (n_res1_in_pairs / n_pairs_total) * (n_res2_in_pairs / n_pairs_total)

        scores[(res1, res2)] = -np.log10(n_obs_res1_res2_pairs / n_exp_res1_res2_pairs)
        obs_exp[(res1, res2)] = [n_obs_res1_res2_pairs, n_exp_res1_res2_pairs, n_res1_in_pairs, n_res2_in_pairs]

    return scores, obs_exp


def num_pairs(x):
    # return len(list(itertools.combinations(sorted(list(x)), 2)))
    # use permutations because
    # using pair[0] and pair[1] in lines 33,36
    return len(list(itertools.permutations(sorted(list(x)), 2)))

def to_pairs(x):
    # return list(itertools.combinations(sorted(list(x)), 2))
    return list(itertools.permutations(sorted(list(x)), 2))

def get_resnames(row):
    if row['resname'] in set(resname_dict.keys()):
        return True
    else:
        return False

