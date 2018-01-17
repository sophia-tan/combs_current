__all__ = ['Analyze', 'refine_df', 'combine_bb_sc']
import time
import pandas as pd
import numpy as np
import os
import prody as pr
from ..apps.constants import one_letter_code
from .. import cluster
#from ..cluster.Interactamer import *
import traceback
import pickle as pkl

class Analyze:
    def __init__(self, csv_directory):
        self._directory = csv_directory
        if self._directory[-1] != '/':
            self._directory += '/'
        self.ifg_atom_density = None
        self.ifg_contact_vdm = None
        self.ifg_hbond_water = None
        self.ifg_ca_hbond_vdm = None
        self.ifg_contact_water = None
        self.ifg_pdb_info = None
        self.ifg_contact_ligand = None
        self.ifg_hbond_ligand = None
        self.vdm_pdb_info = None
        self.ifg_contact_metal = None
        self.ifg_hbond_vdm = None
        self.vdm_sasa_info = None
        for file in (file for file in os.listdir(self._directory) if file[0] != '.'):
            if 'ifg_atom_density' in file:
                try:
                    self.ifg_atom_density = pd.read_csv(self._directory + file)
                    self._header = file.split('ifg_atom_density')[0]
                except:
                    print('ifg_atom_density not loaded')
            elif 'ifg_contact_vdm' in file:
                try:
                    self.ifg_contact_vdm = pd.read_csv(self._directory + file)
                    self._header = file.split('ifg_contact_vdm')[0]
                except:
                    print('ifg_contact_vdm not loaded')
            elif 'ifg_hbond_water' in file:
                try:
                    self.ifg_hbond_water = pd.read_csv(self._directory + file)
                    self._header = file.split('ifg_hbond_water')[0]
                except:
                    print('ifg_hbond_water not loaded')
            elif 'ifg_ca_hbond_vdm' in file:
                try:
                    self.ifg_ca_hbond_vdm = pd.read_csv(self._directory + file)
                    self._header = file.split('ifg_ca_hbond_vdm')[0]
                except:
                    print('ifg_ca_hbond_vdm not loaded')
            elif 'ifg_contact_water' in file:
                try:
                    self.ifg_contact_water = pd.read_csv(self._directory + file)
                    self._header = file.split('ifg_contact_water')[0]
                except:
                    print('ifg_contact_water not loaded')
            elif 'ifg_pdb_info' in file:
                try:
                    self.ifg_pdb_info = pd.read_csv(self._directory + file)
                    self._header = file.split('ifg_pdb_info')[0]
                except:
                    print('ifg_pdb_info not loaded')
            elif 'ifg_contact_ligand' in file:
                try:
                    self.ifg_contact_ligand = pd.read_csv(self._directory + file)
                    self._header = file.split('ifg_contact_ligand')[0]
                except:
                    print('ifg_contact_ligand not loaded')
            elif 'ifg_hbond_ligand' in file:
                try:
                    self.ifg_hbond_ligand = pd.read_csv(self._directory + file)
                    self._header = file.split('ifg_hbond_ligand')[0]
                except:
                    print('ifg_hbond_ligand not loaded')
            elif 'vdm_pdb_info' in file:
                try:
                    self.vdm_pdb_info = pd.read_csv(self._directory + file)
                    self._header = file.split('vdm_pdb_info')[0]
                except:
                    print('vdm_pdb_info not loaded')
            elif 'ifg_contact_metal' in file:
                try:
                    self.ifg_contact_metal = pd.read_csv(self._directory + file)
                    self._header = file.split('ifg_contact_metal')[0]
                except:
                    print('ifg_contact_metal not loaded')
            elif 'ifg_hbond_vdm' in file:
                try:
                    self.ifg_hbond_vdm = pd.read_csv(self._directory + file)
                    self._header = file.split('ifg_hbond_vdm')[0]
                except:
                    print('ifg_hbond_vdm not loaded')
            elif 'vdm_sasa_info' in file:
                try:
                    self.vdm_sasa_info = pd.read_csv(self._directory + file)
                    self._header = file.split('vdm_sasa_info')[0]
                except:
                    print('vdm_sasa_info not loaded')

    def get_distant_vdms(self, seq_distance=10):
        mer = pd.merge(self.ifg_pdb_info, self.vdm_pdb_info, on='iFG_count', suffixes=('_ifg', '_vdm'))
        def func(row):
            try:
                return np.abs(row['resindex_ifg'] - row['resindex_vdm']) > seq_distance
            except:
                return all(np.abs(np.array([int(ri) for ri in row['resindex_ifg'].split()]) - row['resindex_vdm']) > seq_distance)
        return mer[mer.apply(func, axis=1)]

    def get_lonepair_imidazole_ifgs(self, dist_df):
        def lonepair_atoms(row):
            try:
                relres = row['rel_resnums']
                relres = relres.replace('-', '')
                l = range(len(relres))
                inds = [i for i in l if relres[i] == '0']
                if inds:
                    s = row['iFG_atom_names']
                    s = s.split(') (')
                    s = [k.replace('(', '') for k in s]
                    s = [k.replace(')', '') for k in s]
                    s = set([k for m in [s[i].split() for i in inds] for k in m])
                    if ({'NE2'}.issubset(s) and not {'HE2'}.issubset(s)):
                        return True
                    elif ({'ND1'}.issubset(s) and not {'HD1'}.issubset(s)):
                        return True
                    else:
                        return False
                else:
                    return False
            except:
                return False
        hbond_vdms = self.ifg_hbond_vdm[self.ifg_hbond_vdm.apply(lonepair_atoms, axis=1)]
        ca_hbond_vdms = self.ifg_ca_hbond_vdm[self.ifg_ca_hbond_vdm.apply(lonepair_atoms, axis=1)]
        merged = pd.concat([hbond_vdms,ca_hbond_vdms])
        merged = merged[['iFG_count','vdM_count']]
        # merge again with dist_vdms 
        merged = pd.merge(dist_df, merged, on=['iFG_count','vdM_count'])
        return merged
    
    def make_csv(self, df, path=None):
        path = path or self._directory
        df.to_csv(path + self._header + 'distant_vdm.csv', index=False)

    def get_hbonding_vdms(self, df, hbond_seq_pattern, mode='sidechain'):
        mer = pd.merge(df, self.ifg_hbond_vdm, on=['iFG_count', 'vdM_count'])
        mer = mer[mer['rel_resnums'].str.contains(hbond_seq_pattern)]
        if mode == 'sidechain':
            def func(row):
                index = [i for i, rrn in enumerate([n for n in row['rel_resnums'] if n != '-']) if rrn == '0']
                # print(row['rel_resnums'])
                # print(row['vdM_atom_names'])
                # print(row['vdM_atom_names'].strip('()').split(') ('))
                # print(row[['iFG_count','vdM_count']])
                # print(index)
                return any([any(
                    {y for names in [row['vdM_atom_names'].strip('()').split(') (')[ind]] for y in names.split() if
                     y not in ['N', 'O', 'H', 'OXT', 'H1', 'H2', 'H3']}) for ind in index])
            mer = mer[mer.apply(func, axis=1)]
        elif mode == 'mainchain':
            def func(row):
                index = [i for i, rrn in enumerate([n for n in row['rel_resnums'] if n != '-']) if rrn == '0']
                return any([any(
                    {y for names in [row['vdM_atom_names'].strip('()').split(') (')[ind]] for y in names.split() if
                     y in ['N', 'O', 'H']}) for ind in index])
            mer = mer[mer.apply(func, axis=1)]
        else:
            raise NameError("arg *mode* must be 'sidechain' or 'mainchain'")

        return mer #df[(df.iFG_count.isin(mer.iFG_count)) & (df.vdM_count.isin(mer.vdM_count))]

    @staticmethod
    def group_vdms_by_aa(df):
        try:
            return df.groupby('resname')
        except:
            return df.groupby('resname_vdm')

    def parse_vdms(self, df, path_to_vdm=None):
        path = path_to_vdm or self._directory.split('csv')[0] + 'vdM/'
        listdir = os.listdir(path)
        filename_end = '_'.join(listdir[0].split('_')[4:])
        listtags = [[int(file.split('_')[1]), int(file.split('_')[3])] for file in listdir]
        # return [pr.parsePDB(path + 'iFG_' + str(tag[0]) + '_vdM_' + str(tag[1]) + '_' + filename_end)
        #         for tag in df[['iFG_count', 'vdM_count']].values.tolist() if tag in listtags]
        parsed_vdms = []
        for tag in df[['iFG_count', 'vdM_count']].values.tolist():
            if tag in listtags:
                try:
                    parsed_vdms.append(pr.parsePDB(path + 'iFG_' + str(tag[0]) + '_vdM_' + str(tag[1]) + '_'
                                                   + filename_end))
                except Exception:
                    traceback.print_exc()
        return parsed_vdms


    def parse_vdms_by_aa(self, df, subset=None, path_to_vdm=None):
        """subset is a list of residue names"""
        gr_df = self.group_vdms_by_aa(df)
        parsed_by_aa = {}
        subset = subset or gr_df.groups
        for group in set(gr_df.groups).intersection(set(subset)):
            if group in one_letter_code.keys():
                parsed_by_aa[group] = self.parse_vdms(gr_df.get_group(group), path_to_vdm)
        return parsed_by_aa

def parse_vdms_for_interactamers(path_to_csvs):
    an = Analyze(path_to_csvs)
    dist_vdms = an.get_distant_vdms(seq_distance=10)
    hb_dist_vdms = an.get_hbonding_vdms(dist_vdms, '0', mode='sidechain')
    # hb_dist_vdms = pd.merge(dist_vdms, an.ifg_hbond_vdm, on=['iFG_count', 'vdM_count'])
    # hb_dist_vdms = hb_dist_vdms[hb_dist_vdms['rel_resnums'].str.contains('0')]
    # mer2 = pd.merge(hb_mer, an.ifg_pdb_info, on=['iFG_count'])
    # bo = mer2.apply(lambda x: any(
    #     {y for each in x['vdM_atom_names'].strip('()').split(') (') for y in each.split() if y not in 'NOH'}), axis=1)
    return an.parse_vdms_by_aa(hb_dist_vdms)
    # return an.parse_vdms_by_aa(mer2[bo])

def frag_match(seq1, seq2):
    '''Gives % identity match between 2 seq strings of same length'''
    try:
        assert len(seq1) == len(seq2)
        match = 0
        total = len(seq1)
        for i in range(total):
            if seq1[i] == seq2[i]:
                match += 1
        perc = match/total*100
    except:
        perc = 0
    return perc

def repeat_indices():
    '''delete this, was only for troubleshooting'''
    x = []
    with open('leu_leu_indices.txt') as f:
        for line in f:
            x.append(line.strip())
    x = [int(i) for i in x]
    return x

def remove_repeat_proteins(orig_df):
    '''Takes df (ex: dist h-bond vdms). Group by pdbcode, resname_ifg, and
    resname_vdm, because if all 3 are the same, then they're possible repeats.
    ''' 
    df = orig_df[['pdb','resname_ifg', 'resname_vdm', 'sequence_vdm', 'sec_struct_dssp_vdm']]
    grouped = df.groupby(by=['pdb', 'resname_ifg', 'resname_vdm'])
    print(len(orig_df), 'num vdms in original df')
    drop_indices = []
    for pdb_ifg_vdm, group in grouped:
        if len(group) > 1:
            for poss_ix, poss in group.iterrows():
                for ix, row in group.loc[poss_ix+1:].iterrows(): # +1 is to avoid ix and poss_ix 
                # from being in the same row
                    seq_match = frag_match(row['sequence_vdm'], poss['sequence_vdm'])
                    bb_match = frag_match(row['sec_struct_dssp_vdm'], poss['sec_struct_dssp_vdm'])
                    if bb_match > 90 and seq_match > 90:
                        drop_indices.append(ix)
                        group = group.drop(ix, axis=0)
    
    orig_df = orig_df.drop(drop_indices, axis=0)
    print(len(orig_df), 'num vdms after repeats are removed')
    return orig_df

#AAi_database_lookup = pkl.load(open('fake.pkl','rb'))


def drop_rare_aa(df):
    drop_indices = []
    for ix, row in df.iterrows():
        rare = ['MSE', 'SEP', 'TPO', 'CSO']
        if row['resname_vdm'] in rare:
            drop_indices.append(ix)
    df.drop(drop_indices, axis=0,inplace=True)
    return df

def refine_df(path_to_csv, seq_dist, threefive, repeats_already_removed=False, repeats_removed_df=None,lonepair_imidazole=False):
    '''filter df by seq_dist, bb/sc, 3.5.
    Argument seq_dist should be float. threefive should be boolean'''
    an = Analyze(path_to_csv)
    if repeats_already_removed==False:
        dist_vdms = an.get_distant_vdms(seq_distance=seq_dist)
        dist_vdms = remove_repeat_proteins(dist_vdms)
    else:
        dist_vdms = repeats_removed_df

    if lonepair_imidazole:
        dist_vdms = an.get_lonepair_imidazole_ifgs(dist_df = dist_vdms)

    # add info about ifg-vdm dist from contacts csv and drop rare aa
    contactsdf = an.ifg_contact_vdm
    dist_vdms = pd.merge(dist_vdms,contactsdf, on=['iFG_count', 'vdM_count'])
    dist_vdms = drop_rare_aa(dist_vdms)
    vdms_bb = dist_vdms[dist_vdms.apply(cluster.Interactamer.has_bb_or_sc, bb_or_sc='bb',threepfive=threefive,axis=1)]
    vdms_sc = dist_vdms[dist_vdms.apply(cluster.Interactamer.has_bb_or_sc, bb_or_sc='sc',threepfive=threefive,axis=1)]

    return dist_vdms, vdms_bb, vdms_sc, an
    
def combine_bb_sc(bb, sc):
    '''useful for using the refine_df fx to get vdms within 3.5A within a separation dist.
    however, that returns separate dfs for bb/sc, so this combines them and takes out repeats'''

    pd.concat([bb, sc])
    concatenated = pd.concat([bb, sc],ignore_index=True)
    nodup=concatenated.drop_duplicates()
    return nodup
    











# implement redundancy removal feature
# code clustering
# add hydrogens to waters and look for hbonds? or look for hbonds with water based on distances of heavy atoms?

##NOTE figure out why vdms are being parsed multiple times.  probably has to do with pandas merging.

# need to decide distant/local based on connectivity, not res num.  resindex seems ok.


# an = combs.Analysis.Analyze('/Users/npolizzi/Projects/combs/results/carboxamide/20170715/csv/')
# dist_vdms = an.get_distant_vdms(seq_distance=20)
# phe_dist_vdms = dist_vdms[dist_vdms['resname_vdm'] == 'PHE']
# sc_phe_dist_vdms = phe_dist_vdms[phe_dist_vdms.apply(lambda row: any(name in ['CG','CD1','CD2','CE1','CE2','CZ'] for name in row['atom_names_vdm'].split()), axis=1)]
# parse_sc_phe_dist_vdms = an.parse_vdms(sc_phe_dist_vdms)
# ifg_coords_phe,phe_vdms =combs.analysis.make_ifg_coords_array(parse_sc_phe_dist_vdms, cb, 'CE1', 'CE2', 'CZ')
# ifg_coords_phe_flip,phe_vdms_flip =combs.analysis.make_ifg_coords_array(parse_sc_phe_dist_vdms, cb, 'CE2', 'CE1', 'CZ')
# pair_dist = spatial.distance.pdist(np.vstack((ifg_coords_phe,ifg_coords_phe_flip)))
# pair_dist_mat = spatial.distance.squareform(pair_dist)
# mem_phe,cen_phe = combs.cluster(pair_dist_mat)
# all_phe_vdms = phe_vdms
# all_phe_vdms.extend(phe_vdms_flip)
# for i,clust in enumerate(mem_phe):
#     os.mkdir('/Users/npolizzi/Desktop/clusters_phe_flip/' + str(i))
#     for m in clust:
#         if m > 2063:
#             writePDB('/Users/npolizzi/Desktop/clusters_phe_flip/' + str(i) + '/' + repr(all_phe_vdms[m]).split()[1] + '_flip.pdb.gz', all_phe_vdms[m])
#         else:
#             writePDB('/Users/npolizzi/Desktop/clusters_phe_flip/' + str(i) + '/' + repr(all_phe_vdms[m]).split()[1] + '.pdb.gz', all_phe_vdms[m])


# import numpy as np
#    ...: import sys
#    ...: sys.path.append('/Users/npolizzi/Projects/combs/src/')
#    ...: import combs
#    ...:
#
# In [2]: an = combs.Analysis.Analyze('/Users/npolizzi/Projects/combs/results/carb
#    ...: oxamide/20170715/csv/')
#
# In [3]: dist_vdms = an.get_distant_vdms(seq_distance=20)
#
# In [4]: nh_bb_dist_vdms = dist_vdms[dist_vdms.apply(lambda row: any(name in ['C'
#    ...: ,'N','CA'] for name in row['atom_names_vdm'].split()), axis=1)]
#
# # In [5]: parse_co_bb_dist_vdms = an.parse_vdms(co_bb_dist_vdms)
#
# ifg_coords_co_bb,co_bb_vdms =combs.analysis.make_ifg_coords_array(parse_co_bb_dist_vdms, cb, '
#     ...: CA', 'O', 'C')
#
# pair_dist_co_bb = spatial.distance.pdist(ifg_coords_co_bb)
#
# pair_dist_co_bb_mat = spatial.distance.squareform(pair_dist_co_bb)
# mem_co_bb,cen_co_bb = combs.cluster(pair_dist_co_bb_mat)






