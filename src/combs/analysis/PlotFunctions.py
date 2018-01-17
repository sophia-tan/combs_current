__all__ = ['PlotFunctions']

import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
import seaborn as sns
import pickle as pkl
#from AAcodes import *
import traceback

def unique_aa(aa_list):
    '''Takes a list of all amino acids and gets the unique residues.
    Makes a dict that assigns all the amino acids in this group to a number 
    to assign AA to a row or col # of matrix for heatmap.
    Returns dict and labels'''

    aa_list = set(aa_list)
    aa_list = sorted(aa_list)    
    
    aa_dict = {}
    value = 0
    aa_labels = []
    for aa in aa_list: 
        if aa != 'CSO' and aa != 'MSE' and aa != 'SEC' and aa != 'SEP' and aa != 'TPO':
            aa_dict[aa] = value
            value += 1
            aa_labels.append(aa)
    return aa_dict, aa_labels

def plot_heatmap_ifg_vdm(group1, group2, percentage=False):
    '''Plots heatmap for correlation between vdm res and iFG res
    Inputs: data_df columns for the 2 groups of amino acids to be plotted. 
    Param percentage: decide whether to scale as frequency/percentage or not
    Ex: group1 = vdm, group2 = ifgres'''
    
    # get aa dict and labels for both groups 
    group1dict, group1labels = unique_aa(group1)
    group2dict, group2labels = unique_aa(group2)

    # set up matrix to store frequencies of iFG and vdm (iRes) residues 
    matrix = np.zeros([len(group2dict),len(group1dict)])
    
    # fill in matrix with frequencies of finding vdm-ifg interaction  
    for index in range(len(group1)):
        try: # will only work for amino acids in the dictionary, not the weird AAs
            res1 = group1[index]
            res1 = group1dict[res1] # get index of this vdm residue for the matrix 
            res2 = group2[index]
            res2 = group2dict[res2] # get index of this ifg residue for the matrix
            matrix[res2,res1] += 1
        except:
            pass
    
    if percentage: 
        # convert frequencies in matrix to probabilities/percentages 
        matrix_max = matrix.max()
        for row_ix, row in enumerate(matrix):
            for col_ix, c in enumerate(row):
                matrix[row_ix,col_ix] = matrix[row_ix,col_ix] / matrix_max * 100
    
    return matrix, group1labels, group2labels, group1dict, group2dict

def parse_resnums(resnums):
    ''' Takes a list of resnums in string form and converts to list'''
    resnums_list = []
    if resnums.startswith('-'):
        while resnums.startswith('-'):
            resnums_list.append(int(resnums[:2]))
            resnums = resnums[2:]
    while len(resnums) > 0:
        resnums_list.append(int(resnums[0]))
        resnums = resnums[1:]
    return resnums_list    

def inc_counts_dict(counts, num):
    '''Takes a num (bin) and increases the count for that in the counts_dict
    Inputs: counts dictionary, num (bin)'''
    counts[int(num)] += 1
    return counts

def inc_heatmap_dict(heatmapdict, num, vdm0, vdmi):
    '''Takes a num (bin) and goes into that bin for the dict, and adds the 
    vdm0 and vdmi identities (from appropriate csv file!) to the vdm0 and vdmi lists
    Inputs: heatmapdict, num (bin, int), vdm0 (str), vdmi (str)'''
    
    if vdm0 != 'm' and vdmi != 'm' and vdm0!= 'X' and vdmi != 'X':
        heatmapdict[num]['vdm0'].append(vdm0)
        heatmapdict[num]['vdmi'].append(vdmi)
    return heatmapdict

def vdm0_vdmi_names(row, num):
    '''Determines vdm0 and vdmi (where vdm 'i' in bin #)
    identities from a row of a df (ex: of dist hbond vdms)
    containing columns from vdm_pdb_info csv.  
    Uses num (bin) relative resnum to get the correct amino acid from the vdm string
    Inputs: row, num (bin)'''
    # vdm_seq_resnums there are 3 that can't be parsed bc they look like '-1019-4-3-2-1012345'
    vdm_seq_resnums= parse_resnums(row['rel_resnums_x']) # for everything in vdmstring
    vdmfrag_seq = row['sequence_vdm']
    vdmi = get_aa_of_relresnum_from_string(num, vdm_seq_resnums, vdmfrag_seq)
    vdm0 = get_aa_of_relresnum_from_string(0, vdm_seq_resnums, vdmfrag_seq)
    return vdm0, vdmi

def get_aa_of_relresnum_from_string(num, resnum_string, seq_string):
    '''Gets index of a relresnum from fragresnumstring and then pull the AA from 
    the fragment sequence string. Useful for getting the 'ith' 
    resnum from vdmstring to plot correlation b/n vdm0 and i'th vdm
    identities for the 'ith' interaction bin'''
    
    # get index of the num in the context of the whole vdm fragment
    frag_ix = resnum_string.index(num)
    vdm = seq_string[frag_ix]
    return vdm

def normalize_column(matrix, normalization='max'):
    ''' normalize by dividing by col SUM or MAX. max will make all columns on same scale '''
    for colnum in range(len(matrix[0])):
        col = matrix[:,colnum]
        col_sum = np.sum(col)
        col_max = np.max(col)
        for ix_row, row in enumerate(col):
            if normalization == 'sum':
                matrix[ix_row][colnum] = matrix[ix_row][colnum] / col_sum
            if normalization == 'max':
                matrix[ix_row][colnum] = matrix[ix_row][colnum] / col_max

    return matrix

def normalize_whole_db(matrix, group1labels, group1dict):
    ''' need group1dict to correlate AA to a col # in matrix '''
    for colnum in range(len(matrix[0])):
        # take colnum and find AA it correspond to to look it up in lookup table
        for k,v in group1dict.items():
            if v == colnum:
                aa = k
                aa = one_to_three[aa]
        
        ############# NEED TO CHANGE THIS #######################
        AAi_database_lookup = pkl.load(open('fake.pkl','rb'))
        
        raw_counts_dict, freq_dict = AAi_database_lookup
        aafreq = freq_dict[aa] # lookup table loaded from analysis.Analysis
        col = matrix[:,colnum]
        col_sum = np.sum(col)
        #for ix_row, row in enumerate(col):
        #    matrix[ix_row][colnum] = matrix[ix_row][colnum] / col_sum
    return matrix

def plot_heatmap_allAA_allbins(heatmapdict, rows, cols):
    '''cols should always be > rows'''
    
    f, axarr = plt.subplots(2,3,figsize=(12,6.5))
    fig_row = 0
    fig_col = 0
    
    for relresnum, bin_dict in heatmapdict.items():
        if relresnum > 0:
            x = heatmapdict[relresnum]
            group1 = x['vdm0']
            group2 = x['vdmi']
            axarr[fig_row, fig_col].set_title('Position %s' % relresnum)
            matrix, group1labels, group2labels, group1dict, group2dict = plot_heatmap_ifg_vdm(group1, group2)

            #matrix = normalize_column(matrix)
            #matrix = normalize_whole_db(matrix, group1labels, group1dict)
            
            ax = sns.heatmap(matrix, xticklabels=group1labels,
                        yticklabels=group2labels, cbar_kws={'orientation': 'horizontal'}, ax=axarr[fig_row, fig_col])
            
            if fig_col < 2:
                fig_col += 1
            else:
                fig_col = 0
                fig_row += 1
    
    plt.xticks(rotation=0)
    plt.suptitle('Local H bonds in vdm fragment - raw counts')
    f.text(0.5, 0.04, '0th AA vdm', ha='center')
    f.text(0.04, 0.5, 'i\'th AA vdm', va='center', rotation='vertical')
    plt.yticks(rotation=0)
    plt.show()


def plot_single_heatmap(heatmapdict, relresnum, ss=None):
    x = heatmapdict[relresnum]
    group1 = x['vdm0']
    group2 = x['vdmi']
    matrix, group1labels, group2labels, group1dict, group2dict = plot_heatmap_ifg_vdm(group1, group2)

    matrix = normalize_column(matrix)
    #matrix = normalize_whole_db(matrix, group1labels, group1dict)
    
    ax = sns.heatmap(matrix, xticklabels=group1labels,
                yticklabels=group2labels, cbar_kws={'orientation': 'horizontal'})
    
    plt.xticks(rotation=0)
    plt.suptitle('Local H bonds in vdm fragment, position i = %s, all helix' %relresnum)
    plt.yticks(rotation=0)
    plt.show()


def plot_seq_proximity(hbond_counts):
    plt.bar(list(hbond_counts.keys()), hbond_counts.values())
    plt.title('Sequence proximity preference for H-bonds')
    plt.xlabel('relative seq #')
    plt.ylabel('counts')
    plt.show()

def get_vdmfrag_dssp(row, num):
    ''' Get vdm frag of ss between 0th vdm and ith vdm'''
    vdm_frag_dssp = row['sec_struct_dssp_vdm']
    vdm_seq_resnums= parse_resnums(row['rel_resnums_x']) # for everything in vdmstring
    assert len(vdm_frag_dssp) == len(vdm_seq_resnums)
    dssp_string = ''
    start = min(0,num)
    end = max(0,num)
    while start <= end:
        frag_ix = vdm_seq_resnums.index(start)
        ss = vdm_frag_dssp[frag_ix]
        dssp_string += ss
        start += 1
    return dssp_string

def get_vdmfrag_ss(row, num):
    ''' Get vdm frag of ss between 0th vdm and ith vdm'''
    vdm_frag_dssp = row['sec_struct_dssp_vdm']
    vdm_seq_resnums= parse_resnums(row['rel_resnums_x']) # for everything in vdmstring
    assert len(vdm_frag_dssp) == len(vdm_seq_resnums)
    dssp_string = ''
    start = min(0,num)
    end = max(0,num)
    while start <= end:
        frag_ix = vdm_seq_resnums.index(start)
        ss = vdm_frag_dssp[frag_ix]
        dssp_string += ss
        start += 1
    return dssp_string

def inc_ss_df_Hovmoller(df, row, num):
    ss_string = get_vdmfrag_Hovmoller(row, num)
    
    if abs(num) == 1:
        if ss_string == 'HH':
            ss = 'all helix'
        if ss_string == 'SH':
            ss = 'Nterm strand-helix'
        if ss_string == 'HT':
            ss = 'helix-Cterm turn'
        if ss_string == 'HS':
            ss = 'helix-Cterm strand'
        if ss_string == 'SS':
            ss = 'all strand'
        if ss_string == 'TT':
            ss = 'all loop'
        try:
            df.ix[num, ss] += 1
        except:
            pass
        
    elif abs(num) == 2:
        if ss_string == 'HHH':
            ss = 'all helix'
        if ss_string == 'SHH':
            ss = 'Nterm strand-helix'
        if ss_string == 'HHT':
            ss = 'helix-Cterm turn'
        if ss_string == 'HHS':
            ss = 'helix-Cterm strand'
        if ss_string == 'SSS':
            ss = 'all strand'
        if ss_string == 'TTT':
            ss = 'all loop'
        try:
            df.ix[num, ss] += 1
        except:
            pass

    elif abs(num) == 3:
        if ss_string == 'HHHH':
            ss = 'all helix'
        if ss_string == 'SHHH':
            ss = 'Nterm strand-helix'
        if ss_string == 'HHHT':
            ss = 'helix-Cterm turn'
        if ss_string == 'HHHS':
            ss = 'helix-Cterm strand'
        if ss_string == 'SSSS':
            ss = 'all strand'
        if ss_string == 'TTTT':
            ss = 'all loop'
        try:
            df.ix[num, ss] += 1
        except:
            pass
            
    elif abs(num) == 4 or abs(num) == 5: # treat 4 and 5 the same
        if ss_string[:5] == 'HHHHH':
            ss = 'all helix'
        if ss_string[:5] == 'SHHHH':
            ss = 'Nterm strand-helix'
        if ss_string[:5] == 'HHHHT':
            ss = 'helix-Cterm turn'
        if ss_string[:5] == 'HHHHS':
            ss = 'helix-Cterm strand'
        if ss_string[:5] == 'SSSSS':
            ss = 'all strand'
        if ss_string[:5] == 'TTTTT':
            ss = 'all loop'
        try:
            df.ix[num, ss] += 1
        except:
            pass
        
    
    else:
        print(ss_string)
    
    try:
        return df, ss
    except:
        return df

def inc_ss_df_dssp(df, row, num):
    helix = 'GHI'
    strand = 'EB'
    loop = 'STC'
    
    ss_string = get_vdmfrag_dssp(row, num)
    ss_dict = {'helix': 0, 'strand': 0, 'loop': 0, '-': 0} # keeps count of ss in this dssp string
    
    for letter in ss_string:
        if letter in helix:
            ss_dict['helix'] += 1
        elif letter in strand:
            ss_dict['strand'] += 1
        elif letter in loop:
            ss_dict['loop'] += 1
        elif letter == '-':
            ss_dict['-'] += 1

    len_string = len(ss_string)
    for k, v in ss_dict.items():
        if v == len_string:
            df.ix[num, k] += 1
    return df

def plot_ss_seq_proximity(hbond_counts):
    plt.bar(list(hbond_counts.keys()), hbond_counts.values())
    plt.title('Sequence proximity preference for H-bonds')
    plt.xlabel('relative seq #')
    plt.ylabel('counts')
    plt.show()
    
def get_phi_psi_list(bb_string):
    ''' Helper function to convert phi/psi's from string 
    to list'''
    phi_psi_list = []
    bb = bb_string.lstrip('(')
    bb = bb.rstrip(')')
    bb = bb.split(') (')
    for ix, res in enumerate(bb):
        res = res.split(' ')
        phi_psi = []
        for x in res:
            try:
                phi_psi.append(float(x))
            except:
                phi_psi.append(None)
        phi_psi_list.append(phi_psi)
    return phi_psi_list

def get_vdmfrag_Hovmoller(row, num):
    phipsi_dict = {}
    phipsi_dict['H'] = [[-180, 0], [-100,45]] # [phi, psi]
    phipsi_dict['S'] = [[-180,-45], [45,225]]
    phipsi_dict['T'] = [[0, 180], [-90,90]] # this is actually called "turn" in the Hovmoller 2002 paper

    bb_string = row['sec_struct_phi_psi_frag']
    frag = get_phi_psi_list(bb_string)
    vdm_seq_resnums= parse_resnums(row['rel_resnums_x']) # for everything in vdmstring
    
    ss_list = '' # full vdm frag, later get just a subset of it from 0th-ith or ith-0th position
    if len(frag) == len(vdm_seq_resnums): # ignore the 4 instances where it doesn't correlate
        for ix, res in enumerate(frag):
            assigned = 0
            phi, psi = res[0], res[1]
            for k, v in phipsi_dict.items():
                try:
                    if phi >= v[0][0] and phi <= v[0][1]:
                        if psi >= v[1][0] and psi <= v[1][1]:
                            ss_list += k
                            assigned += 1

                except: # will catch None
                    if phi == None:
                      if psi >= v[1][0] and psi <= v[1][1]:
                          ss_list += k
                          assigned += 1
                    if psi == None:
                        if phi >= v[0][0] and phi <= v[0][1]:
                            ss_list += k
                            assigned += 1
                    else:
                        pass
            # if ss wasn't added to ss_list, mark as '-'
            if assigned == 0:
                ss_list += '-'
    
    # get subset that belongs to 0-ith or ith-0 frag. 
    ss_assignment = ''
    start = min(0,num)
    end = max(0,num)
    if ss_list != '':
        while start <= end:
            frag_ix = vdm_seq_resnums.index(start)
            ss = ss_list[frag_ix]
            ss_assignment += ss
            start += 1
    return ss_assignment

def generate_heatmapdict():
    # store AA identity info to plot heatmap
    heatmapdict = {} # key = bin (same as above),val = dict where (key = vdm0 or vdmi, val = list of those)
    for i in [-5,-4,-3,-2,-1,1,2,3,4,5]: # fill in dict with keys
        heatmapdict[i] = {}
        heatmapdict[i]['vdm0'] = []
        heatmapdict[i]['vdmi'] = []
    return heatmapdict
    
