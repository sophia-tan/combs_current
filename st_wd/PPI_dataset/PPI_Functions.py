import pandas as pd

def generate_resdict(rows):
    '''generate dictionary of WT residues that
    were mutated. one dict for each PDB'''
    res_dict = {} 
    mutations = rows['Mutations']
    exp_ddg = rows['Experimental ddG']
    for mut,ddg in zip(mutations,exp_ddg):
        mut = mut.split(';')
        mut = [x.split(' ') for x in mut]
        for m in mut:
            if m[1][-1]!='A':
                # don't use mutations if they're not to alanine
                pass 
            else:
                chID = m[0]
                resname = m[1][0]
                resnum = ''.join(filter(str.isdigit,m[1]))
                if m[1][-2].isdigit() == True:
                    # make sure the resnum is consecutive
                    assert str(resnum) in m[1] 
                    res_dict[chID+resnum] = [resname, m[1][-1], mut]
                else: # resnum has insertion code
                    resnum += m[1][-2]
                    assert str(resnum) in m[1] 
                    res_dict[chID+resnum] = [resname, m[1][-1], mut]
                # dict value is [orig res, mutated res]
    return res_dict

def generate_data_df(rows):
    '''generate dictionary of mutations data (including combination
    mutations and its experimental data. one dict for each PDB'''
    list_for_df = []
    mutations = rows['Mutations']
    exp_ddg = rows['Experimental ddG']
    for mut,ddg in zip(mutations,exp_ddg):
        mut = mut.split(';')
        mut = [x.split(' ') for x in mut]
        ala_mut = [True if (m[1][-1]=='A') else False for m in mut]
        if False not in ala_mut: # all combo mutations are to ala
            list_for_df.append(
                pd.Series([mut,ddg],index=['Mutation','expt val']))
    return pd.DataFrame(list_for_df)
