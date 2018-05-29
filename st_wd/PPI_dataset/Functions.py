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
            if m[1][-2].isdigit() == False or m[1][-1]!='A':
                # 1) ignore weird antibody numbering. bc some are numbered
                #    strangely like PDB 3NPS R60fA and Y60gA
                # 2) don't use mutations if they're not to alanine
                pass 
            else:
                chID = m[0]
                resname = m[1][0]
                resnum = ''.join(filter(str.isdigit,m[1]))
                assert str(resnum) in m[1] # make sure the resnum is consecutive
                res_dict[chID+resnum] = [resname, m[1][-1],ddg]
                # dict value is [orig res, mutated res, ddG]
    return res_dict

def get_pdbdict(pdb_num,variables):
    # this is a hack-y way to get the pdb_dict 
    for k,v in variables:
        # k = name of object
        # v = object 
        if type(v) == list:
            if v[0] == int(pdb_num):
                return v # dict
    
