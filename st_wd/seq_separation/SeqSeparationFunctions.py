import matplotlib.pyplot as plt
import sys                              
import pandas as pd
import pickle as pkl
import itertools

def gen_seq_sep_df(interaction_type):
    seq_effects_lookup = pd.DataFrame(index=[i for i in range(50)], columns=['total counts', 'bb/sc', 'bb/bb', 'sc/sc', 'resnames', 'row_ix'])
    # fill df with zeros
    for ix, row in seq_effects_lookup.iterrows():
        for col in seq_effects_lookup.columns.values:
            if col=='resnames' or col=='row_ix':
                seq_effects_lookup.ix[ix, col] = []
            else:   
                seq_effects_lookup.ix[ix, col] = 0
    return seq_effects_lookup
        
def get_sep(vdms_df, pair):
    vdmresnum_a = vdms_df.ix[pair[0], 'resindex_vdm']
    vdmresnum_b = vdms_df.ix[pair[1], 'resindex_vdm']
    return max(vdmresnum_b ,vdmresnum_a) - min(vdmresnum_b, vdmresnum_a)

def check_type(vdms_df, pair):
    polar_atoms = 'NO'
    vdw_atoms = 'CS'

    int_type = []
    # determine if both vdms make polar or vdw interactions
    vdmatoms_a = vdms_df.ix[pair[0], 'atom_names_vdm'].split(' ')
    vdmatoms_b = vdms_df.ix[pair[1], 'atom_names_vdm'].split(' ')
    
    polar = [False, False] # first element is vdmatoms_a, second is vdmatoms_b
    vdw = [False, False]
    for ind, group in enumerate([vdmatoms_a, vdmatoms_b]): #ind says atoms_a or atoms_b
        for atom in group:
            if atom[0] in polar_atoms:
                polar[ind] = True
            elif atom[0] in vdw_atoms:
                vdw[ind] = True
            else:
                print(atom, 'neither polar nor hydrophobic')
    if polar==[True, True]:
        int_type.append('polar')
    if vdw==[True, True]:
        int_type.append('vdw')
    return int_type
    
def bb_or_sc(vdms_df, pair, interaction_type):
    bb_atoms = ['N', 'CA', 'C', 'O', 'OXT']
    bb_sc_list = []

    if interaction_type != 'hbond':
        vdmatoms_a = vdms_df.ix[pair[0], 'atom_names_vdm'].split(' ')
        vdmatoms_b = vdms_df.ix[pair[1], 'atom_names_vdm'].split(' ')

    else:
        vdmatoms_a = parse_hbond_atoms(vdms_df.ix[pair[0], 'vdM_atom_names'])
        vdmatoms_b = parse_hbond_atoms(vdms_df.ix[pair[1], 'vdM_atom_names'])

    bb = [False, False] # first element is vdmatoms_a, second is vdmatoms_b
    sc = [False, False]

    for ind, group in enumerate([vdmatoms_a, vdmatoms_b]): #ind says atoms_a or atoms_b
        for atoms in group:
            for atom in atoms:
                if atom in bb_atoms:
                    bb[ind] = True
                else:
                    sc[ind] = True
    if bb==[True, True]:
        bb_sc_list.append('bb/bb')
    if sc==[True, True]:
        bb_sc_list.append('sc/sc')
    if (True in bb) and (True in sc):
        bb_sc_list.append('bb/sc')
    return bb_sc_list
    
def parse_hbond_atoms(atomstring):
    ''' Helper function to convert atoms from string 
    to list'''
    atoms_list = []
    atoms = atomstring.lstrip('(')
    atoms = atoms.rstrip(')')
    atoms = atoms.split(') (')
    for atom in atoms:
        atoms_list.append(atom.split(' '))
    return atoms_list
