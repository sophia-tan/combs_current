__all__ = ['one_letter_code', 'ifg_sele_dict', 'resname_dict', 'three_to_one', 'one_to_three', 'AAname','AA_sc_dict','AAname_rev','atoms_dict']

import collections

one_letter_code = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                   'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                   'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                   'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
                   'MSE': 'm', 'ANY': '.', 'APX': '*', 'GNR':'*', 'bb': '?'}


AAname = {'CYS': 'cysteine', 'ASP': 'aspartate', 'SER': 'serine', 'GLN': 'glutamine', 'LYS': 'lysine', 
            'ILE': 'isoleucine', 'PRO': 'proline', 'THR': 'threonine', 'PHE': 'phenylalanine', 
            'ASN': 'asparagine', 'GLY': 'glycine', 'HIS': 'histidine', 'LEU': 'leucine', 'ARG': 'arginine', 
            'TRP': 'tryptophan', 'ALA': 'alanine', 'VAL': 'valine', 'GLU': 'glutamate', 'TYR': 'tyrosine', 
            'MET': 'methionine', 'bb': 'backbone'}

AAname_rev = {}
for k,v in AAname.items():
    AAname_rev[v] = k

AA_sc_dict = {'ALA': ['CB'], 'CYS': ['CB', 'SG'], 'ASP': ['CB', 'CG', 'OD1', 'OD2'],
                    'ASN': ['CB', 'CG', 'OD1', 'ND2'], 'VAL': ['CB', 'CG1', 'CG2'],
                    'GLU': ['CB', 'CG', 'CD', 'OE1', 'OE2'], 'LEU': ['CB', 'CG', 'CD1', 'CD2'],
                    'HIS': ['CB', 'CG', 'ND1', 'CE1', 'NE2', 'CD2'],
                    'ILE': ['CB', 'CG2', 'CG1', 'CD1'], 'MET': ['CB', 'CG', 'SD', 'CE'],
                    'TRP': ['CB', 'CG', 'CD1', 'NE1', 'CE2', 'CD2', 'CE3', 'CZ3', 'CH2', 'CZ2'],
                    'SER': ['CB', 'OG'], 'LYS': ['CB', 'CG', 'CD', 'CE', 'NZ'],
                    'PHE': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'], 'PRO': ['CB', 'CG', 'CD'],
                    'GLY': [], 'THR': ['CB', 'OG1', 'CG2'],
                    'TYR': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
                    'GLN': ['CB', 'CG', 'CD', 'OE1', 'NE2'],
                    'ARG': ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
                    'bb': ['N', 'CA', 'C', 'O']}


atoms_dict = {}
for k,v in AA_sc_dict.items():
    atoms_dict[k] = ['N', 'CA', 'C', 'O'] + v



three_letter_code = {}
for k, v in one_letter_code.items():
    three_letter_code[v] = k

interactamer_atoms = collections.defaultdict(dict)
# format of origin, plane1, plane2, so make sure flipped atoms 
# are in elements 1 and 2
interactamer_atoms['HIS']['Delta'] = ['ND1', 'CG', 'CB']
interactamer_atoms['HIS']['Epsilon'] = ['NE2', 'CD2', 'CE1']
interactamer_atoms['LYS']['Amino'] = ['CE', 'CD', 'NZ']
interactamer_atoms['ASP']['Carboxylate'] = ['CG', 'OD1', 'OD2']
interactamer_atoms['GLN']['Carboxamide'] = ['CD', 'OE1', 'NE2']
interactamer_atoms['GLU']['Carboxylate'] = ['CD', 'OE1', 'OE2']
interactamer_atoms['ASN']['Carboxamide'] = ['CG', 'OD1', 'ND2']
interactamer_atoms['ALA']['Methyl'] = ['C', 'CA', 'CB']
interactamer_atoms['ARG']['Guano'] = ['CZ', 'NH2', 'NH1']
interactamer_atoms['THR']['Alcohol'] = ['OG1', 'CB', 'CG2']
interactamer_atoms['GLY']['MainChain'] = ['CA', 'N', 'C']
interactamer_atoms['TYR']['PhenolOH'] = ['CZ', 'CE1', 'OH']
interactamer_atoms['SER']['Alcohol'] = ['CB', 'OG', 'CA']
interactamer_atoms['TRP']['IndoleNH'] = ['NE1', 'CD1', 'CE2']
interactamer_atoms['CYS']['Thiol'] = ['CA', 'CB', 'SG']
interactamer_atoms['VAL']['Isopropyl'] = ['CA', 'CG1', 'CG2']
interactamer_atoms['LEU']['Isopropyl'] = ['CG', 'CD1', 'CD2']
interactamer_atoms['ILE']['Propyl'] = ['CB', 'CG1', 'CD1']
interactamer_atoms['MET']['Thioether'] = ['CG', 'SD', 'CE']
interactamer_atoms['PHE']['Aryl'] = ['CZ','CE1', 'CE2']
interactamer_atoms['PRO']['Pyrrole'] = ['CB', 'CG', 'CD']

flip_names = {'phenylalanine': [('CE1', 'CE2'), ('CD1', 'CD2')],
              'aspartate': [('OD1', 'OD2')],
              'glutamate': [('OE1', 'OE2')],
              'arginine': [('NH1', 'NH2')],
              'tyrosine': [('CE1', 'CE2'), ('CD1', 'CD2')],
              'valine': [('CG1', 'CG2')],
              'leucine': [('CE1', 'CE2')],
              }


flip_residues = ['PHE', 'ASP', 'GLU', 'ARG', 'VAL', 'LEU'] # TYR is special :(

flip_sets = [{'OD1', 'OD2'}, {'CE1', 'CE2'}, {'NH2', 'NH1'}, {'OE1', 'OE2'}, {'CG1', 'CG2'}]

bb_type_dict = {'N_CA': ['N', 'H', 'CA'], 'C_O': ['C', 'O', 'CA'], 'SC': ['CA', 'N', 'C'],
                'PHI_PSI': ['CA', 'N', 'C']}

residue_sc_names = {'ALA': ['CB'], 'CYS': ['CB', 'SG'], 'ASP': ['CB', 'CG', 'OD1', 'OD2'],
                    'ASN': ['CB', 'CG', 'OD1', 'ND2'], 'VAL': ['CB', 'CG1', 'CG2'],
                    'GLU': ['CB', 'CG', 'CD', 'OE1', 'OE2'], 'LEU': ['CB', 'CG', 'CD1', 'CD2'],
                    'HIS': ['CB', 'CG', 'ND1', 'CE1', 'NE2', 'CD2'],
                    'ILE': ['CB', 'CG2', 'CG1', 'CD1'], 'MET': ['CB', 'CG', 'SD', 'CE'],
                    'TRP': ['CB', 'CG', 'CD1', 'NE1', 'CE2', 'CD2', 'CE3', 'CZ3', 'CH2', 'CZ2'],
                    'SER': ['CB', 'OG'], 'LYS': ['CB', 'CG', 'CD', 'CE', 'NZ'],
                    'PHE': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'], 'PRO': ['CB', 'CG', 'CD'],
                    'GLY': [], 'THR': ['CB', 'OG1', 'CG2'],
                    'TYR': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
                    'GLN': ['CB', 'CG', 'CD', 'OE1', 'NE2'],
                    'ARG': ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2']}

equiv_atoms = {}
equiv_atoms['TYR'] = [['CD1','CD2'],['CE1','CE2']]
equiv_atoms['PHE'] = [['CD1','CD2'],['CE1','CE2']]
equiv_atoms['ARG'] = ['NH1', 'NH2']
equiv_atoms['ASP'] = ['OD1', 'OD2']
equiv_atoms['GLU'] = ['OE1', 'OE2']
equiv_atoms['VAL'] = ['CG1', 'CG2']

all_equiv_atoms = [['CD1', 'CD2'], ['NH1', 'NH2'], ['OD1', 'OD2]'], ['OE1', 'OE2]'], ['CG1', 'CG2'], ['CE1','CE2']]



planar_atoms = {}
planar_atoms['ASN'] = ['CB CG ND2 OD1']
planar_atoms['GLN'] = ['CG CD NE2 OE1']
#planar_atoms['TYR'] = ['CE1 CE2 CZ OH']
planar_atoms['TYR'] = ['CG CD1 CD2 CE1 CE2 CZ OH']
planar_atoms['ARG'] = ['NE CZ NH1 NH2','CG CD NE']
planar_atoms['ASP'] = ['CB CG OD2 OD1']
planar_atoms['ILE'] = ['CG2 CB CG1','CB CG1 CD1']
planar_atoms['GLU'] = ['CG CD OE2 OE1']
planar_atoms['THR'] = ['CG2 CB OG1']
planar_atoms['PHE'] = ['CG CD1 CD2 CE1 CE2 CZ']
planar_atoms['VAL'] = ['CG1 CB CG2']
planar_atoms['GLY'] = ['N CA C','CA C O']
planar_atoms['ALA'] = ['N CA CB','CA CB C']
planar_atoms['CYS'] = ['CA CB SG']
planar_atoms['HIS'] = ['CG ND1 CE1 NE2 CD2']
planar_atoms['LYS'] = ['CD CE NZ']
planar_atoms['MET'] = ['CG SD CE']
planar_atoms['PRO'] = ['CB CG CD']
planar_atoms['SER'] = ['CA CB OG']
planar_atoms['TRP'] = ['CG CD1 NE1 CE2 CD2', 'CE2 CD2 CE3 CZ3 CH2 CZ2']
planar_atoms['LEU'] = ['CD1 CG CD2']
#for k,val in planar_atoms.items():
#    planar_atoms[k] = [v.split(' ') for v in val]
#print(planar_atoms)



ifg_atoms = {}
ifg_atoms['ASN'] = 'CB CG ND2 OD1'
ifg_atoms['GLN'] = 'CG CD NE2 OE1'
ifg_atoms['TYR'] = 'CG CD1 CD2 CE1 CE2 CZ OH'
ifg_atoms['ARG'] = 'NE CZ NH1 NH2'
ifg_atoms['ASP'] = 'CB CG OD2 OD1'
ifg_atoms['ILE'] = 'CD1 CG1 CB CG2'
ifg_atoms['GLU'] = 'CG CD OE2 OE1'
ifg_atoms['THR'] = 'CG2 CB OG1'
ifg_atoms['PHE'] = 'CG CD1 CD2 CE1 CE2 CZ'
ifg_atoms['VAL'] = 'CG1 CB CG2'
ifg_atoms['GLY'] = 'N CA C O'
ifg_atoms['ALA'] = 'N CA CB C'
ifg_atoms['CYS'] = 'CB SG'
ifg_atoms['HIS'] = 'CG ND1 CE1 NE2 CD2'
ifg_atoms['LYS'] = 'CG CD CE NZ'
ifg_atoms['MET'] = 'CB CG SD CE'
ifg_atoms['PRO'] = 'CB CG CD'
ifg_atoms['SER'] = 'CA CB OG'
ifg_atoms['TRP'] = ''
ifg_atoms['LEU'] = 'CB CG CD1 CD2'
for k,v in ifg_atoms.items():
    ifg_atoms[k] = v.split(' ')


ifg_sele_dict = {}
# formatting: if there are multiple ways the atoms can be ordered in the ligand, make a list of dicts (ex: lonepair_imidazole)
# But if there are >1 of the same ifg in the ligand, have the list of ifgatoms be INSIDE the dict (ex: bakboneCO)
ifg_sele_dict['carboxamide'] = {'ASN': 'CB CG ND2 OD1', 'GLN': 'CG CD NE2 OE1', 'APX': 'C10 C11 N3 O1', 'GNR': 'CG CD NE2 OE1'}
ifg_sele_dict['carboxylate'] = {'ASP': 'CB CG OD2 OD1', 'GLU': 'CG CD OE2 OE1', 'GNR': 'CA C OXT O'}
ifg_sele_dict['imidazole'] = {'HIS': 'CG CD2 NE2 CE1 ND1'}
ifg_sele_dict['lonepair_imidazole'] = [{'HIS': 'CG CD2 NE2 CE1 ND1'}, {'APX': 'C13 N1 N6 C10 C12'}, {'APX':'C12 C10 N6 N1 C13'}, {'APX': 'N1 C13 C12 C10 N6'}, {'APX':'C10 C12 C13 N1 N6'}]
ifg_sele_dict['indole'] = {'TRP': 'CG CD1 NE1 CE2 CZ2 CH2 CZ3 CE3 CD2'}
ifg_sele_dict['hydroxyphenyl'] = {'TYR': 'CZ OH'}
#ifg_sele_dict['backboneCO'] = {'APX': ['C8 O3','C19 O2']}
ifg_sele_dict['backboneCO'] = {'ANY':'C O'}
ifg_sele_dict['hydroxyl'] = {'SER':'CB OG'}
ifg_sele_dict['hydroxymethyl'] = {'THR':'CG2 CB OG1'}
ifg_sele_dict['phenyl'] = {'PHE': 'CG CD1 CD2 CE1 CE2 CZ'}
ifg_sele_dict['guanidino'] = {'ARG': 'NE CZ NH1 NH2'}
ifg_sele_dict['amino'] = {'LYS': 'CE NZ', 'GNR': 'CA N'}
ifg_sele_dict['methyl'] = {'ALA': 'CA CB'}
ifg_sele_dict['isopropyl'] = {'VAL': 'CG1 CB CG2', 'LEU': 'CD1, CG, CD2'}
ifg_sele_dict['propyl'] = {'ILE': 'CD1 CG1 CB'}
ifg_sele_dict['thioether'] = {'MET': 'CG SD CE'}


resname_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                   'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                   'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                   'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
                   }

# dictionaries to convert 1-letter AA to 3-letter and vice versa

three_to_one = {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', \
    'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L', 'MET':'M', 'ASN':'N', \
    'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'}

one_to_three = {}
for key,val in three_to_one.items():
    one_to_three[val] = key

