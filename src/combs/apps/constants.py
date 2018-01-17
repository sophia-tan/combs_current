__all__ = ['one_letter_code', 'ifg_sele_dict', 'resname_dict', 'three_to_one', 'one_to_three']

import collections

one_letter_code = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                   'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                   'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                   'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
                   'MSE': 'm', 'ANY': '.', 'APX': '*'}

three_letter_code = {}
for k, v in one_letter_code.items():
    three_letter_code[v] = k

#interactamer atoms are listed such that the flippable atoms are at indices 1 and 2
# interactamer_atoms = {'HIS': ['NE2', 'CD2', 'CE1', 'ND1', 'CG'],
#                       'LYS': ['CE', 'CD', 'NZ'],
#                       'ASP': ['CG', 'OD1', 'OD2', 'CB'],
#                       'PHE': ['CZ', 'CE1', 'CE2', 'CD1', 'CD2', 'CG'],
#                       'ASN': ['CG', 'OD1', 'ND2', 'CB'],
#                       'GLN': ['CD', 'OE1', 'NE2', 'CG'],
#                       'ALA': ['C', 'CA', 'CB'],
#                       'ARG': ['CZ', 'NH2', 'NH1', 'NE'],
#                       'THR': ['OG1', 'CB', 'CG2'],
#                       'GLY': ['CA', 'N', 'C'],
#                       'TYR': ['CZ', 'CE1', 'CE2', 'CD1', 'CD2', 'CG', 'OH'],
#                       'LEU': ['CG', 'CD1', 'CD2', 'CB'],
#                       'VAL': ['CB', 'CG1', 'CG2'],
#                       'GLU': ['CG', 'OE1', 'OE2', 'CD', 'CB'],
#                       'PRO': ['CB', 'CG', 'CD'],
#                       'SER': ['CB', 'OG', 'CA'],
#                       'CYS': ['CB', 'SG', 'CA'],
#                       'MET': ['SD', 'CG', 'CE', 'CB'],
#                       'TRP': ['NE1', 'CD1', 'CE2', 'CB', 'CD2', 'CZ2', 'CE3', 'CZ3', 'CH2'],
#                       'ILE': ['CG1', 'CD1', 'CB', 'CG2'],
#                       'BB_CO': ['C', 'O', 'CA'],
#                       'BB_NH': ['CA', 'N', 'C']}

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


dict_ = {'LYS': ['CE', 'CD', 'NZ'],
          'ASP': ['CG', 'OD1', 'OD2'],
          'PHE': ['CZ', 'CE1', 'CE2'],
          'ASN': ['CG', 'OD1', 'ND2'],
          'GLN': ['CD', 'OE1', 'NE2'],
          'ALA': ['C', 'CA', 'CB'],
          'ARG': ['CZ', 'NH2', 'NH1'],
          'THR': ['OG1', 'CB', 'CG2'],
          'GLY': ['CA', 'N', 'C'],
          'TYR': ['CZ', 'CE1', 'OH'],
          'LEU': ['CG', 'CD1', 'CD2'],
          'VAL': ['CB', 'CG1', 'CG2'],
          'GLU': ['CD', 'OE1', 'OE2'],
          'PRO': ['CB', 'CG', 'CD'],
          'SER': ['CB', 'OG', 'CA'],
          'CYS': ['CB', 'SG', 'CA'],
          'MET': ['SD', 'CG', 'CE'],
          'TRP': ['NE1', 'CD1', 'CE2'],
          'ILE': ['CG1', 'CD1', 'CG2'],
          }


flip_names = {'PHE': [('CE1', 'CE2'), ('CD1', 'CD2')],
              'ASP': [('OD1', 'OD2')],
              'GLU': [('OE1', 'OE2')],
              'ARG': [('NH1', 'NH2')],
              'TYR': [('CE1', 'CE2'), ('CD1', 'CD2')],
              'VAL': [('CG1', 'CG2')],
              'LEU': [('CE1', 'CE2')],
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

ifg_sele_dict = {}
# formatting: if there are multiple ways the atoms can be ordered in the ligand, make a list of dicts (ex: lonepair_imidazole)
# But if there are >1 of the same ifg in the ligand, have the list of ifgatoms be INSIDE the dict (ex: bakboneCO)
ifg_sele_dict['carboxamide'] = {'ASN': 'CB CG ND2 OD1', 'GLN': 'CG CD NE2 OE1', 'APX': 'C10 C11 N3 O1'}
ifg_sele_dict['carboxylate'] = {'ASP': 'CB CG OD2 OD1', 'GLU': 'CG CD OE2 OE1'}
ifg_sele_dict['imidazole'] = {'HIS': 'CG CD2 NE2 CE1 ND1'}
ifg_sele_dict['lonepair_imidazole'] = [{'APX': 'C13 N1 N6 C10 C12'}, {'APX':'C12 C10 N6 N1 C13'}, {'APX': 'N1 C13 C12 C10 N6'}, {'APX':'C10 C12 C13 N1 N6'}]
ifg_sele_dict['indole'] = {'TRP': 'CG CD1 NE1 CE2 CZ2 CH2 CZ3 CE3 CD2'}
ifg_sele_dict['tyrCOH'] = {'TYR': 'CZ OH'}
#ifg_sele_dict['backboneCO'] = {'APX': ['C8 O3','C19 O2']}
ifg_sele_dict['backboneCO'] = {'ANY':'C O'}
ifg_sele_dict['serineCOH'] = {'SER':'CB OG'}
ifg_sele_dict['thrCOH'] = {'THR':'CB OG1'}
ifg_sele_dict['phenyl'] = {'PHE': 'CG CD1 CD2 CE1 CE2 CZ', 'TYR': 'CG CD1 CD2 CE1 CE2 CZ'}
ifg_sele_dict['guanidino'] = {'ARG': 'NE CZ NH1 NH2'}
ifg_sele_dict['amino'] = {'LYS': 'CE NZ'}


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

