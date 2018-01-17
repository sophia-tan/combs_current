import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
import combs

# stores input dicts for Combs objects
all_ifgs = {}
all_ifgs['carboxamide'] = {'ASN': 'CB CG OD1 ND2', 'GLN': 'CG CD OE1 NE2'}
all_ifgs['carboxylate'] = {'ASP': 'CB CG OD1 OD2', 'GLU': 'CG CD OE1 OE2'}
all_ifgs['amino'] = {'LYS': 'CE CD NZ'}
all_ifgs['imidazole'] = {'HIS': 'CG CD2 NE2 CE1 ND1'}
all_ifgs['indole'] = {'TRP': 'CE2 NE1 CD1'}
all_ifgs['tyrCOH'] = {'TYR': 'CZ OH'}
all_ifgs['backboneCO'] = {'ANY': 'C O'}
all_ifgs['serCOH'] = {'SER':'CB OG'}
all_ifgs['thrCOH'] = {'THR':'CB OG1'}


