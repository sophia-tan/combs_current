import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
from ScoringFunctions import *
#from PyrosettaScores import *
from residues_integrin import *
import prody as pr
import numpy as np, pickle as pkl, pandas as pd
from pprint import pprint
from itertools import *
script, lookup_dir = sys.argv

atoms = {}
atoms['carboxylate'] = 4
atoms['guanidino'] = 4
atoms['phenyl'] = 6
atoms['carboxamide'] = 4
atoms['backboneCO'] = 2
atoms['thrCOH'] = 2
atoms['tyrCOH'] = 2
atoms['amino'] = 2



ifg_naming = {}
ifg_naming['Carboxamide'] = 'carboxamide'
ifg_naming['Carboxylate'] = 'carboxylate'
ifg_naming['PhenolOH'] = 'tyrCOH'
ifg_naming['Alcohol'] = 'thrCOH'
ifg_naming['Guano'] = 'guanidino'
ifg_naming['Amino'] = 'amino'
ifg_naming['Aryl'] = 'phenyl'

def get_interacting_atoms(parsed, resi, resn):    
    target = parsed.select('chain %s and resnum %s'%(resi[0], resi[1:]))
    assert len(list(set(target.getResnames()))) == 1
    assert constants.one_letter_code[target.getResnames()[0]] == resn

    # find out if interaction is BB or SC
    other_chain = parsed.select('not chain %s'%resi[0])
    bb = ['C', 'O', 'OXT', 'CA', 'N']
    polar = ['O', 'N']
    interacting_atoms = []
    #print(resi, resn)
    for atom in target:
        radius = 3.5
        #if atom.getName()[0] in polar:
        #    radius = 3.5
        #else:
        #    radius = 4.8
        for nbr in pr.findNeighbors(atom, radius, other_chain):
            ifgatom, vdmatom, dist = nbr
            ifgindex, vdmindex = ifgatom.getIndex(), vdmatom.getIndex()
            ifgatomelem, vdmatomelem = ifgatom.getName(), vdmatom.getName()
            if dist <= 3.5:
                interacting_atoms.append((ifgindex,vdmindex))
            #else:
            #    if ifgatomelem[0]=='C' and vdmatomelem[0]=='C':
            #        interacting_atoms.append(vdmindex)
            #        print(ifgatom, vdmatom, vdmatom.getResname())
    return list(set(interacting_atoms))

def get_ifg_vdm(targetresi, targetresn, int_resn, int_res_atoms, target_atoms):
    '''logic for selecting where the ifg is, and which residue it came from.
    this logic only works for this selection! not strict enough for other things
    look at target residue/atoms and the interacting residue/atoms it interacts with 
    (target = mutated res)
    if Y594, vdm is phenyl and ifg is the int_res's iFG
    if target is making backboneCO interaction, take the interacting residue's iFG
    else take target's ifg and see if the int_res_atoms are bbCO
    '''
    intresatoms = [parsed.select('index %s'%x).getNames()[0] for x in int_res_atoms]
    targetatoms = [parsed.select('index %s'%x).getNames()[0] for x in target_atoms]
    if targetresi == 'B594':
        vdm = 'phenyl'
        ifg = ifg_naming[list(constants.interactamer_atoms[int_res].keys())[0]]
        vdminfo = ['B', 594]
        ifginfo = ['A', parsed.select('index %s'%int_res_atoms[0]).getResnums()[0]]
        ifgres = int_res
        vdmres = targetresn
    elif 'O' in targetatoms:
        ifg = ifg_naming[list(constants.interactamer_atoms[int_res].keys())[0]]
        vdm = 'backboneCO'
        ifgres = int_res
        vdmres = targetresn
        ifginfo = [parsed.select('index %s'%int_res_atoms[0]).getChids()[0], parsed.select('index %s'%int_res_atoms[0]).getResnums()[0]]
        vdminfo = [parsed.select('index %s'%target_atoms[0]).getChids()[0], parsed.select('index %s'%target_atoms[0]).getResnums()[0]]
    else:
        ifg = ifg_naming[list(constants.interactamer_atoms[targetresn].keys())[0]]
        ifginfo = [parsed.select('index %s'%target_atoms[0]).getChids()[0], parsed.select('index %s'%target_atoms[0]).getResnums()[0]]
        ifgres = targetresn
        if 'O' in intresatoms:
            vdm = 'backboneCO'
        elif 'CD1' in intresatoms or 'CE1' in intresatoms:
            vdm = 'phenyl'
        else:
            vdm = ifg_naming[list(constants.interactamer_atoms[int_res].keys())[0]]
        vdminfo = [parsed.select('index %s'%int_res_atoms[0]).getChids()[0], parsed.select('index %s'%int_res_atoms[0]).getResnums()[0]]
        vdmres = int_res
    return ifg, vdm, ifginfo, vdminfo, ifgres, vdmres

def get_rmsdlist(ifg, vdm, ifginfo, vdminfo, ifgres, vdmres):
    rmsdlist = []
    normsd = []
    lookup = pkl.load(open(lookup_dir+'parsedvdms/parsedvdms_%s_%s.pkl'%(ifg,vdm),'rb'))
    query_atoms = [] # need to add each atom one by one bc if you ask prody to select a list of atoms, the order 
    # might not be consistent, depending on how the authors deposited pdb
    for typ in [[ifg, ifgres, ifginfo], [vdm,vdmres,vdminfo]]:
        if typ[0] != 'backboneCO':
            for atom in constants.ifg_sele_dict[typ[0]][typ[1]].split(' '):
                selection = parsed.select('chain %s and resnum %s and name %s'%(typ[2][0], str(typ[2][1]), atom))
                query_atoms.append(selection.getCoords()[0])
        else:
            for atom in ['C', 'O']:
                selection = parsed.select('chain %s and resnum %s and name %s'%(typ[2][0], str(typ[2][1]), atom))
                query_atoms.append(selection.getCoords()[0])

    query_atoms = np.array(query_atoms)
    # get rmsd for every parsedvdm
    for parsed_vdm in lookup:
        lookup_atoms = []
        parsedifg = parsed_vdm.select('chain Y and resnum 10').getResnames()[0]
        parsedvdm = parsed_vdm.select('chain X and resnum 10').getResnames()[0]
        
        try:
            for zipped in [[ifg,'Y',parsedifg], [vdm,'X',parsedvdm]]:
                if zipped[0] != 'backboneCO':
                    for atom in constants.ifg_sele_dict[zipped[0]][zipped[2]].split(' '):
                        selection = parsed_vdm.select('chain %s and resnum 10 and name %s'%(zipped[1], atom))
                        lookup_atoms.append(selection.getCoords()[0])
                else:
                    for atom in ['C', 'O']:
                        selection = parsed_vdm.select('chain %s and resnum 10 and name %s'%(zipped[1], atom))
                        lookup_atoms.append(selection.getCoords()[0])
            lookup_atoms = np.array(lookup_atoms)
            moved = pr.superpose(lookup_atoms, query_atoms)[0]
            rmsd = pr.calcRMSD(moved, query_atoms)
            rmsdlist.append(rmsd)
        except:
            normsd.append('1')
    return rmsdlist, normsd


#### START CODE ####
bb = ['C', 'O', 'OXT', 'CA', 'N']
parsed = pr.parsePDB('integrin.pdb')
geomdict = {}
for resi, resn in integrin_res.items():
    # for each target res, find the residues w/in 3.5 and 4.8A
    interacting_atoms = get_interacting_atoms(parsed, resi, resn)
    # list of list where inner list is [targetindex, int_resindex]
    # can get unique interacting residues by resname bc no repeat interacting resname/resnum per target res
    interacting_res = set(list([parsed.select('index %s'%x[1]).getResnames()[0] for x in interacting_atoms]))
    if len(interacting_res) > 0:
        geomdict[(resi,resn)] = {}
        print('--------------------')
        print('Target Res', resi, resn)
        for int_res in interacting_res: 
            # interacting_atoms has all the indices for all the atoms the target resn interacts with, but we only
            # want to look at the ones for this int_res
            int_res_atoms = [x[1] for x in interacting_atoms if parsed.select('index %s'%x[1]).getResnames()[0] == int_res]
            target_res_atoms = [x[0] for x in interacting_atoms if parsed.select('index %s'%x[1]).getResnames()[0] == int_res]
            print(int_res)
            ifg, vdm, ifginfo, vdminfo, ifgres, vdmres = get_ifg_vdm(resi, constants.three_letter_code[resn], int_res, int_res_atoms, target_res_atoms)
            # get num of atoms rmsd is taken over
            num_atoms = atoms[ifg] + atoms[vdm]

            rmsdlist,normsd = get_rmsdlist(ifg, vdm, ifginfo, vdminfo, ifgres, vdmres)
            rmsdlist = [x/num_atoms for x in rmsdlist]
            print(len(rmsdlist), len(normsd))
            geomdict[(resi,resn)][int_res] = [rmsdlist, normsd]

pkl.dump(geomdict, open('newgeomdict.pkl','wb'))
