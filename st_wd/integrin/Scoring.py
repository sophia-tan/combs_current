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

def get_interacting_atoms(parsed, resi, resn):    
    target = parsed.select('chain %s and resnum %s'%(resi[0], resi[1:]))
    assert len(list(set(target.getResnames()))) == 1
    assert constants.one_letter_code[target.getResnames()[0]] == resn

    # find out if interaction is BB or SC
    other_chain = parsed.select('not chain %s'%resi[0])
    bb = ['C', 'O', 'OXT', 'CA', 'N']
    polar = ['O', 'N']
    interacting_atoms = []
    for atom in target:
        radius = 3.5
        if atom.getName()[0] in polar:
            radius = 3.5
        else:
            radius = 4.8
        for nbr in pr.findNeighbors(atom, radius, other_chain):
            ifgatom, vdmatom, dist = nbr
            ifgindex, vdmindex = ifgatom.getIndex(), vdmatom.getIndex()
            ifgatomelem, vdmatomelem = ifgatom.getName(), vdmatom.getName()
            if dist <= 3.5:
                interacting_atoms.append((ifgindex,vdmindex))
            else:
                if ifgatomelem[0]=='C' and vdmatomelem[0]=='C':
                    interacting_atoms.append((ifgindex,vdmindex))
    return list(set(interacting_atoms))

def get_ifg_vdm(parsed,ifgresn, vdmresn, ifg_contact_atoms, vdm_contact_atoms):
    '''determine if the ifg is making bb or sc intrxns, and if the vdm is making 
    bb or sc intrxns'''

    ifgcontactatoms= [parsed.select('index %s'%x).getNames()[0] for x in ifg_contact_atoms]
    vdmcontactatoms= [parsed.select('index %s'%x).getNames()[0] for x in vdm_contact_atoms]
    
    #ifg = ifg_naming[list(constants.interactamer_atoms[targetresn].keys())[0]]
    
    bb = ['C', 'O', 'CA', 'N'] # excluding OXT bc don't want to deal with that 
    ifgtype = []
    vdmtype = []
    if len(set(bb).intersection(set(ifgcontactatoms))) > 0: # if there are bb atoms
        ifgtype.append('bb')
    if len(set(bb).intersection(set(ifgcontactatoms))) < len(set(ifgcontactatoms)): # sc atoms
        ifgtype.append('sc')
    if len(set(bb).intersection(set(vdmcontactatoms))) > 0: # if there are bb atoms
        vdmtype.append('bb')
    if len(set(bb).intersection(set(vdmcontactatoms))) < len(set(vdmcontactatoms)): # sc atoms
        vdmtype.append('sc')
    
    def assign(typ,res):
        if typ=='bb':
            return 'backbone'
        elif typ=='sc':
            return constants.AAname[res]
    
    ifgtype = [assign(typ,ifgresn) for typ in ifgtype]
    vdmtype = [assign(typ,vdmresn) for typ in vdmtype]

    ifginfo = [parsed.select('index %s'%ifg_contact_atoms[0]).getChids()[0], parsed.select('index %s'%ifg_contact_atoms[0]).getResnums()[0]]
    vdminfo = [parsed.select('index %s'%vdm_contact_atoms[0]).getChids()[0], parsed.select('index %s'%vdm_contact_atoms[0]).getResnums()[0]]
    return ifgtype, vdmtype, ifginfo, vdminfo

def getcoords_and_dump(parsed,ifgtype,vdmtype,ifginfo,vdminfo,lookup_dir, db_dir, pickled=True):
    ''' get coords of ifgs and vdms, and then dumps in pickle file'''
    rmsd_lists = [] # multiple lists bc it's for every ifg and every vdm in an intrxn
    for ifg in ifgtype:
        for vdm in vdmtype:
            if pickled==False:
                try:
                    x=pkl.load(open('./output_data/ifg_{}_vdm_{}_coords.pkl'.format(ifg,vdm),'rb'))
                except:
                    print(ifg,vdm, 'getting combed vdms')
                    get_coordslist(parsed,ifg, vdm, ifginfo, vdminfo, lookup_dir, db_dir)
            else: # already pickled, so can load
                pklf=pkl.load(open('./output_data/ifg_{}_vdm_{}_coords.pkl'.format(ifg,vdm),'rb'))
                ifg_rmsds, vdm_rmsds = calc_integrin_rmsds(parsed, ifg, vdm, ifginfo, \
                    vdminfo, lookup_dir, db_dir, pklf)
                rmsd_lists.append([ifg,vdm,ifg_rmsds,vdm_rmsds])
    return rmsd_lists

def calc_integrin_rmsds(parsed, ifg, vdm, ifginfo, vdminfo, lookup_dir, db_dir, pklf):
    '''for every bb/sc interaction, get coords for all ifgs and vdms'''
    # need to add each atom one by one bc if you ask prody to select a list of atoms
    #, the order might not be consistent, depending on how the authors deposited pdb
    query_ifg_atoms = [] # align on ifg atoms
    query_vdm_atoms = [] # get rmsd of vdm atoms to see how many vdms have same geom

    # get coordinates of ifg. superpose onto combed vdms. then, record rmsds of 'ifg' and vdm.
    # ifg in this case could be backbone, or residue sc
    ifg_num_atoms = len(constants.AA_sc_dict[constants.AAname_rev[ifg]])
    vdm_num_atoms = len(constants.AA_sc_dict[constants.AAname_rev[vdm]])

    for atom in constants.AA_sc_dict[constants.AAname_rev[ifg]]:
        selection = parsed.select('chain {} and resnum {} and name {}'.format(ifginfo[0],ifginfo[1],atom))
        query_ifg_atoms.append(selection.getCoords()[0])
    for atom in constants.AA_sc_dict[constants.AAname_rev[vdm]]:
        selection = parsed.select('chain {} and resnum {} and name {}'.format(vdminfo[0],vdminfo[1],atom))
        query_vdm_atoms.append(selection.getCoords()[0])
    query_total = query_ifg_atoms + query_vdm_atoms
    query_ifg_atoms = np.array(query_ifg_atoms)
    query_vdm_atoms = np.array(query_vdm_atoms)
    query_total = np.array(query_total)

    # get combed vdms from pklf
    lookup_ifgs, lookup_vdms, lookup_ifgandvdms, nocoords = pklf
    
    # superpose ifgs and get rmsds per atom
    ifg_rmsds = []
    vdm_rmsds = []
    
    assert len(lookup_ifgs)==len(lookup_vdms)==len(lookup_ifgandvdms)
    assert len(nocoords)/len(lookup_ifgs) < 0.01 # less than 1% of combed vdms have no coords info
    for i in range(len(lookup_ifgs)):
        movedifg,transf = pr.superpose(lookup_ifgs[i], query_ifg_atoms) # move lookup to query
        ifg_rmsds.append(pr.calcRMSD(movedifg,query_ifg_atoms)/ifg_num_atoms)
        movedvdm = pr.applyTransformation(transf,lookup_vdms[i]) # move lookup vdm to query's frame
        vdm_rmsds.append(pr.calcRMSD(movedvdm,query_vdm_atoms)/vdm_num_atoms)
        return None,None
    return ifg_rmsds, vdm_rmsds

def get_coordslist(parsed, ifgori, vdmori, ifginfo, vdminfo, lookup_dir, db_dir):

    if ifgori == 'backbone':
        ifg = 'glycine'
    else:
        ifg=ifgori
    if vdmori == 'backbone':
        vdm = 'glycine'
    else:
        vdm=vdmori
    # get coords for all combed vdms
    lookup = pkl.load(open(lookup_dir+'refinedvdms/vdms_of_{}.pkl'.format(ifg), 'rb'))
    lookup = lookup[lookup['resname_vdm']==constants.AAname_rev[vdm]]
    #lookup = lookup[:10] ###delete###
    nocoords = []
    lookup_ifgs = []
    lookup_vdms = []
    lookup_ifgandvdms = []

    for ix, row in lookup.iterrows():
        try:
            try:
                par = pr.parsePDB(db_dir+row['pdb']+'H.pdb')
            except:
                db_dir = '/home/gpu/Sophia/STcombs/20171118/database/reduce/'
                par = pr.parsePDB(db_dir+row['pdb']+'H.pdb')
            lookup_ifg_atoms = []
            lookup_vdm_atoms = []
            if 'segi_ifg' in lookup.columns:
                ifgchid, ifgresnum, ifgseg = row['chid_ifg'], row['resnum_ifg'], row['segi_ifg']
                vdmchid, vdmresnum, vdmseg = row['chid_vdm'], row['resnum_vdm'], row['segi_vdm']
                for atom in constants.AA_sc_dict[constants.AAname_rev[ifgori]]:
                    selection = par.select('chain {} and segment {} and resnum {} and name {}'.format(ifgchid, \
                        ifgseg, ifgresnum, atom))
                    lookup_ifg_atoms.append(selection.getCoords()[0])
                for atom in constants.AA_sc_dict[constants.AAname_rev[vdmori]]:
                    selection = par.select('chain {} and segment {} and resnum {} and name {}'.format(vdmchid, \
                        vdmseg, vdmresnum, atom))
                    lookup_vdm_atoms.append(selection.getCoords()[0])
            else: # old combing and no segment info
                ifgchid, ifgresnum = row['chid_ifg'], row['resnum_ifg']
                vdmchid, vdmresnum = row['chid_vdm'], row['resnum_vdm']
                for atom in constants.AA_sc_dict[constants.AAname_rev[ifgori]]:
                    selection = par.select('chain {} and resnum {} and name {}'.format(ifgchid, \
                        ifgresnum, atom))
                    lookup_ifg_atoms.append(selection.getCoords()[0])
                for atom in constants.AA_sc_dict[constants.AAname_rev[vdmori]]:
                    selection = par.select('chain {} and resnum {} and name {}'.format(vdmchid, \
                        vdmresnum, atom))
                    lookup_vdm_atoms.append(selection.getCoords()[0])
            lookup_total = lookup_ifg_atoms + lookup_vdm_atoms
            lookup_ifg_atoms = np.array(lookup_ifg_atoms)
            lookup_vdm_atoms = np.array(lookup_vdm_atoms)
            lookup_total = np.array(lookup_total)
            lookup_ifgs.append(lookup_ifg_atoms)
            lookup_vdms.append(lookup_vdm_atoms)
            lookup_ifgandvdms.append(lookup_total)
        except:
            nocoords.append('1')
    print('# unsucessfully processed combed vdms: ', len(nocoords))
    megalist = [lookup_ifgs, lookup_vdms, lookup_ifgandvdms, nocoords]
    pkl.dump(megalist, open('./output_data/ifg_{}_vdm_{}_coords.pkl'.format(ifgori,vdmori),'wb'))
