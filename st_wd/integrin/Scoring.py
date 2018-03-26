import sys, traceback, copy
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
#from ScoringFunctions import *
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
    if type(ifgresn)==list:
        ifgresn=ifgresn[1]
    if type(vdmresn)==list:
        vdmresn=vdmresn[1]

    ifgcontactatoms= [parsed.select('index %s'%x).getNames()[0] for x in ifg_contact_atoms]
    vdmcontactatoms= [parsed.select('index %s'%x).getNames()[0] for x in vdm_contact_atoms]
    
    #ifg = ifg_naming[list(constants.interactamer_atoms[targetresn].keys())[0]]
    
    bb = ['C', 'O', 'CA', 'N'] # excluding OXT bc don't want to deal with that 
    ifgtype = []
    vdmtype = []

    if len(set(constants.ifg_atoms[ifgresn]).intersection(set(ifgcontactatoms))) > 0: # if interacting through FG
        ifgtype.append('fg')
    if len(set(bb).intersection(set(ifgcontactatoms))) > 0: # if there are bb atoms
        ifgtype.append('bb')
    if len(set(bb).intersection(set(ifgcontactatoms))) < len(set(ifgcontactatoms)): # sc atoms
        ifgtype.append('sc')
    if len(set(constants.ifg_atoms[vdmresn]).intersection(set(vdmcontactatoms))) > 0: # if interacting through FG
        vdmtype.append('fg')
    if len(set(bb).intersection(set(vdmcontactatoms))) > 0: # if there are bb atoms
        vdmtype.append('bb')
    if len(set(bb).intersection(set(vdmcontactatoms))) < len(set(vdmcontactatoms)): # sc atoms
        vdmtype.append('sc')
    
    def assign(typ,res):
        if typ=='bb':
            return 'backbone'
        elif typ=='sc':
            return constants.AAname[res]
        elif typ=='fg':
            return [constants.AAname[res],constants.ifg_atoms[res]]
    
    ifgtype = [assign(typ,ifgresn) for typ in ifgtype]
    vdmtype = [assign(typ,vdmresn) for typ in vdmtype]
    ifginfo = [parsed.select('index %s'%ifg_contact_atoms[0]).getChids()[0], parsed.select('index %s'%ifg_contact_atoms[0]).getResnums()[0]]
    vdminfo = [parsed.select('index %s'%vdm_contact_atoms[0]).getChids()[0], parsed.select('index %s'%vdm_contact_atoms[0]).getResnums()[0]]
    return ifgtype, vdmtype, ifginfo, vdminfo

def getcoords_and_dump(parsed,ifgtype,vdmtype,ifginfo,vdminfo,lookup_dir, db_dir, pickled=True, ifgs=False):
    ''' get coords of ifgs and vdms, and then dumps in pickle file.
    Use ifg=True if you want to align just the iFG groups and not whole sc'''
    rmsd_lists = [] # multiple lists bc it's for every ifg and every vdm in an intrxn
    for origifg in ifgtype:
        for origvdm in vdmtype:
            ifgFG=False
            vdmFG=False
            ifgatoms=None
            vdmatoms=None
            ifg,vdm=origifg,origvdm
            ifgname,vdmname=ifg,vdm
            if type(ifg)==list:
                ifg,ifgatoms=origifg[0],origifg[1]
                ifgFG=True
                ifgname=ifg+'FG'
            if type(vdm)==list:
                vdm,vdmatoms=origvdm[0],origvdm[1]
                vdmFG=True
                vdmname=vdm+'FG'
            if pickled==False:
                try:
                    x=pkl.load(open('./output_data/ifg_{}_vdm_{}_coords.pkl'.format(ifg,vdm),'rb'))
                except:
                    print(ifg,vdm, 'getting combed vdms')
                    get_coordslist(parsed,ifg, vdm, ifginfo, vdminfo, lookup_dir, db_dir)
            else: # already pickled, so can load
                pklf=pkl.load(open('./output_data/ifg_{}_vdm_{}_coords.pkl'.format(ifg,vdm),'rb'))
                rmsds = calc_integrin_rmsds(parsed, ifg, vdm, ifginfo, \
                    vdminfo, lookup_dir, db_dir, pklf, ifgFG=ifgFG,vdmFG=vdmFG,ifgatoms=ifgatoms,vdmatoms=vdmatoms)
                #
                #ifg_num_atoms = len(constants.AA_sc_dict[constants.AAname_rev[ifg]])
                #vdm_num_atoms = len(constants.AA_sc_dict[constants.AAname_rev[vdm]])
                #num_atoms = ifg_num_atoms+vdm_num_atoms
                #
                #rmsd_lists.append([ifgname,vdmname,rmsds,num_atoms])
    return rmsd_lists

def calc_integrin_rmsds(parsed, ifg, vdm, ifginfo, vdminfo, lookup_dir, db_dir, pklf, ifgFG=False,vdmFG=False, \
    ifgatoms=None,vdmatoms=None):
    '''for every bb/sc interaction, get coords for all ifgs and vdms. include ifgatoms and vdmatoms if interaction
    is through FG'''
    # need to add each atom one by one bc if you ask prody to select a list of atoms
    #, the order might not be consistent, depending on how the authors deposited pdb
    query_atoms = []

    # ifg/vdms in this case could be backbone, or residue sc
    ifg_full_atoms = constants.AA_sc_dict[constants.AAname_rev[ifg]]
    vdm_full_atoms = constants.AA_sc_dict[constants.AAname_rev[vdm]]
    for atom in ifg_full_atoms:
        selection = parsed.select('chain {} and resnum {} and name {}'.format(ifginfo[0],ifginfo[1],atom))
        query_atoms.append(selection.getCoords()[0])
    for atom in vdm_full_atoms:
        selection = parsed.select('chain {} and resnum {} and name {}'.format(vdminfo[0],vdminfo[1],atom))
        query_atoms.append(selection.getCoords()[0])
    #query_atoms = np.array(query_atoms)

    # get combed vdms from pklf. we ony care about lookup_ifgandvdms
    lookup_ifgs, lookup_vdms, lookup_ifgandvdms, nocoords = pklf
    
    # superpose interactamers 
    rmsds = []
    # choose just the vdms that are interacting through its FG and count how many
    
    
    assert len(lookup_ifgs)==len(lookup_vdms)==len(lookup_ifgandvdms)
    assert len(nocoords)/len(lookup_ifgs) < 0.01 # less than 1% of combed vdms have no coords info
    #for ind in range(len(lookup_ifgandvdms)):
    #    result = get_rmsds(ifgFG,vdmFG,ifgatoms,vdmatoms,ifg_full_atoms,vdm_full_atoms,query_atoms,\
    #        lookup_ifgandvdms[ind])
    #    rmsds.append(result)
    return rmsds

def get_rmsds(ifgFG,vdmFG,ifgatoms,vdmatoms,ifg_full_atoms,vdm_full_atoms,query_atoms,\
            lookup_atoms):
    '''looks really messy :( but doing this because need to try this with flipped\
    residues that have ambiguous atoms, like TYR'''
    keep_indices = []
    # do ifg first
    equivifglist=[ifg_full_atoms]
    equivvdmlist=[vdm_full_atoms]
    for pair in constants.all_equiv_atoms:
        if set(pair).issubset(ifg_full_atoms):
            equivinds = [i for i,n in enumerate(ifg_full_atoms) if n in pair]
            a,b=equivinds
            newifg=ifg_full_atoms.copy()
            newifg[a],newifg[b]=newifg[b],newifg[a]

        if set(pair).issubset(vdm_full_atoms):
            equivinds = [i for i,n in enumerate(vdm_full_atoms) if n in pair]
            a,b=equivinds
            newvdm=vdm_full_atoms.copy()
            newvdm[a],newvdm[b]=newvdm[b],newvdm[a]
    try:
        equivifglist.append(newifg)
    except:
        pass
    try:
        equivvdmlist.append(newvdm)
    except:
        pass

    # see which of the flipped residues will give best rmsd
    compare_rmsds = []
    for ifg_full_atoms in equivifglist:
        for vdm_full_atoms in equivvdmlist:
            if ifgFG==True:
                keep_indices += [i for i,c in enumerate(ifg_full_atoms) if c in ifgatoms]
            elif ifgFG==False:
                keep_indices += [i for i,c in enumerate(ifg_full_atoms)]
            if vdmFG==True:
                keep_indices += [i+len(ifg_full_atoms) for i,c in enumerate(vdm_full_atoms) if c in vdmatoms]
            elif vdmFG==False:
                keep_indices += [i+len(ifg_full_atoms) for i,c in enumerate(vdm_full_atoms)]
            query_indices = sorted(keep_indices) # otherwise, order of query atoms will change 
            # too and the rmsd's b/n flipped residues won't change...
            queryatoms = np.array([query_atoms[i] for i in query_indices])
            lookupatoms = np.array([lookup_atoms[i] for i in keep_indices])
            moved,transf = pr.superpose(lookupatoms, queryatoms) 
            compare_rmsds.append(pr.calcRMSD(moved,queryatoms))
    return min(compare_rmsds)

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

def get_clusters(targetresi,parsed,ifgtype,vdmtype,ifginfo,vdminfo,lookup_dir, db_dir, ifgs=False):
    ''' get coords of ifgs and vdms, 
    Use ifg=True if you want to align just the iFG groups and not whole sc'''
    rmsd_lists = [] # multiple lists bc it's for every ifg and every vdm in an intrxn
    for origifg in ifgtype:
        for origvdm in vdmtype:
            ifgFG=False
            vdmFG=False
            ifgatoms=None
            vdmatoms=None
            ifg,vdm=origifg,origvdm
            ifgname,vdmname=ifg,vdm
            if type(ifg)==list:
                ifg,ifgatoms=origifg[0],origifg[1]
                ifgFG=True
                ifgname=ifg+'FG'
            if type(vdm)==list:
                vdm,vdmatoms=origvdm[0],origvdm[1]
                vdmFG=True
                vdmname=vdm+'FG'


            info = [ifgFG,vdmFG,ifgatoms,vdmatoms]

            get_coordslistforcluster(targetresi,parsed,ifg, vdm, ifginfo, vdminfo, lookup_dir, db_dir, info)
    return rmsd_lists

def getcoords(row,info):
    ifgFG,vdmFG,ifgatoms,vdmatoms,ifgori,vdmori,ifg_full_atoms,vdm_full_atoms,query_atoms,lookup,targetresi = info
    try:
        try:
            db_dir = '/home/gpu/Sophia/combs/st_wd/20180207_db_molprobity_biolassem/'
            par = pr.parsePDB(db_dir+row['pdb']+'H.pdb')
        except:
            db_dir = '/home/gpu/Sophia/STcombs/20171118/database/reduce/'
            par = pr.parsePDB(db_dir+row['pdb']+'H.pdb')
        lookup_ifg_atoms = []
        lookup_vdm_atoms = []
        if 'segi_ifg' in lookup.columns:
            ifgchid, ifgresnum, ifgseg = row['chid_ifg'], row['resnum_ifg'], row['segi_ifg']
            vdmchid, vdmresnum, vdmseg = row['chid_vdm'], row['resnum_vdm'], row['segi_vdm']
            printout = copy.deepcopy(par)
            printout = printout.select('(chain {} and segment {} and resnum {}) or (chain {} and segment {} and resnum {})'.format(ifgchid,ifgseg,ifgresnum,vdmchid,vdmseg,vdmresnum))
            printout.select('chain {} and segment {} and resnum {}'.format(ifgchid,ifgseg,ifgresnum)).setChids('Y')
            printout.select('chain {} and segment {} and resnum {}'.format(vdmchid,vdmseg,vdmresnum)).setChids('X')
            printout.select('all').setResnums(10)
            for atom in constants.AA_sc_dict[constants.AAname_rev[ifgori]]:
                selection = par.select('chain {} and segment {} and resnum {} and name {}'.format(ifgchid,ifgseg,ifgresnum,atom))
                lookup_ifg_atoms.append(selection.getCoords()[0])
            for atom in constants.AA_sc_dict[constants.AAname_rev[vdmori]]:
                selection = par.select('chain {} and segment {} and resnum {} and name {}'.format(vdmchid,vdmseg,vdmresnum,atom))
                lookup_vdm_atoms.append(selection.getCoords()[0])
        else: # old combing and no segment info
            ifgchid, ifgresnum = row['chid_ifg'], row['resnum_ifg']
            vdmchid, vdmresnum = row['chid_vdm'], row['resnum_vdm']
            printout = copy.deepcopy(par)
            printout = printout.select('(chain {} and resnum {}) or (chain {} and resnum {})'.format(ifgchid,ifgresnum,vdmchid,vdmresnum))
            printout.select('chain {} and resnum {}'.format(ifgchid,ifgresnum)).setChids('Y')
            printout.select('chain {} and resnum {}'.format(vdmchid,vdmresnum)).setChids('X')
            printout.select('all').setResnums(10)
            for atom in constants.AA_sc_dict[constants.AAname_rev[ifgori]]:
                selection = par.select('chain {} and resnum {} and name {}'.format(ifgchid,ifgresnum,atom))
                lookup_ifg_atoms.append(selection.getCoords()[0])
            for atom in constants.AA_sc_dict[constants.AAname_rev[vdmori]]:
                selection = par.select('chain {} and resnum {} and name {}'.format(vdmchid,vdmresnum,atom))
                lookup_vdm_atoms.append(selection.getCoords()[0])
        lookup_total = lookup_ifg_atoms + lookup_vdm_atoms
        #can print out selection, but check rmsd first to integrins
        rmsd = get_rmsds(ifgFG,vdmFG,ifgatoms,vdmatoms,ifg_full_atoms,vdm_full_atoms,query_atoms,\
            lookup_total)
        if rmsd <= 0.5:
            outdir = '/home/gpu/Sophia/combs/st_wd/integrin/output_data/pdbfiles/'
            pr.writePDB(outdir+'cutoff/{}_{}_{}_{}.pdb'.format(targetresi,ifgori,vdmori,row.name),printout)
            pr.writePDB(outdir+'{}_{}_{}_{}.pdb'.format(targetresi,ifgori,vdmori,row.name),printout)
            return [np.array(lookup_total),rmsd]
        else:
            outdir = '/home/gpu/Sophia/combs/st_wd/integrin/output_data/pdbfiles/'
            pr.writePDB(outdir+'{}_{}_{}_{}.pdb'.format(targetresi,ifgori,vdmori,row.name),printout)
            return [np.array(lookup_total),rmsd]
    except Exception:
        traceback.print_exc()
        return np.nan

def get_coordslistforcluster(targetresi,parsed, ifgori, vdmori, ifginfo, vdminfo, lookup_dir, db_dir,info):
    if ifgori == 'backbone':
        ifg = 'glycine'
    else:
        ifg=ifgori
    if vdmori == 'backbone':
        vdm = 'glycine'
    else:
        vdm=vdmori

    query_atoms = []

    # ifg/vdms in this case could be backbone, or residue sc
    ifg_full_atoms = constants.AA_sc_dict[constants.AAname_rev[ifg]]
    vdm_full_atoms = constants.AA_sc_dict[constants.AAname_rev[vdm]]
    for atom in ifg_full_atoms:
        selection = parsed.select('chain {} and resnum {} and name {}'.format(ifginfo[0],ifginfo[1],atom))
        query_atoms.append(selection.getCoords()[0])
    for atom in vdm_full_atoms:
        selection = parsed.select('chain {} and resnum {} and name {}'.format(vdminfo[0],vdminfo[1],atom))
        query_atoms.append(selection.getCoords()[0])
    #query_atoms = np.array(query_atoms)

    # get coords for all combed vdms
    lookup = pkl.load(open(lookup_dir+'refinedvdms/vdms_of_{}.pkl'.format(ifg), 'rb'))
    lookup = lookup[lookup['resname_vdm']==constants.AAname_rev[vdm]]
    #lookup = lookup[:10] ###delete###

    info += [ifgori,vdmori,ifg_full_atoms,vdm_full_atoms,query_atoms, lookup, targetresi]
    lookup['coords'] = lookup.apply(getcoords,info=info,axis=1)
    pkl.dump(lookup, open('./output_data/pdbfiles/ifg_{}_vdm_{}_coords.pkl'.format(ifgori,vdmori),'wb'))
