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
    other_chain = parsed.select('not (chain {} and resnum {})'.format(resi[0],resi[1:]))
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
    
    bb = ['C', 'O', 'OXT', 'CA', 'N']
    
    ifgfg = len(set(constants.ifg_atoms[ifgresn]).intersection(set(ifgcontactatoms)))
    ifgbb = len(set(bb).intersection(set(ifgcontactatoms))) 
    vdmbb = len(set(bb).intersection(set(vdmcontactatoms))) 
    vdmfg = len(set(constants.ifg_atoms[vdmresn]).intersection(set(vdmcontactatoms)))

    if ifgbb > ifgfg:
        ifgcontactatoms = ['N','CA','C','O']
    elif ifgfg>=ifgbb:
        ifgcontactatoms=constants.ifg_atoms[ifgresn]
    if vdmbb > vdmfg:
        vdmcontactatoms = ['N','CA','C','O']
    elif vdmfg>=vdmbb:
        vdmcontactatoms=constants.ifg_atoms[vdmresn]
    
    print(ifgcontactatoms,vdmcontactatoms)


    def assign(atomlist,res):
        return [constants.AAname[res],atomlist]
    
    ifgtype = [assign(ifgcontactatoms,ifgresn)]
    vdmtype = [assign(vdmcontactatoms,vdmresn)]
    ifginfo = [parsed.select('index %s'%ifg_contact_atoms[0]).getChids()[0], parsed.select('index %s'%ifg_contact_atoms[0]).getResnums()[0]]
    vdminfo = [parsed.select('index %s'%vdm_contact_atoms[0]).getChids()[0], parsed.select('index %s'%vdm_contact_atoms[0]).getResnums()[0]]
    return ifgtype, vdmtype, ifginfo, vdminfo

def get_clusters(targetresi,parsed,ifgtype,vdmtype,ifginfo,vdminfo,lookup_dir, db_dir, int_res):
    ''' get coords of ifgs and vdms, 
    Use ifg=True if you want to align just the iFG groups and not whole sc'''
    df_lists = []# multiple lists bc it's for every ifg and every vdm in an intrxn
    for origifg in ifgtype:
        for origvdm in vdmtype:
            ifgres,ifgatoms=origifg
            vdmres,vdmatoms=origvdm
            if ifgres!='cysteine':
                info = [ifgres,vdmres,ifgatoms,vdmatoms]
                df = get_coordslistforcluster(targetresi,parsed,ifgres, vdmres, ifgatoms,vdmatoms,lookup_dir,db_dir,ifginfo,vdminfo,info,int_res)
                df_lists.append(df)
    return df_lists

def getcoords(row,info):
    ifgres,vdmres,ifgatoms,vdmatoms,query_atoms,lookup,targetresi,ignore = info
    try:
        try:
            db_dir = '/home/gpu/Sophia/STcombs/20171118/database/reduce/'
            par = pr.parsePDB(db_dir+row['pdb']+'H.pdb')
        except:
            db_dir = '/home/gpu/Sophia/combs/st_wd/20180207_db_molprobity_biolassem/'
            par = pr.parsePDB(db_dir+row['pdb']+'H.pdb')
        lookup_ifg_atoms = []
        lookup_vdm_atoms = []
        
        ifgchid, ifgresnum = row['chid_ifg'], row['resnum_ifg']
        vdmchid, vdmresnum = row['chid_vdm'], row['resnum_vdm']
        printout = copy.deepcopy(par)
        printout = printout.select('(chain {} and resnum {}) or (chain {} and resnum {})'.format(ifgchid,ifgresnum,vdmchid,vdmresnum))
        printout.select('chain {} and resnum {}'.format(ifgchid,ifgresnum)).setChids('Y')
        printout.select('chain {} and resnum {}'.format(vdmchid,vdmresnum)).setChids('X')
        printout.select('all').setResnums(10)
        
        ifgatoms = list(set(ifgatoms))
        vdmatoms = list(set(vdmatoms))
        
        # account for flipped residues  
        # 1) do ifg
        lookupifg_list = []
        lookupifg_list.append(ifgatoms)
        if ifgres in constants.flip_names:
            copyls = ifgatoms.copy()
            for ix, atom in enumerate(ifgatoms):
                for pair in constants.flip_names[ifgres]:
                    if pair[0] not in ifgatoms and pair[1] not in ifgatoms:
                        pass
                    elif pair[0] in ifgatoms or pair[1] in ifgatoms:
                        if atom in pair:
                            if atom==pair[0]:
                                copyls[ix] = pair[1]
                            elif atom==pair[1]:
                                copyls[ix] = pair[0]
                    if copyls not in lookupifg_list:
                        lookupifg_list.append(copyls)
        # 2) do vdm
        lookupvdm_list = []
        lookupvdm_list.append(vdmatoms)
        if vdmres in constants.flip_names:
            copyls = vdmatoms.copy()
            for ix, atom in enumerate(vdmatoms):
                for pair in constants.flip_names[vdmres]:
                    if pair[0] not in vdmatoms and pair[1] not in vdmatoms:
                        pass
                    elif pair[0] in vdmatoms or pair[1] in vdmatoms:
                        if atom in pair:
                            if atom==pair[0]:
                                copyls[ix] = pair[1]
                            elif atom==pair[1]:
                                copyls[ix] = pair[0]
                    if copyls not in lookupvdm_list:
                        lookupvdm_list.append(copyls)
        
        #outdir = '/home/gpu/Sophia/combs/st_wd/integrin/output_data/pdbfiles/'
        #pr.writePDB(outdir+'{}_{}_{}_{}.pdb'.format(targetresi,ifgres,vdmres,row.name),printout)

        def get_lookup_relevant_atoms(par,ifgatoms,ifgchid,ifgresnum,vdmatoms,vdmchid,vdmresnum):
            coords = []
            for atom in ifgatoms:
                selection = par.select('chain {} and resnum {} and name {}'.format(ifgchid,ifgresnum,atom))
                coords.append(selection.getCoords()[0])
            for atom in vdmatoms:
                selection = par.select('chain {} and resnum {} and name {}'.format(vdmchid,vdmresnum,atom))
                coords.append(selection.getCoords()[0])
            return np.array(coords)
        
        compare_rmsds = []
        query_atoms = np.array(query_atoms)
        for ifgls in lookupifg_list:
            for vdmls in lookupvdm_list:
                lookupatoms = get_lookup_relevant_atoms(par,ifgls,ifgchid,ifgresnum,vdmls,vdmchid,vdmresnum)
                moved,transf = pr.superpose(lookupatoms, query_atoms) 
                compare_rmsds.append(pr.calcRMSD(moved,query_atoms))
        return [min(compare_rmsds),len(query_atoms)]

    except Exception:
        traceback.print_exc()
        return np.nan

def get_coordslistforcluster(targetresi,parsed,ifg, vdm, ifgatoms,vdmatoms,lookup_dir,db_dir,ifginfo,vdminfo,info,int_res):
    
    try:
        pl = pkl.load(open('./output_data/{}_{}_{}_ifg_{}_vdm_{}_coords.pkl'.format(targetresi,int_res[0],int_res[1],ifg,vdm),'rb'))
        if len(pl) < 100:
            raise Exception
        return []
    except:
        query_atoms = []
        ifgatoms = list(set(ifgatoms))
        vdmatoms = list(set(vdmatoms))
        for atom in ifgatoms:
            selection = parsed.select('chain {} and resnum {} and name {}'.format(ifginfo[0],ifginfo[1],atom))
            query_atoms.append(selection.getCoords()[0])
        for atom in vdmatoms:
            selection = parsed.select('chain {} and resnum {} and name {}'.format(vdminfo[0],vdminfo[1],atom))
            query_atoms.append(selection.getCoords()[0])

        # get coords for all combed vdms
        lookup = pkl.load(open(lookup_dir+'refinedvdms/vdms_of_{}.pkl'.format(ifg), 'rb'))
        lookup = lookup[lookup['resname_vdm']==constants.AAname_rev[vdm]]
        #lookup = lookup[:5] ###delete###
        info += [query_atoms, lookup, targetresi,int_res[0]]
        lookup['rmsds'] = lookup.apply(getcoords,info=info,axis=1)
        pkl.dump(lookup, open('./output_data/{}_{}_{}_ifg_{}_vdm_{}_coords.pkl'.format(targetresi,int_res[0],int_res[1],ifg,vdm),'wb'))
        return lookup

