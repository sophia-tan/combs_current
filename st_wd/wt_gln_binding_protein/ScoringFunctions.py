# gets cb sasa scores for each vdm of a FG for a design
# commandline arguments: pdb of design, txtfile of vdm_resnum vdm_aa 
# make sure there's a SINGLE SPACE b/n resnum and aa 
# and aa is 3 letters

import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
import prody as pr
import freesasa, math
import numpy as np, pickle as pkl, pandas as pd
from combs import analysis
from combs import parse
from combs import apps

#script, design_pdb, fg_vdm_txt, ligand_text= sys.argv

############################################################################################
################### make pandas df for vdms and their scores ###############################
############################################################################################
def make_df(fg_vdm_txt):
    vdms_aa = []
    vdms_resnum = []
    ifg = []
    with open(fg_vdm_txt) as inF:
        for line in inF:
            line = line.strip()
            line = line.split(' ')
            vdms_resnum.append(line[0])
            vdms_aa.append(line[1])
            ifg.append(line[2])
    columns = ['iFG', 'AA', 'BB or SC', 'sasa_cb_3a', 'burial score', 
                'dunbrack', 'rama', 'rosetta hbonds', 'f(AAi) bb', 'f(AAi) sc']
    df = pd.DataFrame(columns=columns)
    #df = pd.DataFrame(columns=columns, index=pd.Series(vdms_resnum, name='resnum'))
    # fill in df
    alphabet='abcdefg'
    for ix in range(len(vdms_aa)):
        if len(ifg[ix].split(','))>1:
            alphabet_ix = 0
            for fg in ifg[ix].split(','):
                index = str(vdms_resnum[ix])+alphabet[alphabet_ix]
                df = fill_df(df, vdms_resnum, index, ix, fg, vdms_aa)
                alphabet_ix += 1
        else: 
            index = str(vdms_resnum[ix])
            df = fill_df(df, vdms_resnum, index, ix, ifg[ix], vdms_aa)
    return df

def fill_df(df,vdms_resnum, index, ix, fg, vdms_aa):
    df.loc[index, 'iFG'] = fg
    df.loc[index, 'AA'] = vdms_aa[ix]
    return df
############################################################################################
########################### SASAs at cb 3A probe ###########################################
############################################################################################

# 1) get dictionary with vdm AAs, resnums, and sasas at cb 3A probe ################
def cb_sasas(design_pdb, score_df):
    sasa_dict = {} # key is resnum, value is list [aa, cbsasa]
    
    # parse design and do freesasa calc
    prody_parsed = pr.parsePDB(design_pdb, altloc='A', model=1)
    fs_struct = freesasa.Structure(design_pdb) # more atoms 
    fs_result = freesasa_cb(prody_parsed, probe_radius=3) # less atoms bc this is Cb cutoff 

    # get sasa
    for resnum, aa in zip(score_df.index, score_df['AA']): # shouldn't have to worry about 
            #neg resnums in designed proteins
        # get Cb atoms 
        try: 
            resnum = int(resnum)
        except:
            resnum = int(resnum[:-1])
        prody_pdb_bb_cb_atom_ind = prody_parsed.select('protein and (backbone or name CB) and \
            not element H D').getIndices()
        # get Cb atoms for resnum
        sele = prody_parsed.select('protein and (backbone or name CB) and resnum ' + str(resnum) \
            + ' and not element H D')
        bb_cb_atom_ind = sele.getIndices() 
        sasa_3A_probe = '{0:.2f}'.format(sum(fs_result.atomArea(i) for i in \
            np.where(np.in1d(prody_pdb_bb_cb_atom_ind,bb_cb_atom_ind))[0]))
        sasa_dict[int(resnum)] = [aa, float(sasa_3A_probe)]
    return sasa_dict

def freesasa_cb(prody_parsed, probe_radius=1.4):
    cb_sele = prody_parsed.select('protein and (backbone or name CB) and not element H D')
    coords = list(x for y in cb_sele.getCoords() for x in y)
    radii = list(freesasa.Classifier().radius(x, y) for x, y in zip(cb_sele.getResnames(), \
        cb_sele.getNames()))
    return freesasa.calcCoord(coords, radii, freesasa.Parameters({'probe-radius': probe_radius}))

# 2) 'find sasa bin for that vdm and score' ###########################

def scoring_sasa(sasadict, score_df, lookup_dir):
    # load lookup table
    lookup = pkl.load(open(lookup_dir+'scores_largeprobe_sasa_db_lookup.pkl','rb'))
    # take -np.log10 for all values in df, and drop rows 150-180
    lookup = lookup.drop([160, 170, 180, 190])
    lookup = lookup.applymap(lambda g: -np.log10(g))
    for ix, row in score_df.iterrows():
        try:
            ix = int(ix)
        except:
            ix = int(ix[:-1])
        value = sasadict[int(ix)]
        aa, sasa = value[0], value[1]
        if sasa > 150: # this site is exposed
            score = -1
        else:
            if sasa == 0:
                sasabin = 10
            else:
                sasabin = int(math.ceil(sasa / 10.0)) * 10 
            score = np.round(lookup.loc[sasabin,aa],2)
        row['sasa_cb_3a'] = sasa
        row['burial score'] = score

    return score_df

############################################################################################
############################ get freq_aai score ############################################
############################################################################################i
def get_ifgatoms(ligand_text):
    print(ligand_text)
    ifgdict = {}
    with open(ligand_text) as inF:
        for line in inF:
            line = line.strip()
            line = line.split(' ')
            ifg, resnum, atoms= line[0], line[1], line[2]
            atoms = atoms.split('+')
            atoms = ' '.join(atoms)
            ifgdict[ifg] = [resnum, atoms]
    return ifgdict

def freqaai(score_df, parsed_des, ligand_text, lookup_dir):
    # get ifgatoms from ifg_text
    ifgdict = get_ifgatoms(ligand_text)

    for ix, row in score_df.iterrows():
        vdm_name, res_num, ifg_name = row['AA'], ix, row['iFG']
        try:
            res_num = int(res_num)
        except:
            res_num = int(res_num[:-1])
        v = parsed_des.select('resnum '+str(res_num))
        v = v.select('not element H D')
        if v.getResnames()[0] != vdm_name:
            raise Exception
        
        # find out if interaction is BB or SC
        bb = ['C', 'O', 'OXT', 'CA', 'N']
        polar = ['O', 'N']
        ifgresnum, ifgatoms = ifgdict[ifg_name]
        ifg = parsed_des.select('resnum %s and name %s'%(ifgresnum, ifgatoms))

        vdmatoms = []
        for atom in ifg:
            if atom.getName()[0] in polar:
                radius = 3.5
            else:
                radius = 4.8
            for nbr in pr.findNeighbors(atom, radius, v):
                ifgatom, vdmatom, dist = nbr
                ifgatom, vdmatom = ifgatom.getName(), vdmatom.getName()
                if dist <= 3.5:
                    vdmatoms.append(vdmatom)
                else:
                    if ifgatom[0]=='C' and vdmatom[0]=='C':
                        vdmatoms.append(vdmatom)
        vdmatoms = list(set(vdmatoms))
        bbinteraction = sum([x in bb for x in vdmatoms])
        scinteraction = sum([x not in bb for x in vdmatoms])

        # get lookup info for ifg
        pklfile = lookup_dir + 'AAi_freq_combed_%s.pkl'%ifg_name
        lookup = pkl.load(open(pklfile,'rb'))
        
        # get database frequencies 
        db_dict = analysis.EnergyTerms.AAi_db_lookup(lookup_dir)
        if bbinteraction > 0:
            score = lookup.ix[vdm_name, 'vdm_freq_bb']
            row['f(AAi) bb'] = -np.log10(score / db_dict[vdm_name])
        if scinteraction > 0: 
            score = lookup.ix[vdm_name, 'vdm_freq_sc']
            row['f(AAi) sc'] = -np.log10(score / db_dict[vdm_name])
    return score_df


