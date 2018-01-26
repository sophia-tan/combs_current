# gets cb sasa scores for each vdm of a FG for a design
# commandline arguments: pdb of design, txtfile of vdm_resnum vdm_aa 
# make sure there's a SINGLE SPACE b/n resnum and aa 
# and aa is 3 letters

import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
import prody as pr
import freesasa, math
import numpy as np, pickle as pkl, pandas as pd

script, design_pdb, fg_vdm_txt= sys.argv

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
    columns = ['iFG', 'AA', 'sasa_cb_3a', 'burial score', 'omega', 
                'dunbrack', 'rama']
    df = pd.DataFrame(columns=columns, index=pd.Series(vdms_resnum, name='resnum'))
    # fill in df
    for ix in range(len(vdms_aa)):
        df.loc[vdms_resnum[ix], 'iFG'] = ifg[ix]
        df.loc[vdms_resnum[ix], 'AA'] = vdms_aa[ix]
    return df

############################################################################################
########################### SASAs at cb 3A probe ###########################################
############################################################################################

# 1) get dictionary with vdm AAs, resnums, and sasas at cb 3A probe ################
def cb_sasas(design_pdb, fg_vdm_txt):
    sasa_dict = {} # key is resnum, value is list [aa, cbsasa]
    
    # get the vdm AAs and resnums from txtfile
    vdms_aa = []
    vdms_resnum = []
    with open(fg_vdm_txt) as inF:
        for line in inF:
            line = line.strip()
            line = line.split(' ')
            vdms_resnum.append(line[0])
            vdms_aa.append(line[1])

    # parse design and do freesasa calc
    prody_parsed = pr.parsePDB(design_pdb, altloc='A', model=1)
    fs_struct = freesasa.Structure(design_pdb) # more atoms 
    fs_result = freesasa_cb(prody_parsed, probe_radius=3) # less atoms bc this is Cb cutoff 

    # get sasa
    for resnum, aa in zip(vdms_resnum, vdms_aa): # shouldn't have to worry about 
            #neg resnums in designed proteins
        # get Cb atoms 
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

sasadict = cb_sasas(design_pdb, fg_vdm_txt)

# 2) 'find sasa bin for that vdm and score' ###########################

def scoring_sasa(sasadict):
    # load lookup table
    directory = '/home/gpu/Sophia/combs/st_wd/sasa/'
    lookup = pkl.load(open(directory+'scores_largeprobe_sasa_db_lookup.pkl','rb'))
    #lookup = pkl.load(open(directory+'scores_largeprobe_sasa_db_lookup.pkl'%ifg,'rb'))
    # take -np.log10 for all values in df, and drop rows 150-180
    lookup = lookup.drop([160, 170, 180, 190])
    lookup = lookup.applymap(lambda g: -np.log10(g))
    for key, value in sasadict.items(): # key is resnum
        aa, sasa = value[0], value[1]
        if sasa > 150: # this site is exposed
            score = -1
        else:
            if sasa == 0:
                sasabin = 10
            else:
                sasabin = int(math.ceil(sasa / 10.0)) * 10 
    
            score = np.round(lookup.loc[sasabin,aa],2)


    
