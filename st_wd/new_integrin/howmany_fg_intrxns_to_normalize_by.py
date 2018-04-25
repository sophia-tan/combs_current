import sys, os
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis, apps
#from ScoringFunctions import *
#from PyrosettaScores import *
from residues_integrin import *
import prody as pr
import numpy as np, pickle as pkl, pandas as pd
from pprint import pprint
from itertools import *

direc = '/home/gpu/Sophia/combs/st_wd/integrin/output_data/'
csv_1 = '/home/gpu/Sophia/STcombs/20171118/{}/csv/'
csv_2 = '/home/gpu/Sophia/combs/st_wd/20180207db_combed_csvs/{}/'

def fg_contacting(row,vdm,ifg):
    ifgcontacts = []
    vdmcontacts = []
    dist = row['dist_info'].lstrip('(')
    dist = dist.rstrip(')')
    dist = dist.split(') (')
    for d in dist:
        d=d.split(' ')
        i,v,ang=d
        if float(ang) <= 3.5:
            vdmcontacts.append(v)
            ifgcontacts.append(i)
        elif float(ang) <= 4.8 and i[0]=='C' and v[0]=='C':
            vdmcontacts.append(v)
            ifgcontacts.append(i)
    vdminter = len(set(apps.constants.ifg_atoms[vdm]).intersection(set(vdmcontacts)))
    ifginter = len(set(apps.constants.ifg_atoms[ifg]).intersection(set(ifgcontacts)))
    if vdminter > 0 and ifginter>0:
        return True
    else:
        return False

for pklf in os.listdir(direc):
    ifg,vdm = pklf.split('_')[1],pklf.split('_')[3]
    if ifg=='backbone':
        ifg='glycine'
    if vdm=='backbone':
        vdm='glycine'
    try:
        an = analysis.Analyze(csv_1.format(ifg))
    except:
        an = analysis.Analyze(csv_2.format(ifg))
    lookup = pkl.load(open(direc+pklf,'rb'))[2]
    vdmsdir = '/home/gpu/Sophia/combs/st_wd/Lookups/refinedvdms/'
    refinedvdms = pkl.load(open(vdmsdir+'vdms_of_{}.pkl'.format(ifg),'rb'))
    vdm=apps.constants.AAname_rev[vdm]
    ifg=apps.constants.AAname_rev[ifg]
    refinedvdms = refinedvdms[refinedvdms['resname_vdm']==vdm]
    merged = refinedvdms.merge(an.ifg_contact_vdm,left_index=True,right_index=True)
    #merged = merged[:10] # delete
    unsuccessful = len(merged)-len(lookup) # couldn't get coords extracted
    # count how many of those vdms have FG groups interacting
    fg_intrxns = merged.apply(fg_contacting,vdm=vdm,ifg=ifg,axis=1)
    print(ifg,vdm)
    print(sum(fg_intrxns))
