# print out counts for how many interactamers are within 0.5 rmsd

import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
import pickle as pkl, numpy as np
import collections

script, cutoff=sys.argv
cutoff_float = float(cutoff)

geomdict=pkl.load(open('ifgs_integrin_geom_dict.pkl','rb'))

geodic = [(k,v) for k,v in geomdict.items()]
geodic = collections.OrderedDict(sorted(geodic))

def apply_cutoff(rmsds,cutoff):
    ls=[]
    for ix, num in enumerate(rmsds):
        if num<=cutoff:
            ls.append(ix)
    return [rmsds[i] for i in ls]
    

dic = {}
for tag in geodic:
    print('----------------------------')
    k,v=tag,geodic[tag]
    targetresi,targetresn = k
    targetresn = constants.three_letter_code[targetresn]
    # sort order of v
    sortedv = [(k,v) for k,v in v.items()]
    sortedv = collections.OrderedDict(sortedv)
    dic[tag]= {}

    for int_res, rmsdlists in sortedv.items():
        dic[tag][int_res] = []
        # there are 2 lists in rmsdlists. first is with combed vdMs of chA as ifg,
        # second is with combed vdMs of chB as ifg
        int_resn, intrxns = rmsdlists[0], rmsdlists[1]
        for ix_pair,pair in enumerate(intrxns):
            ifg,vdm=pair[0],pair[1]
            rmsds=rmsdlists[1][ix_pair][2]+rmsdlists[2][ix_pair][2]
            rmsds=apply_cutoff(rmsds,cutoff_float)
            
            if ifg=='backbone':
                ifg='bb'
            elif 'FG' in ifg:
                ifg='fg'
            else:
                ifg='sc'
            if vdm=='backbone':
                vdm='bb'
            elif 'FG' in vdm:
                vdm='fg'
            else:
                vdm='sc'
            print(targetresi+' '+targetresn+ifg+'-'+int_resn+vdm)
            print(len(rmsds))

