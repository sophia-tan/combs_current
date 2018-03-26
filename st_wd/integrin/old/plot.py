# abandoned this bc just print out numbers
import sys
import matplotlib.pyplot as plt
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
import pickle as pkl, numpy as np
import collections

script, cutoff=sys.argv
cutoff_float = float(cutoff)

geomdict=pkl.load(open('integrin_geom_dict.pkl','rb'))

geodic = [(k,v) for k,v in geomdict.items()]
geodic = collections.OrderedDict(sorted(geodic))

def apply_cutoff(rmsds,cutoff):
    ls=[]
    for ix, num in enumerate(rmsds):
        if num<=cutoff:
            ls.append(ix)
    return [rmsds[i] for i in ls]
    
f, axarr = plt.subplots(8, sharex='col', sharey='row')
ix = 0

for tag in geodic:
    k,v=tag,geodic[tag]
    targetresi,targetresn = k
    targetresn = constants.three_letter_code[targetresn]
    # sort order of v
    sortedv = [(k,v) for k,v in v.items()]
    sortedv = collections.OrderedDict(sortedv)

    for int_res, rmsdlists in sortedv.items():
        # there are 2 lists in rmsdlists. first is with combed vdMs of chA as ifg,
        # second is with combed vdMs of chB as ifg
        for rmsdlist in rmsdlists:
            rmsds = []
            for pair in rmsdlist:
                ifg,vdm=pair[0],pair[1]
                rmsds += pair[2] # pair[3] is num_atoms
            rmsds=apply_cutoff(rmsds,cutoff_float)
            if ifg=='backbone':
                ifg='bb'
            else:
                ifg='sc'
            if vdm=='backbone':
                vdm='bb'
            else:
                vdm='sc'
            axarr[ix].hist(rmsds,bins=30,label=targetresi+' '+targetresn+ifg+'-'+int_res+vdm,histtype='step',lw=2)
            axarr[ix].legend(scatterpoints=30,fontsize=8,loc='upper left')
    ix += 1

f.text(0.5, 0.04, 'interactamer rmsd',ha='center')

plt.suptitle('integrin hotspot interactamers with {}A cutoff'.format(cutoff))
plt.show()
