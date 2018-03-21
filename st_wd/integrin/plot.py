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

def apply_cutoff(ifgrmsds,vdmrmsds,cutoff):
    '''alignment of query onto lookup ifg may be off, so throw away
    data points if it's > this cutoff'''
    ls=[]
    for ix, num in enumerate(ifgrmsds):
        if num<=cutoff:
            ls.append(ix)
    return [ifgrmsds[i] for i in ls], [vdmrmsds[i] for i in ls]
    
f, axarr = plt.subplots(7, 2, sharex='col', sharey='row')
row = 0
ix = 0

for tag in geodic:
    k,v=tag,geodic[tag]
    targetresi,targetresn = k
    targetresn = constants.three_letter_code[targetresn]
    # sort order of v
    sortedv = [(k,v) for k,v in v.items()]
    sortedv = collections.OrderedDict(sortedv)

    for int_res, rmsdlists in sortedv.items():
        for col,typ in enumerate(rmsdlists): # first typ is chA as ifg, second is chB as ifg
            for rmsdlist in typ:
                ifg,vdm,ifgrmsds,vdmrmsds=rmsdlist
                ifgrmsds, vdmrmsds=apply_cutoff(ifgrmsds,vdmrmsds,cutoff_float)
                if ifg=='backbone':
                    ifg='bb'
                else:
                    ifg='sc'
                if vdm=='backbone':
                    vdm='bb'
                else:
                    vdm='sc'
                axarr[row,col].hist(vdmrmsds,bins=100,label=targetresi+' '+targetresn+ifg+'-'+int_res+vdm,histtype='step',lw=2)
            if ix%2==1:
                axarr[row,col].legend(scatterpoints=100,fontsize=8,loc='upper right')
            ix += 1
    row += 1

f.text(0.5, 0.04, 'vdM rmsd/atom', ha='center')

#plt.xlabel('interactamer rsmd/atom')
plt.suptitle('integrin hotspot interactamers with {}A cutoff'.format(cutoff))
plt.show()
