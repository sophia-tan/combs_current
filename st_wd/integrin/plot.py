import sys
import matplotlib.pyplot as plt
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
import pickle as pkl, numpy as np

geomdict=pkl.load(open('newgeomdict.pkl','rb'))

f, axarr = plt.subplots(8,sharex=True)
ix=0
for (targetresi,targetresn), v  in geomdict.items():
    targetresn = constants.three_letter_code[targetresn]
    for int_res, rmsdlist in v.items():
        axarr[ix].hist(rmsdlist[0],bins=50,label=targetresi+targetresn+'-'+int_res,histtype='step',lw=3)
        #axarr[ix].hist(rmsdlist[0],bins=50,label=targetresi+targetresn+'-'+int_res,histtype='step',lw=3,log=True)
        axarr[ix].legend()
    ix += 1
#            
##plt.suptitle('%s contacts - raw counts at %s A^2 cutoff' % (ifg,sasa))
#for target,int_res in geomdict.items():
#     print(int_res)
#
plt.xlabel('interactamer rsmd/atom')
plt.suptitle('integrin hotspots')
plt.show()
