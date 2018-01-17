import pickle as pkl
import matplotlib.pyplot as plt

d = [1, 1.2, 1.4, 1.6, 1.8]
f, axarr = plt.subplots(7)
ix = 0
geomdict = pkl.load(open('geomdict.pkl','rb'))
for resi, v  in geomdict.items():
    if (len(v) > 0):
        if resi != 'A753':
            for int_res, nnlist in v:
                axarr[ix].bar(d,nnlist,label='%s with %s'%(resi, int_res), lw=4, alpha=0.4, width=0.2)
            
            axarr[ix].legend()
            ix+=1
#    axarr[subplot].bar(x,hist[0],color=colors[cr],label=labels[cr],lw=4,alpha=0.4,width=1)
#plt.suptitle('%s contacts - raw counts at %s A^2 cutoff' % (ifg,sasa))
plt.show()


