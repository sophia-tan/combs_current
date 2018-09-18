import seaborn as sns
import pandas as pd
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
from scipy import stats
from sys import argv
sns.set_style("white")

f,(ax1,ax2) = plt.subplots(1,2)

colordict = {'chA':'orange', 'chB':'mediumturquoise'}
for method,ax in zip(['rosetta','pdb'],[ax1,ax2]):
    AAs=['R671', 'I673', 'N753','F755','S758','V760', 'E785','H787','R900', 'D552', 'Y594', 'T603','H626','K658', 'V664']
    chain=['chA', 'chA', 'chA', 'chA', 'chA', 'chA', 'chA', 'chA', 'chA', 
            'chB', 'chB', 'chB', 'chB', 'chB', 'chB']

    colors = [colordict[x] for x in chain]
    activation = [0.54, 0.23, 0.37, 0.4, 0.64, 0.73, 0.47, 0.17, 0.31, 0.18, 0.46, 0.82, 0.16, 0.44, 0.54]

    #sortedactivation = sorted(activation)
    
    if method == 'pdb':
        #new_order = [np.where(activation==i)[0][0] for ix,i in enumerate(sortedactivation)]
        
        calc = [0.020514348923318413, 0.0029861888764464353, 0.012941572602034023, 0.020118966976494766, 0.0, 0.01819351281815676, 0.011952744745410703, 0.0, 0.006519967400162999, 0.0, 0.0067537532663049585, 0.030234725967370472, 0.0, 0.006134969325153374, 0.027692705862054093]

        xaxis = 'E(h)'
    
    if method == 'rosetta':
        #new_order = [np.where(activation==i)[0][0] for ix,i in enumerate(sortedactivation)]
        xaxis = 'Rosetta ddG (kcal/mol)'
        calc = [1.72, 0.23, 0.8, 1.08, -0.04, 0.38, 1.57, -0.17, 1.49,  0.03, 1.71, 1.5, -0.0, 0.38, -0.09]
    
    #calc = [calc[ix] for ix in new_order]
    #AAs = [AAs[ix] for ix in new_order]
    #colors = [colors[ix] for ix in new_order]
    #activation = [activation[ix] for ix in new_order]

    def best_fit_slope(xs,ys):
        xs,ys=np.array(xs),np.array(ys)
        m = (((np.mean(xs)*np.mean(ys)) - np.mean(xs*ys)) / \
                 ((np.mean(xs)*np.mean(xs)) - np.mean(xs*xs)))
        b = np.mean(ys) - m*np.mean(xs)
        print('m',m)
        print('b',b)
        return m,b
    
    print('############')

    activation_no_ser  = [n for aa, n in zip(AAs, activation) if aa != 'S758']
    print('activation no ser')
    print(activation_no_ser)
    calc_no_ser  = [n for aa, n in zip(AAs, calc) if aa != 'S758']
    print('calc no ser')
    print(calc_no_ser)
    corr = stats.pearsonr(activation_no_ser, calc_no_ser)
    print(corr[0],corr[0]**2)
    
    m,b=best_fit_slope(calc_no_ser,activation_no_ser)
    regression_line = [(m*x)+b for x in calc_no_ser]
    ran = max(calc_no_ser)-min(calc_no_ser)
    x=[min(calc_no_ser)-(ran/10), max(calc_no_ser)+(ran/6)]
    y=[m*n+b for n in x]
    
    def text(label,x,y):
        xytext=(5,-4.5)
        if ax==ax2 and label=='H626':
            xytext=(-32,-4)
        if ax==ax1 and label=='R671':
            xytext=(-32,-4)
        if ax==ax2 and label=='Y594':
            xytext=(-32,-4)
        if ax==ax1 and label=='E785':
            xytext=(5,2)
        if ax==ax2 and label=='D552':
            xytext=(5,2)
        if ax==ax1 and label=='H787':
            xytext=(-12,7)

        ax.annotate(
            label,xy=(x,y),
            xytext=xytext,textcoords='offset points',size=12)
    
    sns.despine(right=True)

    if method=='pdb':
        ax.plot(x,y,ls=':')
    ax.set_xlabel(xlabel=xaxis,fontdict={'size':12})
    
    for label, x, y in zip(AAs, calc,activation):
        text(label,x,y)
    
    ax.scatter(calc,activation,marker='o',edgecolors='black',facecolors=colors,
            lw=1)
    ax.scatter(calc[4],activation[4],marker='o',edgecolors='black',facecolors='r',
        lw=1) # overwriting with diff color

######### below includes a way to display p-val ###########
f.text(0.015, 0.5, 'activation index', ha='center', va='center', rotation='vertical',fontsize=12)
f.text(0.52, 0.5, 'activation index', ha='center', va='center', rotation='vertical',fontsize=12)
#f.text(0.33,0.15,'$R^2=0.52$\np-value = 0.04',size=12,multialignment='left',bbox={'edgecolor':'orange','pad':3,'facecolor':'none','lw':2})
#f.text(0.8,0.15,'$R^2=0.74$\np-value = 0.001',size=12,multialignment='left',bbox={'edgecolor':'orange','pad':3,'facecolor':'none','lw':2})
f.text(0.39,0.15,'$R^2=0.21$' ,size=12,multialignment='left',bbox={'edgecolor':'black','pad':3,'facecolor':'none','lw':1})
f.text(0.86,0.15,'$R^2=0.73$',size=12,multialignment='left',bbox={'edgecolor':'black','pad':3,'facecolor':'none','lw':1})
f.text(0.58,0.73,'(not included\nin correlation)',size=8,multialignment='left')
f.text(0.11,0.73,'(not included\nin correlation)',size=8,multialignment='left')

f.subplots_adjust(wspace=0)
plt.tight_layout()
plt.show()
