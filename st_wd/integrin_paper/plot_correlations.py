import seaborn as sns
import pandas as pd
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
from scipy import stats
from sys import argv
sns.set_style("white")

f,(ax1,ax2) = plt.subplots(1,2,sharey=True)

#for method,ax in zip(['rosetta','pdb'],[ax1,ax2]):
for method,ax in zip(['pdb'],[ax1,ax2]):
    #AAs=['R671', 'I673', 'N753','F755','S758','V760', 'E785','H787','R900','L959', \
    #     'D552', 'Y594', 'T603','H626','K658', 'V664']
    #activation = np.array([0.54,0.23,0.37,0.4,0.64,0.83,0.47,0.17,0.31,0.05, \
    #    .12, .44, .399, .10, .20, .42])
    #
    
    AAs=['R671', 'I673', 'N753','F755','S758','V760', 'E785','H787','R900','L959', \
         'D552', 'Y594', 'T603','H626','K658', 'V664', 'E534']
    activation = np.array([0.54,0.23,0.37,0.4,0.64,0.83,0.47,0.17,0.31,0.05, \
        .12, .44, .399, .10, .20, .42, 0.76])
    
    
    
    sortedactivation = sorted(activation)
    
    
    if method == 'pdb':
        new_order = [np.where(activation==i)[0][0] for ix,i in enumerate(sortedactivation)]

        calc = [-2.4329709040267118, -3.4090147995618105, -2.6006553576794031, -2.8004132052090203, -4.0, -3.381630716842734, -2.7327972824759992, -4.0, -2.9473636913326939, -4.0, -4.0, -3.0250766057272291, -2.5602510798264042, -4.0, -3.031105362355941, -3.5854607295085006, -3.6240241246478613]
        xaxis = 'interaction score'
        #yerr = [.03,.03,.03,.13,.08,.12,.04,.02,.06,.001,0,0,0,0,0,0] 
        #yerr = [yerr[ix] for ix in new_order]
        col = 'c'
    
    if method == 'rosetta':
        new_order = [np.where(activation==i)[0][0] for ix,i in enumerate(sortedactivation)]
        xaxis = 'Rosetta DDG (kcal/mol)'
        calc = [1.53,.12,1.32,1.76,.56,.71,1.2,.06,.54,0,\
        0, .6, 2.89, 0, .63, .34]
        col = 'lavender'
    
    calc = [calc[ix] for ix in new_order]
    AAs = [AAs[ix] for ix in new_order]
    activation = [activation[ix] for ix in new_order]
    
    def best_fit_slope(xs,ys):
        xs,ys=np.array(xs),np.array(ys)
        m = (((np.mean(xs)*np.mean(ys)) - np.mean(xs*ys)) / \
                 ((np.mean(xs)*np.mean(xs)) - np.mean(xs*xs)))
        b = np.mean(ys) - m*np.mean(xs)
        print('m',m)
        print('b',b)
        return m,b
    
    m,b=best_fit_slope(calc,activation)
    regression_line = [(m*x)+b for x in calc]
    #x=[0.001,11.875]
    #y = [0.18405186840597873,0.67299933809206247]
    x=[-1.1,2]
    y=[0.17258020395350002, 0.582692403831]

    corr = stats.pearsonr(activation,calc)
    print(corr[0],corr[1])
    
    def text(label,x,y):
        xytext=(5,-4.5)
        if label=='H626':
            xytext=(5,-9)
        if ax==ax1 and label=='V664':
            xytext=(5,-9)
        #if ax==ax2 and label=='V760':
        #    xytext=(0,-15)
        if ax==ax1 and label=='T603':
            xytext=(-10,-15)

        ax.annotate(
            label,xy=(x,y),
            xytext=xytext,textcoords='offset points',size=12)
    
    sns.despine(right=True)

    if method=='pdb':
        #plt.errorbar(calc,activation,yerr=yerr,fmt='o')
        ax.plot(x,y,ls=':')
    ax.set_xlabel(xlabel=xaxis,fontdict={'size':12})
    
    for label, x, y in zip(AAs, calc,activation):
        text(label,x,y)
    
    ax.scatter(calc,activation,marker='o',edgecolors='gray',facecolors=col,
            lw=1)
    #if ax==ax2:
    ax.scatter(calc[14],activation[14],marker='o',edgecolors='black',facecolors='r',
        lw=1) # overwriting with diff color

f.text(0.015, 0.5, 'activation index', ha='center', va='center', rotation='vertical',fontsize=12)
f.text(0.33,0.15,'$R^2=0.42$\np-value = 0.12',size=12,multialignment='left',bbox={'edgecolor':'orange','pad':3,'facecolor':'none','lw':2})
f.text(0.8,0.15,'$R^2=0.85$\np-value = 0.000007',size=12,multialignment='left',bbox={'edgecolor':'orange','pad':3,'facecolor':'none','lw':2})
f.text(0.55,0.755,'(not included\nin correlation)',size=8,multialignment='left')
f.text(0.15,0.755,'(not included\nin correlation)',size=8,multialignment='left')

f.subplots_adjust(wspace=0)
plt.tight_layout()
plt.show()
