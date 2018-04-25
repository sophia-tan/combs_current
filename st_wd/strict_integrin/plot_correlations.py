import seaborn as sns
import pandas as pd
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
from scipy import stats
from sys import argv

method = 'pdb'
#method = 'rosetta'

AAs=['R671', 'I673', 'N753','F755','S758','V760', 'E785','H787','R900','L959']
activation = np.array([0.54,0.23,0.37,0.4,0.64,0.83,0.47,0.17,0.31,0.05])
sortedactivation = sorted(activation)
new_order = [np.where(activation==sortedactivation[ix])[0][0] for ix,i in enumerate(activation)]

if method == 'pdb':
    #calc = [0.5876814386243279, 0.4032804566044927, 0.22163804152906882, 0.4291674761782889, 0.5912912001834979, 0.9126884158047551, 0.18211558178879497, 0.2621639109492408, 0.2681702900399585, 0.05829681858255777] # direct_intrxn .5
    calc = [-4.2771133874742526, -5.7181631524917833, -4.8999582956655283, -4.0109998525647814, -4.2159223456810651, -2.6056650167654216, -5.807820281775248, -4.1060426799463441, -5.8992296015819585, -7.373757685249247]
    
    xaxis = 'interaction geometry score'
    xaxis = 'raw counts'
if method == 'rosetta':
    calc = [1.53,0.12,1.32,1.76,0.56,0.71,1.20,0.06,0.54,0]
    
    xaxis = 'Rosetta DDG'
calc = [calc[ix] for ix in new_order]

AAs = [AAs[ix] for ix in new_order]
activation = sortedactivation

def best_fit_slope(xs,ys):
    xs,ys=np.array(xs),np.array(ys)
    m = (((np.mean(xs)*np.mean(ys)) - np.mean(xs*ys)) / \
             ((np.mean(xs)*np.mean(xs)) - np.mean(xs*xs)))
    b = np.mean(ys) - m*np.mean(xs)
    return m,b

m,b=best_fit_slope(calc,activation)
regression_line = [(m*x)+b for x in calc]
plt.plot(calc,regression_line)
corr = stats.pearsonr(activation,calc)
print(corr[0])
#corr = stats.spearmanr(activation,calc)
#print(corr[0])

def plot(label,x,y):
    plt.annotate(
        label,
        xy=(x, y), xytext=(10, 10),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.1', fc='white', alpha=0.2),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

for label, x, y in zip(AAs, calc,activation):
    plot(label,x,y)

plt.ylabel('activation index')
plt.xlabel(xaxis)

plt.scatter(calc,activation)
plt.show()
