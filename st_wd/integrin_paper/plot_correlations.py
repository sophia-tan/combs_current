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
    calc = [128, 43, 115, 108, 0.001, 193, 45, 0.001, 33, 0.001]
    calc = [-0.56442745504304404, 0.33243845991560533, 0.74842252727210479, -0.36632714168594738, -3.0, 0.70577371239096354, -0.98238370794390373, -3.0, 0.31439395722196267, -3.0]

    xaxis = 'interaction geometry score'
    xaxis = 'raw counts'
if method == 'rosetta':
    
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
