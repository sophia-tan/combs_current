import seaborn as sns
import pandas as pd
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
from scipy import stats
from sys import argv

method = 'pdb'
#method = 'rosetta'

AAs=['R671', 'I673', 'N753','F755','V760', 'E785','R900']
activation = np.array([0.61,0.23,0.54,0.42,0.83,0.47,0.31])
sortedactivation = sorted(activation)
new_order = [np.where(activation==sortedactivation[ix])[0][0] for ix,i in enumerate(activation)]

if method == 'pdb':
    calc = [0.25992319688762933, 0.13085399449035812, 0.25539677713590758, 0.11906665955820937, 1.6150566346526551, 0.2936227772095939, 0.59271803556308211]
    #calc = [0.224,0.151,0.423,0.085,1.04,0.339,0.424]
    xaxis = 'interaction geometry score (%)'
if method == 'rosetta':
    calc = [1.53,0.12,1.32,1.76,0.71,1.20,0.54]
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
#corr = stats.pearsonr(activation,calc)
corr = stats.spearmanr(activation,calc)
print(corr[0])

for label, x, y in zip(AAs, calc,activation):
    plt.annotate(
        label,
        xy=(x, y), xytext=(10, 10),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.1', fc='white', alpha=0.2),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

plt.ylabel('activation index')
plt.xlabel(xaxis)

plt.scatter(calc,activation)
plt.show()
