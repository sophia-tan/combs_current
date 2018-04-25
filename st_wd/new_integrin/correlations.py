import seaborn as sns
import pandas as pd
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
from scipy import stats
from sys import argv

AAs=['R671', 'I673', 'N753','F755','V760', 'E785','R900']
activation = np.array([0.61,0.23,0.54,0.42,0.83,0.47,0.31])

#if method == 'rosetta':
#    calc = [1.53,0.12,1.32,1.76,0.71,1.20,0.54]


for calc in [one,two,three]:
    corr = stats.pearsonr(activation,calc)
    print(corr[0])
    corr = stats.spearmanr(activation,calc)
    print(corr[0])

