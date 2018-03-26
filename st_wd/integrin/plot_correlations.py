import seaborn as sns
import pandas as pd
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
from scipy import stats
from sys import argv

activation = [0.61,0.23,0.54,0.42,0.83,0.47,0.31]
#data1 = [114,20,137,66,391,99,69]
#data2 = [0.1,0.065,0.371,0.047,0.823,0.152,0.211]
#data3 = [114,20,89,62,391,99,69]
data4 = [0.224,0.151,0.423,0.085,1.04,0.339,0.424]

#corr = stats.pearsonr(activation,data1)
#print(corr[0])
#corr = stats.pearsonr(activation,data2)
#print(corr[0])
#corr = stats.pearsonr(activation,data3)
#print(corr[0])
corr = stats.pearsonr(activation,data4)
print(corr[0])

AAs=['R671', 'I673', 'N753','F755','V760', 'E785','R900']
for label, x, y in zip(AAs, activation, data4):
    plt.annotate(
        label,
        xy=(x, y), xytext=(10, 10),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.1', fc='yellow', alpha=0.2),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

plt.xlabel('activation index')
plt.ylabel('interaction geometry score')


plt.scatter(activation,data4)
plt.show()
