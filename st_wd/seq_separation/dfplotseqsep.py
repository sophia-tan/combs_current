from pprint import pprint
import matplotlib.pyplot as plt
import pickle as pkl
import pandas as pd
import matplotlib
from sys import argv
import seaborn as sns

script, df_path = argv
df = pkl.load(open(df_path,'rb'))

# turn 0's to 1's for log scale
for ix, row in df.iterrows():
    for col in df.columns.values:
        if df.ix[ix, col] == 0:
            df.ix[ix,col] = 1

title = df_path.split('.')[0]
df.plot(kind='bar', logy=True)
plt.xlim(1,49)
plt.title(title)
#plt.yscale('log', nonposy='clip')
plt.show()
