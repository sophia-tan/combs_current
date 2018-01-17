from pprint import pprint
import matplotlib.pyplot as plt
import pickle as pkl
import pandas as pd
import matplotlib
from sys import argv
import seaborn as sns
plt.style.use('ggplot')

script, df_path = argv
df = pkl.load(open(df_path,'rb'))

# turn 0's to 1's for log scale
for ix, row in df.iterrows():
    for col in df.columns.values:
        if df.ix[ix, col] == 0:
            df.ix[ix,col] = 1

title = df_path.split('.')[0]
for i in range(4):
    index = [i for i in range(50)]
    plt.bar(index, df['total counts'], width=0.2)
    index = [x+0.2 for x in index]
    plt.bar(index, df['bb/sc'], width=0.2)
    index = [x+0.2 for x in index]
    plt.bar(index, df['sc/sc'], width=0.2)
    index = [x+0.2 for x in index]
    plt.bar(index, df['bb/bb'], width=0.2)
    index = [x+0.2 for x in index]

plt.xlim(0.7,49.9)
plt.ylabel('counts')
plt.title(title)
plt.xticks([i for i in range(50)][1:])
plt.xlabel('green=total counts, pink=bb/sc, red=sc/sc, blue=bb/bb')
plt.yscale('log', nonposy='clip')
plt.show()
