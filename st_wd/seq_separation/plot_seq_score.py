import pickle as pkl
import matplotlib.pyplot as plt
from pprint import pprint
import numpy as np
from sys import argv
import pandas as pd

script, carbox_df = argv

lengths = pd.read_pickle(open('db_chain_lengths.pkl','rb'))
dic = {} # key = length, val = counts

for i in lengths:
    try: 
        dic[i] += 1
    except:
        dic[i] = 1

##### plot distribution of chain lengths ###
#x = []
#y = []
#for key, val in dic.items():
#    x.append(key)
#    y.append(val)
#plt.bar(x,y)
#plt.show()

#### store how many windows exist for each sep # for db chains ########
all_windows = {} # key = sep # , val = counts of all windows summed over every chain
for sep in range(1, 101):
    all_windows[sep] = 0
    for key, val in dic.items():
        if sep > key:
            pass
        else:
            windows = key - sep
            windows = windows * val
            all_windows[sep] += windows

###### plot # windows in db ##########
#sep = []
#windows = []
#for k, v in all_windows.items():
#    sep.append(k) 
#    windows.append(v)
#    
#plt.xlabel('separation #')
#plt.ylabel('# windows')
#plt.bar(sep, windows)
#plt.show()
########################################################################
##### compare # windows at each sep # to reference state (sep # 50)
f_inc = all_windows
for key, val in f_inc.items():
    f_inc[key] = val / 1526818 # num windows at sep 50

##### plot distribution of # windows compared to ref state
sep = []
frac = []
for k, v in f_inc.items():
    sep.append(k) 
    frac.append(v)
#plt.xlabel('separation #')
#plt.ylabel('f(included)')
#plt.bar(sep, frac)
#plt.show()
########################################################################
##### good stuff, now look at carboxamide: divide f(obs)/f(exp) 
df = pd.read_pickle(open(carbox_df, 'rb'))
obs = df['total counts']
exp = [f*obs[49] for f in frac]
plt.title('obs/exp interactions, no skip hbond vdms, cbrt correction')

score = []
for i in range(49):
    correction = np.cbrt(i) / np.cbrt(50)
    score.append(obs[i]/exp[i] * correction)

plt.bar([i for i in range(49)], score)
plt.xlabel('D (separation #)')
plt.ylabel('enrichment ratio')
plt.show()
