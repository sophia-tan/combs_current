import pandas as pd
import pickle as pkl
from pprint import pprint
from sys import argv
import matplotlib.pyplot as plt, numpy as np

script, combed_aa_path = argv
db_lookup = pkl.load(open('AAi_database_lookups.pkl', 'rb'))
#AAi_lookup = pkl.load(open('AAi_combed_carboxamide_lookup.pkl','rb')) # old combs with 267k
combed_counts = pkl.load(open(combed_aa_path, 'rb'))

######### make df for all freq ######### 
names = [x for x in combed_counts['has_sc freq'].keys()] + ['total counts']
aa_source = ['bb_sc freq', 'has_sc freq', 'db']
display = pd.DataFrame(index=names, columns=aa_source)

######### add db freq to df ######### 
for k, v in db_lookup[1].items():
    display.ix[k, 'db'] = round(v, 2)

######### how to frequencies of combed AAs compare with db frequencies? ######### 
bb_sc_freq = combed_counts['bb_sc freq']
has_sc_freq = combed_counts['has_sc freq']
for x in ['bb_sc freq', 'has_sc freq']:
    dic = combed_counts[x]
    name = x.split(' ')[0] + ' raw_counts' # get name of raw counts dic to get total counts
    display.ix['total counts', x] = sum([val for key, val in combed_counts[name].items()])
    for k, v in dic.items():
        display.ix[k, x] = round(v/db_lookup[1][k], 2)

########### PLOT ############
#sort 
display = display.sort_values(['has_sc freq'], axis=0, ascending=False)

bar_x = np.arange(20)
bar_bb_sc = []
bar_has_sc = []
bar_labels = []
for aa in display.index:
    if aa not in ['total counts', 'MSE', 'CSO', 'SEP', 'SEC', 'TPO']:
        bar_labels.append(aa)
        bar_bb_sc.append(display.loc[aa,'bb_sc freq'])
        bar_has_sc.append(display.loc[aa,'has_sc freq'])
display = display.drop(['db'], axis=1)
display = display.drop(['total counts', 'MSE', 'CSO', 'SEP', 'SEC', 'TPO'])

display = display.rename(index=str, columns={'has_sc freq':'sc interactions'})
display = display.rename(index=str, columns={'bb_sc freq':'bb+sc interactions'})
print(display)
display[['sc interactions', 'bb+sc interactions']].plot.bar()

plt.title('Amino acid preferences for carboxamide vdMs')
plt.ylabel('enrichment ratio')
plt.show()
