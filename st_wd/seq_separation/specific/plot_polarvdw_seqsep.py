import matplotlib.pyplot as plt
import sys                              
import pandas as pd
import pickle as pkl
import itertools

df = pkl.load(open('disttempseqeffects.pkl','rb'))
plt.title('raw counts for separations between distant vdms making polar interactions w/ same carboxamide')

counts = df['polar']
plt.bar([i for i in range(50)], counts)

separation_num = 5
skip=10
constant=0.001

#y_axis = []
#for separation_num in range(50):
#    exp_val = 1959620 - skip * 2 * 8667 - separation_num * 8667
#    y_axis.append(exp_val*constant)
#
#plt.plot([i for i in range(50)], y_axis)
#
#
#plt.xlim(0,10)


plt.yscale('log', nonposy='clip')
plt.show()

