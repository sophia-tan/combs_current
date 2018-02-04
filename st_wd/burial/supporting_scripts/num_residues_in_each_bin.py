import matplotlib.pyplot as plt, numpy as np, pandas as pd
import pickle as pkl

# load df with burial distances of each residue in db
df = pkl.load(open('../database_mindistconvexhull.pkl','rb'))
df = df['burialdist']
df = df[:20000]

def in_bin(val, upperbound):
    if val >= upperbound-1 and val < upperbound:
        return True
    else:
        return False

total = len(df)
index = [x for x in range(33)][1:]
cum_sum = 0
ls = []
ls_perc = []
for ind in index:
    num_res=df.apply(in_bin, upperbound=ind).sum()
    cum_sum += num_res
    ls.append(num_res)
    ls_perc.append(np.round(cum_sum/total*100,1))

bar = plt.bar(height=ls,left=index)

def autolabel(rects):
    """
    Attach a text label above each bar displaying its height
    """
    for ix, rect in enumerate(rects):
        height = rect.get_height()
        perc = ls_perc[ix]
        plt.text(rect.get_x() + rect.get_width()/2., 1.02*height,
                '%d' % perc,  ha='center', va='bottom')

autolabel(bar)
plt.ylabel('counts')
plt.xlabel('upper bound of burial bin')
plt.title('# of database residues in burial bin')
plt.show()
