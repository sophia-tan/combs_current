# determine what cut-off to use to determine buried vs. not buried

import pandas as pd
import matplotlib.pyplot as plt

for aa in ['alanine', 'lysine', 'aspartate', 'phenylalanine']:
    csvfile = '~/Sophia/combs/st_wd/20180207db_combed_csvs/{}/{}_ifg_atom_density.csv'.format(aa, aa)
    df = pd.read_csv(csvfile)
    convexhull = df['min_hull_dist_CB']
    plt.title(aa)
    plt.hist(convexhull, bins=int(max(convexhull)))
    plt.xlabel('min dist to convex hull')
    #plt.show()
