# after getting all the coord pkl files for multiples of 20,
# concatenate them to a megalist

import sys, math
import pickle as pkl

script, aa = sys.argv

look_dir = '/home/gpu/Sophia/combs/st_wd/Lookups/refinedvdms/'

vdm_df = pkl.load(open(look_dir+'vdms_of_{}.pkl'.format(aa),'rb'))
megalist = []
for n in range(1,math.ceil(len(vdm_df)/20000+1)):
    N = n*20
    ls = pkl.load(open(look_dir+'coords_of_{}_multiple_{}.pkl'.format(aa, N), 'rb'))
    megalist += ls

assert len(megalist) == len(vdm_df)
assert len(vdm_df) == len(set([i[0] for i in megalist]))

no_coords = 0
for item in megalist:
    if len(item) != 3:
        no_coords += 1
assert no_coords/len(vdm_df)*100 < 1.5 # make sure it's < 1.5%

pkl.dump(megalist, open(look_dir+'coords_of_{}.pkl'.format(aa),'wb'))
