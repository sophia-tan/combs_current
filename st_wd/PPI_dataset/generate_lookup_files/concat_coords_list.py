# after getting all the coord pkl files for multiples of 20,
# concatenate them to a megalist

import os

script, aa = sys.argv

look_dir = '/home/gpu/Sophia/combs/st_wd/Lookups/refinedvdms/'


for i in range(1,10):
    print(i*20)
#megalist = []
#lookup = pkl.load(open(look_dir+'vdms_of_{}.pkl'.format(aa),'rb'))
#
#
#
#pkl.dump(coords, open(look_dir+'coords_of_{}_multiple_{}.pkl'.format(aa, multiple_of_20),'wb'))
