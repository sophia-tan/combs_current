import matplotlib.pyplot as plt
import sys                              
import pandas as pd
import pickle as pkl
import itertools

seq_effects_lookup = pd.DataFrame(index=[i for i in range(50)], columns=['total_counts', 'bb/sc', 'bb/bb', 'ss/ss', 'resnames', 'ss', 'row_ix'])
# fill df with zeros
for ix, row in seq_effects_lookup.iterrows():
    for col in seq_effects_lookup.columns.values:
        if col=='resnames' or col=='ss' or col=='row_ix':
            seq_effects_lookup.ix[ix, col] = []
        else:   
            seq_effects_lookup.ix[ix, col] = 0

all_distances = []
dist_vdms = pkl.load(open('dist_hbonds_no_repeats.pkl','rb'))

for unique_ifg in set(dist_vdms['iFG_count']):
    vdms = dist_vdms[dist_vdms['iFG_count'] == unique_ifg]
    row_indices = vdms.index.values
    if len(vdms) > 1:
        pairs_vdms = itertools.combinations(row_indices, 2)
        for pair in pairs_vdms: # these are pairs of row indices
            if dist_vdms.ix[pair[0], 'chid_vdm'] == dist_vdms.ix[pair[1], 'chid_vdm']: 

                # calculate separation distance 
                vdmresnum_a = dist_vdms.ix[pair[0], 'resnum_vdm']
                vdmresnum_b = dist_vdms.ix[pair[1], 'resnum_vdm']
                separation = max(vdmresnum_b ,vdmresnum_a) - min(vdmresnum_b, vdmresnum_a)
                if separation < 50:
                    all_distances.append(separation)
                    # add pair of row indices to df for this separation #
                    seq_effects_lookup.ix[separation, 'row_ix'].append(pair)
                    vdmresname_a = dist_vdms.ix[pair[0], 'resname_vdm']
                    vdmresname_b = dist_vdms.ix[pair[1], 'resname_vdm']
                    seq_effects_lookup.ix[separation, 'resnames'].append([vdmresname_a, vdmresname_b])

                    #if 'N' in vdmatoms_a or 'CA' in vdmatoms_a or 'C ' in vdmatoms_a: # bb atoms
                    #    print(vdmatoms_a)
pkl.dump(all_distances, open('seq_separation_dist_hbonds_vdms_list.pkl','wb'))
pkl.dump(seq_effects_lookup, open('temphbondseqeffects.pkl','wb'))
