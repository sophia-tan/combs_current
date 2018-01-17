import matplotlib.pyplot as plt
import sys                              
import pandas as pd
import pickle as pkl
import itertools

###########seq_effects_lookup = pd.DataFrame(index=[i for i in range(13)], columns=['polar', 'vdw', 'bb/sc', 'bb/bb', 'ss/ss', 'resnames', 'ss', 'row_ix'])
############ fill df with zeros
###########for ix, row in seq_effects_lookup.iterrows():
###########    for col in seq_effects_lookup.columns.values:
###########        seq_effects_lookup.ix[ix, col] = 0
###########
###########all_distances = []
###########dist_vdms = pkl.load(open('distant_vdms_seq_dist_10.pkl','rb'))
###########
###########polar_atoms = 'NO'
###########vdw_atoms = 'CS'
###########
###########for unique_ifg in set(dist_vdms['iFG_count']):
###########    vdms = dist_vdms[dist_vdms['iFG_count'] == unique_ifg]
###########    row_indices = vdms.index.values
###########    if len(vdms) > 1:
###########        pairs_vdms = itertools.combinations(row_indices, 2)
###########        for pair in pairs_vdms: # these are pairs of row indices
###########            if dist_vdms.ix[pair[0], 'chid_vdm'] == dist_vdms.ix[pair[1], 'chid_vdm']: 
###########
###########                # calculate separation distance 
###########                vdmresnum_a = dist_vdms.ix[pair[0], 'resnum_vdm']
###########                vdmresnum_b = dist_vdms.ix[pair[1], 'resnum_vdm']
###########                separation = max(vdmresnum_b ,vdmresnum_a) - min(vdmresnum_b, vdmresnum_a)
###########                if separation < 50:
###########                    all_distances.append(separation)
###########
###########                # determine if this is a polar or vdw interaction
###########                #if vdm_atom
###########                vdmatoms_a = dist_vdms.ix[pair[0], 'atom_names_vdm'].split(' ')
###########                vdmatoms_b = dist_vdms.ix[pair[1], 'atom_names_vdm'].split(' ')
###########                vdmatoms = vdmatoms_a + vdmatoms_b
###########                for atom in vdmatoms:
###########                    if atom[0] in polar_atoms:
###########                        pass
###########                    elif atom[0] in vdw_atoms:
###########                        pass
###########                    else:
###########                        print(atom)
###########
###########                
###########                #if 'N' in vdmatoms_a or 'CA' in vdmatoms_a or 'C ' in vdmatoms_a: # bb atoms
###########                #    print(vdmatoms_a)
   
all_distances = pkl.load(open('seq_separation_all_dist_vdms.pkl','rb'))
plt.hist(all_distances, bins=49)
plt.title('raw counts for separations between vdms interacting w/ same carboxamide')


separation_num = 5
skip=10
constant=0.001

###y_axis = []
###for separation_num in range(50):
###    exp_val = 1959620 - skip * 2 * 8667 - separation_num * 8667
###    y_axis.append(exp_val*constant)
###
####plt.plot([i for i in range(50)], y_axis)
###
###


#plt.xlim(0,50)




plt.yscale('log', nonposy='clip')
plt.show()

## need to keep track of bb, sc, polar/vdw, keep resnames for heatmap and ss
