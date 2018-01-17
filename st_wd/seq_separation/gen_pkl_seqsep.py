# command line arguments:
# df_path = path to vdms (csv output converted to df pkl)
# interaction_type: can be hbond, polar, vdw, or all

from sys import argv
from SeqSeparationFunctions import *

script, df_path, interaction_type = argv

vdms_df = pkl.load(open(df_path, 'rb'))

# are these all interactions, or hbond interactions?

seqsep_df = gen_seq_sep_df(interaction_type)

for unique_ifg in set(vdms_df['iFG_count']):
    vdms = vdms_df[vdms_df['iFG_count'] == unique_ifg]
    row_indices = vdms.index.values
    if len(vdms) > 1:
        pairs_vdms = itertools.combinations(row_indices, 2)
        for pair in pairs_vdms: # these are pairs of row indices
            if vdms_df.ix[pair[0], 'chid_vdm'] == vdms_df.ix[pair[1], 'chid_vdm']: 
                # calculate separation distance 
                separation = get_sep(vdms_df, pair)
                if separation < 50:
                    # check interaction type
                    if interaction_type == 'hbond' or interaction_type == 'all':
                        pass # because all rows in hbond df are hbond
                    elif interaction_type == 'polar' or interaction_type == 'vdw':
                        int_type = check_type(vdms_df, pair)
                        if interaction_type not in int_type:
                            break

                    # now we know these cooperating vdms are < 50 and the right int type!
                    # add pair of row indices to df for this separation #
                    seqsep_df.ix[separation, 'row_ix'].append(pair)
                    vdmresname_a = vdms_df.ix[pair[0], 'resname_vdm']
                    vdmresname_b = vdms_df.ix[pair[1], 'resname_vdm']
                    seqsep_df.ix[separation, 'resnames'].append([vdmresname_a, vdmresname_b])
                    seqsep_df.ix[separation, 'total counts'] += 1
                        
                    # determine whether the interaction is bb/bb, bb/sc, or sc/sc
                    bb_sc_list = bb_or_sc(vdms_df, pair, interaction_type)
                    for p in bb_sc_list:
                        seqsep_df.ix[separation, p] += 1

fname = df_path.split('.')[0] 
pkl.dump(seqsep_df, open('%s_%s_seqeffects.pkl' % (interaction_type, fname),'wb'))
