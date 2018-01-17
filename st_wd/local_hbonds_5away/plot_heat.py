from pprint import pprint
import matplotlib.pyplot as plt
import pickle as pkl
from AAcodes import *
import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.analysis import Analysis
from combs.analysis.Analysis import *
from combs.analysis import PlotFunctions
import pandas as pd
import traceback
import matplotlib
#matplotlib.style.use('ggplot')

######df = pkl.load(open('dist_hbonds.pkl','rb'))
#######ind = Analysis.repeat_indices()
#######df = pd.DataFrame(df, index=ind)
######
######df = Analysis.remove_repeat_proteins(df)
######
######with open('newdist_hbonds_no_repeats.pkl','wb') as f:
######    pkl.dump(df, f)

df = pkl.load(open('dist_hbonds_no_repeats.pkl','rb'))


counts = {} # dict where key = bin, value = counts
for i in [-5,-4,-3,-2,-1,1,2,3,4,5]: # fill in dict with keys
    counts[i] = 0

col_names = ['all helix', 'Nterm strand-helix', 'helix-Cterm turn', 'helix-Cterm strand', 'all strand', 'all loop', 'no ss', '% no ss']

ss_df = pd.DataFrame(index=[-5,-4,-3,-2,-1,1,2,3,4,5], columns=col_names)

# fill with zeros
for ix, row in ss_df.iterrows():
    for col in col_names:
        row[col] = 0

# choose to plot heatmap for all ss, or specific ss
plot_ss = 'all helix'
#plot_ss = 'all'

heatmapdict = PlotFunctions.generate_heatmapdict()

for ix, row in df.iterrows():
    hb_resnums = PlotFunctions.parse_resnums(row['rel_resnums_y']) # for resnums in vdmfrag ifg interacts with 
    for num in hb_resnums:
        if num != 0:
            try:
                # add to counts dict for histogram plot
                counts = PlotFunctions.inc_counts_dict(counts, num)
                # add to heatmapdict for heatmap
                vdm0, vdmi = PlotFunctions.vdm0_vdmi_names(row, num)
                if plot_ss == 'all':
                    heatmapdict = PlotFunctions.inc_heatmap_dict(heatmapdict, num, vdm0, vdmi)
                
                try:
                    ss_df, ss = PlotFunctions.inc_ss_df_Hovmoller(ss_df, row, num)
                    if ss == plot_ss:
                        heatmapdict = PlotFunctions.inc_heatmap_dict(heatmapdict, num, vdm0, vdmi)
                except:
                    ss_df = PlotFunctions.inc_ss_df_Hovmoller(ss_df, row, num)
                
                
                
                #if 
            except Exception:
                traceback.print_exc()
                pass
#print(ss_df.ix[1:5,:])
pkl.dump(ss_df.ix[1:5,:], open('dist_hbonds_ss.pkl','wb'))

PlotFunctions.plot_single_heatmap(heatmapdict, 4)
