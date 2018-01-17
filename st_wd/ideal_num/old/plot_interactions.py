# assumes that ifg_interaction_types.pkl is in working dir
import matplotlib.pyplot as plt, pickle as pkl, pandas as pd
from matplotlib import cm
import numpy as np
from sys import argv
from scipy import stats

script, csv_dir, ifg, sasa = argv
#sasa is the cutoff in A**2

sasacsvpath=csv_dir+'%s_ifg_atom_density.csv'%ifg

colors = cm.Dark2.colors

sasacsv = pd.read_csv(sasacsvpath,index_col='iFG_count')

contactsdf, megalist = pkl.load(open('%s_interaction_types.pkl'%ifg,'rb'))

labels = ['# hbonds iFG makes', '# vdMs iFGs have vdw contacts with', '# vdMs iFGs hbond with','# vdMs iFGs have polar contacts with', '# atoms on vdMs that hbond with iFG', '# Ca atoms on vdMs that hbond with iFG', '# atoms on vdMs that make vdw contacts with iFG', '# atoms on vdMs that make polar contacts with iFG']  
#labels = ['# hbonds iFG makes', '# vdMs iFGs have vdw contacts with', '# vdMs iFGs hbond with','# vdMs iFGs have polar contacts with', '# atoms on vdMs that hbond with iFG', '# Ca atoms on vdMs that hbond with iFG', '# atoms on vdMs that make vdw contacts with iFG', '# atoms on vdMs that make polar contacts with iFG']  

########## ONLY GET ELEMENTS IN LIST THAT ARE BELOW SASA CUTOFF ##################
# get iFGs that are below cutoff 
# poss_ifgs has more ifgs than is in the contacts pkl file bc contacts pkl file 
# has repeats removed
poss_ifgs = sasacsv.index[sasacsv['iFG_sasa_CB_3A_probe'] < float(sasa)].tolist()
# for every ifg in contactsdf (that has repeats removed), if it's also in poss_ifgs, then keep 
# the index number
contactsdf_ifgs = list(set(contactsdf['iFG_count']))
contactsdf_ifgs = np.array(contactsdf_ifgs)


indices = np.isin(contactsdf_ifgs, np.array(poss_ifgs))
indices = list(np.where(indices==True)[0])
newlist = []
for ls in megalist: 
    newlist.append([ls[i] for i in indices])
megalist = newlist
[water_contacts_ifglevel , hb_contacts_vdmlevel , vdw_contacts_vdmlevel,
    polar_contacts_vdmlevel, hb_contacts_atomlevel, ca_hb_contacts_atomlevel, 
    vdw_contacts_atomlevel, polar_contacts_atomlevel, totalhbond ] = megalist 

########### only keep ifgs that pass cutoff so indices match ifg values
##########contactsdf_ifgs = contactsdf_ifgs[indices] 
########### get indices for ifgs that have 4 water contacts
##########unfiltered = [i for i, x in enumerate(water_contacts_ifglevel) if x==4]
########### use those indices to get ifg numbers
##########print(contactsdf_ifgs[unfiltered]) 

#f, axarr = plt.subplots(4)
#def plot(subplot,ix, each, cr):
#    hist = np.histogram(each, bins=int(max(each))+1)
#    x = np.arange(len(hist[0]))
#    x = [i+ix/5 for i in x] # offset so histograms don't overlap
#    axarr[subplot].bar(x,hist[0],color=colors[cr],label=labels[cr],lw=4,alpha=0.4,width=1)
#    axarr[subplot].legend()
## top first 2 subplots, for each ifg (vdmlevel)
#for ix, each in enumerate([totalhbond, vdw_contacts_vdmlevel]):
#    plot(0,ix,each,ix)
#for ix, each in enumerate([hb_contacts_vdmlevel, polar_contacts_vdmlevel]):
#    plot(1,ix,each,ix+2)
# bottom subplot, for each vdm (atomlevel)
#for ix, each in enumerate([hb_contacts_atomlevel, ca_hb_contacts_atomlevel]):
#    flat = [item for sublist in each for item in sublist]
#    plot(2,ix,flat,ix+4)
#for ix, each in enumerate([vdw_contacts_atomlevel, polar_contacts_atomlevel]):
#    flat = [item for sublist in each for item in sublist]
#    plot(3,ix,flat,ix+6)
#plt.suptitle('%s contacts - raw counts at %s A^2 cutoff' % (ifg,sasa))
#plt.show()


#### for scoring ####
# at the 0th # of bonds (so same as the bin edge). the LAST ONE is the sum of the 
# last 2 bin edges
scores = {'hbonds':{}, 'vdws':{}, 'polars':{}}
lbl = ['hbonds', 'vdws', 'polars']
for ix, ls in enumerate([totalhbond]):
#for ix, ls in enumerate([totalhbond, vdw_contacts_vdmlevel, polar_contacts_vdmlevel]):
    plt.figure()
    score_dict = scores[lbl[ix]]
    mean = np.mean(ls)
    print(mean)
    std = np.std(ls)
    for i in range(int(max(ls))+1):
        score_dict[i] = -(i-mean)/std
    logscores=[]
    for key, values in score_dict.items():
        logscores.append(values)
    
    # find N to know where to draw upper bound of n
    hist = np.histogram(ls, bins=int(max(ls)+1))
    N = np.cumsum(hist[0])
    N = sum([x<=0.95*sum(hist[0]) for x in N])
    
    
    plt.plot(logscores)
    plt.xlim(0, N)
    #plt.ylim(min(logscores), max(logscores))
    plt.ylim(min(logscores[:N+1]), max(logscores[:N+1]))
    plt.ylabel('score')
    plt.xlabel('# of interactions')
    plt.title('score for ideal # of hbonds %s makes (%s A^2 cutoff)'%(ifg, sasa))
    # er it craps out the same name for all the lists 
    #plt.show()

