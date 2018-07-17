import sys, os
from scipy import stats
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
from residues_integrin import *
import pickle as pkl, numpy as np, pandas as pd
import matplotlib.pyplot as plt, seaborn as sns
sns.set_style("white")

### PLOTTING DETAILS ###
ax = plt.subplot(111)

def best_fit_slope(xs,ys):
    xs,ys=np.array(xs),np.array(ys)
    m = (((np.mean(xs)*np.mean(ys)) - np.mean(xs*ys)) / \
             ((np.mean(xs)*np.mean(xs)) - np.mean(xs*xs)))
    b = np.mean(ys) - m*np.mean(xs)
    return m,b

ind = 0
    
AAs=['R671', 'I673', 'N753','F755','S758','V760', 'E785','H787','R900','L959', \
     'D552', 'Y594', 'T603','H626','K658', 'V664', 'E534']
activation = np.array([0.54,0.23,0.37,0.4,0.64,0.83,0.47,0.17,0.31,0.05, \
    .12, .44, .399, .10, .20, .42, 0.76])
sortedactivation = sorted(activation)
new_order = [np.where(activation==i)[0][0] for ix,i in enumerate(sortedactivation)]
AAs = [AAs[ix] for ix in new_order]
activation = [activation[ix] for ix in new_order]

### ANALYZING DETAILS ###
def rmsd_filter(method,cutoff,score,nonmembrane,buried):
    print('--------------------------')
    print('method, nonmembrane, buried, cutoff, score =',
        method, nonmembrane, buried, cutoff, score)
    
    norm_by_num_vdms = []
    norm_by_direct_vdms = []
    norm_by_avg_num_NNs = []
    norm_by_avg_num_NNs_nosing = []
    norm_by_med_num_NNs = []
    norm_by_med_num_NNs_nosing = []

    for targetres, resn in integrin_res.items():
    
        #print('----------------')
        #print('Target Res: ', targetres)
        rawcounts = [] # num obs
        num_vdms = [] # option to normalize by
        num_clustered = [] # option to normalize by
        avgclus = [] # option to normalize by
        avgclus_wo_sing = [] # option to normalize by
        medclus = [] # option to normalize by
        medclus_wo_sing = [] # option to normalize by
        for pklf in os.listdir('./output_data/'):
            if pklf.startswith(targetres[:4]) and pklf.endswith(
                'matches_{}_nonmembrane_{}_buried_{}_{}.pkl'.format(method,
                nonmembrane, buried, cutoff )):
                matches = pkl.load(open('./output_data/'+pklf,'rb')) 
                [rmsd, resi, ifgresn, ifgresi, vdmresn, vdmresi,
                        num_nn, norm_metrics] = matches
                num_all_vdms, num_direct, [avgsize, avgsize_no_sing,
                    medsize, medsize_no_sing] = norm_metrics

                rawcounts.append(num_nn+1)
                num_vdms.append(num_all_vdms + 1)
                num_clustered.append(num_direct + 1) 
                avgclus.append(avgsize + 1)
                avgclus_wo_sing.append(avgsize_no_sing + 1)
                medclus.append(medsize + 1)
                medclus_wo_sing.append(medsize_no_sing + 1)

                #print('Interacting residue: ', pklf.split('_')[2])
                #print(num_nn, num_all_vdms, num_direct)
        if rawcounts == []:
            rawcounts = [1]
            num_vdms = [10000]
            num_clustered = [10000]
            for ls in [avgclus, avgclus_wo_sing, medclus, medclus_wo_sing]:
                ls.append(10) 

        ######################################################################
        # Diff methods of calculating score for list of observed counts 'a' 
        # and list of expected counts 'b'
        ######################################################################
        def list_to_ndarray(a, b):
            return [np.array(a), np.array(b)]
        
        def calc(a,b): 
            a, b = list_to_ndarray(a, b)
            # all variations of just rawcounts, no normalization
            if score == 'a':
                return sum(a)
            if score == 'b':
                return sum(np.log10(a))
            if score == 'c':
                return np.log10(sum(a))

            # all variations of arrays 'a' and 'b' treated equally
            if score == 'd':
                return sum(a)/sum(b)
            if score == 'e':
                return sum(a/b)
            if score == 'f':
                return np.log10(sum(a/b))
            if score == 'g':
                return np.log10(sum(a)/sum(b))
            if score == 'h':
                return sum(np.log10(a/b))
            if score == 'i':
                return sum(np.log10(a)/np.log10(b))
            if score == 'j':
                return np.log10(sum(a))/np.log10(sum(b))
            if score == 'k':
                return sum(np.log10(a))/sum(np.log10(b))

            # all variations with log of only array 'a' and not 'b'
            if score == 'l':
                return sum(np.log10(a)/b)
            if score == 'm':
                return np.log10(sum(a))/sum(b)
            if score == 'n':
                return sum(np.log10(a))/sum(b)
            
            # all variations with log of only array 'b' and not 'a'
            if score == 'o':
                return sum(a/np.log10(b))
            if score == 'p':
                return sum(a)/np.log10(sum(b))
            if score == 'q':
                return sum(a)/sum(np.log10(b))

        #print(rawcounts)
        #print(avgclus)
        #print('rawcounts')
        #print(rawcounts)
        #print('num vdms')
        #print(num_vdms)
        norm_by_num_vdms.append(calc(rawcounts, num_vdms))
        norm_by_direct_vdms.append(calc(rawcounts, num_clustered))
        norm_by_avg_num_NNs.append(calc(rawcounts, avgclus))
        norm_by_avg_num_NNs_nosing.append(calc(rawcounts, avgclus_wo_sing))
        norm_by_med_num_NNs.append(calc(rawcounts, medclus))
        norm_by_med_num_NNs_nosing.append(calc(rawcounts, medclus_wo_sing))

    #return [norm_by_num_vdms, norm_by_direct_vdms, norm_by_avg_num_NNs,
    #    norm_by_avg_num_NNs_nosing, norm_by_med_num_NNs, norm_by_med_num_NNs_nosing]
    return [norm_by_num_vdms]
    
def get_correlation(megalist):
    AAs=['R671', 'I673', 'N753','F755','S758','V760', 'E785','H787','R900','L959', \
         'D552', 'Y594', 'T603','H626','K658', 'V664', 'E534']
    activation = np.array([0.54,0.23,0.37,0.4,0.64,0.83,0.47,0.17,0.31,0.05, \
        .12, .44, .399, .10, .20, .42, 0.76])
    #AAs=['R671', 'I673', 'N753','F755','S758','V760', 'E785','H787','R900','L959', \
    #     'D552', 'Y594', 'T603','H626','K658', 'V664']
    #activation = np.array([0.54,0.23,0.37,0.4,0.64,0.83,0.47,0.17,0.31,0.05, \
    #    .12, .44, .399, .10, .20, .42])

    ''' ughhhhhh taking out the funky serine'''
    activation = [activation[x] for x in range(len(activation)) if x != 4] # delete 
    mega_calc = []
    for origcalc in megalist:
        ''' ughhhhhh taking out the funky serine'''
        calc = [origcalc[x] for x in range(len(origcalc)) if x != 4] # delete
        spearcorr = stats.spearmanr(activation,calc)
        pearscorr = stats.pearsonr(activation,calc)
        mega_calc.append([pearscorr, origcalc])
    return mega_calc

#for method in ['sc_only_ifg']:
for method in ['planar_group_no_bb']:
    for nonmembrane in [False]:
        for buried in [False]:
            for score in ['j']:
            #for score in ['e','f','g','h','i','j','k','l','m','n','o','p','q']:
                for cut in [.3]:
                #for cut in [.3,.4,.5,.6,.7,.8,.9,1]:
                    megalist = rmsd_filter(method,cut,score,nonmembrane,buried)
                    
                    mega_calc = get_correlation(megalist)
                    for element in mega_calc:
                        pearscorr, calc = element
                        if pearscorr[0] > .7:
                            # plotting details

                            calc = [calc[ix] for ix in new_order]
                            
                            sns.despine(right=True)
                            
                            print(pearscorr)
                            ax.scatter(calc,activation,marker='o',edgecolors='gray',facecolors='c',
                                       lw=1)
                            ax.scatter(calc[14],activation[14],marker='o',edgecolors='black',facecolors='r',
                                lw=1) # overwriting with diff color
                    
                    
                            m,b=best_fit_slope(calc,activation)
                            regression_x = [x for x in calc]
                            regression_line = [(m*x)+b for x in calc]
                    
                            ax.plot(regression_x,regression_line)
                        
                            def text(label,x,y):
                                xytext=(5,-4.5)
                                if label=='H626':
                                    xytext=(5,-9)
                                ax.annotate(
                                    label,xy=(x,y),
                                    xytext=xytext,textcoords='offset points',size=12)
                            
                            for label, x, y in zip(AAs, calc,activation):
                                text(label,x,y)
                            
                            ind += 1
plt.title('geometric matches normalized by # vdms')
plt.xlabel('interaction score')
plt.ylabel('experimental activation index')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.show()
                    #    
                    #    plt.title('normalized by # directly interacting vdms, r^2=0.18')
                    #    plt.tight_layout()
                    #    plt.show()
