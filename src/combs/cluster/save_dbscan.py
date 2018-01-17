import scipy, scipy.spatial
from sklearn.cluster import DBSCAN
from scipy.cluster import hierarchy as hier
from scipy.cluster.hierarchy import fcluster
from scipy.spatial import distance as ssd
from sklearn.decomposition import PCA as sklearnPCA
from AtomsDictionary import *

def dbcluster(iFGcoords_array,pdbs,rmsd=1.1, min_core_num=6,print_clusters_to_file = False,outputdir = None):
    db = DBSCAN(eps=rmsd, min_samples=min_core_num).fit(iFGcoords_array)
    labels = db.labels_
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    ## to get cluster sizes
    clustersize_db=[]
    for i in range(0,n_clusters):
        clustersize_db.append(len(labels[labels==i]))
    sorted_labels = np.zeros(len(labels))
    indices=np.arange(0,len(labels))
    sortedlabelclustersize = sorted(zip(clustersize_db,range(0,n_clusters)),reverse=True)
    j=1
    for csize,i in sortedlabelclustersize:
        sorted_labels[labels==i] = j
        j+=1
    sorted_clustersize_db = sorted(clustersize_db,reverse=True)
    parsed_clusters = [[] for x in range(n_clusters)]
    j=0
    for csize,i in sorted(zip(clustersize_db,range(0,n_clusters)),reverse=True):
        for each in indices[labels==i]:
            parsed_clusters[j].append([pdbs[each][0],'cluster'+str(j)+'_'+str(csize)+'mems_'+pdbs[each][1]])       
        j+=1
    pair_dist = scipy.spatial.distance.pdist(iFGcoords_array)
    pair_dist_matrix = scipy.spatial.distance.squareform(pair_dist)
    # pair_dist_matrix is the rmsd matrix
    cluster_centroids = []
    for k in range(1,n_clusters+1):
        ind = np.argmin(np.mean(pair_dist_matrix[sorted_labels==k,:][:,sorted_labels==k],0))
        cluster_centroids.append([parsed_clusters[k-1][ind][0],'centroid_cluster'+str(k)+'_'+parsed_clusters[k-1][ind][1]])

    ##To print clusters to file:
    if print_clusters_to_file:
        for clus_num, each in enumerate(parsed_clusters):
            clus_num = str(clus_num+1)
            for item in each: 
                name = item[1].split('_')
                writePDB(outputdir+'/'+'cluster%s_'%clus_num+name[3],item[0])

    return n_clusters,sorted_clustersize_db,sorted_labels, parsed_clusters, cluster_centroids

def do_per_cluster_pca_analysis(parsed_clusters,iFGatomlist,originatomlist,vec1atomlist,vec2atomlist,stdcond):
    param_arrays = []
    pca_analysis =[]
    for cluster in parsed_clusters:
        param_array = make_6param_iFG_array(cluster,iFGatomlist,originatomlist,vec1atomlist,vec2atomlist)
        pca_components,pca_explained_var_ratio,pca_covar_matrix = do_pca(param_array,std_cond=stdcond)
        param_arrays.append(param_array)
        pca_analysis.append([pca_components,pca_explained_var_ratio,pca_covar_matrix])
    return param_arrays, pca_analysis

def do_db_clustering_and_pca(iFGcoords_array,pdbs,iFGatomlist,originatomlist,vec1atomlist,vec2atomlist,rmsdcutoff=1.1, mincorenum=6, print_clusters = False,outdir = None,stdcond=False):
    n_clusters,sorted_clustersize_db,sorted_labels, parsed_clusters, cluster_centroids = dbcluster(iFGcoords_array,pdbs,rmsd=rmsdcutoff, min_core_num=mincorenum,print_clusters_to_file = print_clusters,outputdir = outdir)
    
    param_arrays, pca_analysis = do_per_cluster_pca_analysis(parsed_clusters,iFGatomlist,originatomlist,vec1atomlist,vec2atomlist,stdcond)
    return n_clusters,sorted_clustersize_db,sorted_labels, parsed_clusters, cluster_centroids, param_arrays, pca_analysis 
    
# you can only run this after you generate the pkl file!
def count_vdm(AA, two_vdm=False, second_AA=None):
    ''' Find out how many vdms of this AA are in the whole dataset, and how many are clustered.
    Input: string of 3-letter AA code (ex: 'ASN')
    Returns: # vdms of that AA in dataset, # clustered, and ratio'''
    
    vdm_count = 0
    clustered_vdm_count = 0
    with open('/home/gpu/Nick/output/comb_database_carboxamide.csv') as infile:
        for line in infile:
            line = line.strip()
            split = line.split(',')
            if split[2] == AA:
                vdm_count += 1
    all_data = pic.load(open('data_%s.pkl' % AA,'rb'))
    param_arrays = all_data[7]
    for cluster in param_arrays:
        clustered_vdm_count += len(cluster)
    ratio = clustered_vdm_count / vdm_count

    return clustered_vdm_count, vdm_count, ratio

def make_rmsd_matrix(iFGcoords_array):
    pair_dist = scipy.spatial.distance.pdist(iFGcoords_array)
    pair_dist_matrix = scipy.spatial.distance.squareform(pair_dist)
    # pair_dist_matrix is the rmsd matrix
    return pair_dist_matrix

def plot_pc_space(param_arrays,title):
    #### to plot boring 2D scatterplot in PC space###
    fig, axarr = plt.subplots(5,10,sharex=True, sharey=True)
    fig_row = 0
    fig_column = 0
    for i in range(len(param_arrays)): # i is cluster number
        sklearn_pca = sklearnPCA(n_components=2)
        sklearn_transf = sklearn_pca.fit_transform(param_arrays[i])  
        x_values = sklearn_transf[:,0]
        y_values = sklearn_transf[:,1]
        axarr[fig_row,fig_column].scatter(x_values, y_values, s=0.5)
        axes = plt.gca()
        if fig_column < 9:
            fig_column += 1
        else:
            fig_column = 0
            fig_row += 1
    
    fig.text(0.5, 0.01, 'PC1', ha='center')
    fig.text(0.04, 0.5, 'PC2', va='center', rotation='vertical')
    plt.suptitle('PCA for carboxamide interactomers with ASP+GLU vandermers')
    plt.show()

def plot_original_parameters_space(param_arrays, pca_analysis):
    ## limit param_arrays
    param_arrays = param_arrays[:20]
    pca_analysis = pca_analysis[:20]
    radians = []
    for i, clus in enumerate(param_arrays):
        radians_clus = np.array([])
        for ix, c in enumerate(clus):
            p1 = c[3]/1.33
            p2 = c[4]/1.33
            p3 = c[5]/1.33
            array = np.array([c[0], c[1], c[2], p1, p2, p3])
            if ix == 0:
                radians_clus = array
            else:
                radians_clus = np.vstack((radians_clus, array))

        radians.append(radians_clus)
    param_arrays = np.copy(radians)

    sns.set_context("notebook", font_scale=0.8)
    params_dict = {}
    params_dict[0] = 'r'
    params_dict[1] = 'phi'
    params_dict[2] = 'theta'
    params_dict[3] = 'x'
    params_dict[4] = 'y'
    params_dict[5] = 'z'
    f, axarr = plt.subplots(4,7,figsize=(12,6.5))
    sns.despine(left=True)
    fig_row = 0
    fig_column = 0
    for clusnum, cluster in enumerate(pca_analysis):
        status = 0 # when status changes to 1, it means only top 1 pc vector. if status is 2, top 2 pca
        exp_var = 0 # will be updated
        explained_var = cluster[1]
        important_parameters = [] # integers that represent parameters 
        if explained_var[0] > 0.7: # only look at pc1
            weights = cluster[0][0] # for pc1 only
            indices = np.where(abs(weights)>0.6)[0]
            for ix in indices:
                important_parameters.append(ix)
            if len(indices) == 0: 
                indices = np.where(abs(weight)>0.5)[0]
                for ix in indices:
                    important_parameters.append(ix)
            status = 1
            exp_var = explained_var[0]
        else: # look at top 2 pc
            count = 0
            for weights in (cluster[0][0], cluster[0][1]):
                indices = np.where(abs(weights)>0.6)[0]
                for ix in indices:
                    important_parameters.append(ix)
                    count += 1
            if count==0:
                for weights in (cluster[0][0],cluster[0][1]):
                    indices = np.where(abs(weights)>0.5)[0]
                    for ix in indices:
                        important_parameters.append(ix)
            status = 2 
            exp_var = explained_var[0] + explained_var[1]
        important_parameters = list(set(important_parameters))
        if len(important_parameters) == 1:
            sns.distplot(param_arrays[clusnum][:,important_parameters[0]], ax=axarr[fig_row,fig_column])
            axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[0]]), verticalalignment='top', horizontalalignment='right',
                            transform=axarr[fig_row,fig_column].transAxes,
                            fontsize=10, weight='bold')
            if fig_column < 6:
                fig_column += 1
            else:
                fig_column = 0
                fig_row += 1
        if len(important_parameters) == 2:
            axarr[fig_row,fig_column].scatter(param_arrays[clusnum][:,important_parameters[0]], param_arrays[clusnum][:,important_parameters[1]], s=1.3)
            axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s vs %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[0]],params_dict[important_parameters[1]]), verticalalignment='top', horizontalalignment='right',
                    transform=axarr[fig_row,fig_column].transAxes,
                    fontsize=10, weight='bold')
    
            if fig_column < 6:
                fig_column += 1
            else:
                fig_column = 0
                fig_row += 1
        if len(important_parameters) == 3:
            axarr[fig_row,fig_column].scatter(param_arrays[clusnum][:,important_parameters[0]], param_arrays[clusnum][:,important_parameters[1]], s=1.3)
            axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s vs %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[0]],params_dict[important_parameters[1]]), verticalalignment='top', horizontalalignment='right',
                    transform=axarr[fig_row,fig_column].transAxes,
                    fontsize=10, weight='bold')
            if fig_column < 6:
                fig_column += 1
            else:
                fig_column = 0
                fig_row += 1
            axarr[fig_row,fig_column].scatter(param_arrays[clusnum][:,important_parameters[1]], param_arrays[clusnum][:,important_parameters[2]], s=1.3)
            axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s vs %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[1]],params_dict[important_parameters[2]]), verticalalignment='top', horizontalalignment='right',
                    transform=axarr[fig_row,fig_column].transAxes,
                    fontsize=10, weight='bold')
            if fig_column < 6:
                fig_column += 1
            else:
                fig_column = 0
                fig_row += 1
            axarr[fig_row,fig_column].scatter(param_arrays[clusnum][:,important_parameters[0]], param_arrays[clusnum][:,important_parameters[2]], s=1.3)
            axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s vs %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[0]],params_dict[important_parameters[2]]), verticalalignment='top', horizontalalignment='right',
                    transform=axarr[fig_row,fig_column].transAxes,
                    fontsize=10, weight='bold')
            if fig_column < 6:
                fig_column += 1
            else:
                fig_column = 0
                fig_row += 1
    plt.show()


def plot_contours_original_parameters_space(param_arrays, pca_analysis):
    ## limit param_arrays
    param_arrays = param_arrays[:2]
    pca_analysis = pca_analysis[:2]
    f, axarr = plt.subplots(1,2,figsize=(12,6.5))
    for clusnum, cluster in enumerate(pca_analysis):
        explained_var = cluster[1]
        important_parameters = [] # integers that represent parameters 
        if explained_var[0] > 0.7: # only look at pc1
            weights = cluster[0][0] # for pc1 only
            indices = np.where(abs(weights)>0.6)[0]
            for ix in indices:
                important_parameters.append(ix)
            if len(indices) == 0: 
                indices = np.where(abs(weight)>0.5)[0]
                for ix in indices:
                    important_parameters.append(ix)
        else: # look at top 2 pc
            count = 0
            for weights in (cluster[0][0], cluster[0][1]):
                indices = np.where(abs(weights)>0.6)[0]
                for ix in indices:
                    important_parameters.append(ix)
                    count += 1
            if count==0:
                for weights in (cluster[0][0],cluster[0][1]):
                    indices = np.where(abs(weights)>0.5)[0]
                    for ix in indices:
                        important_parameters.append(ix)
        important_parameters = list(set(important_parameters))
        
        if len(important_parameters) == 1:
            pass
        if len(important_parameters) == 2:
            counts,xbins,ybins=np.histogram2d(param_arrays[clusnum][:,important_parameters[0]],param_arrays[clusnum][:,important_parameters[1]])
            # get the avg of each xbin and ybin 
            x_vals = []
            y_vals = []
            for i in range(len(xbins)-1):
                x_vals.append((xbins[i]+xbins[i+1])/2)
                y_vals.append((ybins[i]+ybins[i+1])/2)
            
            CS = axarr[clusnum].contour(x_vals,y_vals,counts)
            #axarr[clusnum].clabel(CS, inline=1, fontsize=10)
            axarr[clusnum].scatter(param_arrays[clusnum][:,important_parameters[0]], param_arrays[clusnum][:,important_parameters[1]], s=1.3)
            #plt.clabel(CS, inline=1, fontsize=10)
            #plt.title('ASN/GLN cluster 1')
            #axarr[fig_row,fig_column].scatter(param_arrays[clusnum][:,important_parameters[0]], param_arrays[clusnum][:,important_parameters[1]], s=1.3)
            #axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s vs %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[0]],params_dict[important_parameters[1]]), verticalalignment='top', horizontalalignment='right',
            #        transform=axarr[fig_row,fig_column].transAxes,
            #        fontsize=10, weight='bold')
    plt.show()
'''
        if len(important_parameters) == 3:
            axarr[fig_row,fig_column].scatter(param_arrays[clusnum][:,important_parameters[0]], param_arrays[clusnum][:,important_parameters[1]], s=1.3)
            axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s vs %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[0]],params_dict[important_parameters[1]]), verticalalignment='top', horizontalalignment='right',
                    transform=axarr[fig_row,fig_column].transAxes,
                    fontsize=10, weight='bold')
            if fig_column < 6:
                fig_column += 1
            else:
                fig_column = 0
                fig_row += 1
            axarr[fig_row,fig_column].scatter(param_arrays[clusnum][:,important_parameters[1]], param_arrays[clusnum][:,important_parameters[2]], s=1.3)
            axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s vs %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[1]],params_dict[important_parameters[2]]), verticalalignment='top', horizontalalignment='right',
                    transform=axarr[fig_row,fig_column].transAxes,
                    fontsize=10, weight='bold')
            if fig_column < 6:
                fig_column += 1
            else:
                fig_column = 0
                fig_row += 1
            axarr[fig_row,fig_column].scatter(param_arrays[clusnum][:,important_parameters[0]], param_arrays[clusnum][:,important_parameters[2]], s=1.3)
            axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s vs %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[0]],params_dict[important_parameters[2]]), verticalalignment='top', horizontalalignment='right',
                    transform=axarr[fig_row,fig_column].transAxes,
                    fontsize=10, weight='bold')
            if fig_column < 6:
                fig_column += 1
            else:
                fig_column = 0
                fig_row += 1







    
            plt.tight_layout()
            plt.show()

'''


# this is for carboxamide only
from PCA_functions import *
from AtomsDictionary import * 

    
#################################
##
###do pca and pickle
##pca_data = do_per_cluster_pca_analysis(sorted_clusters, iFG_sel_lists, iFG_origin_list, iFG_x_list, iFG_inplane_list,stdcond=False)
##pic.dump(pca_data,open('pca_data_ASP_GLU.pkl','wb'))
#
######################################################################################################
## load pca data, just use param_arrays pca_data[0] to do PCA again (sklearn_transf) and NOT pca_analysis pca_data[1]
#pca_data = pic.load(open('pca_data_ASP_GLU.pkl','rb'))
#param_arrays, pca_analysis = pca_data
#
### limit param_arrays
#param_arrays = param_arrays[:18]
#pca_analysis = pca_analysis[:18]
#
######################################################################################################
##plot_original_parameters_space(param_arrays, pca_analysis)
##plot_contours_original_parameters_space(param_arrays, pca_analysis)
######################################################################################################
#
#
#ex to run!!
##iFGcoords, pdbs = parsePDBs_and_make_iFG_array('/home/gpu/Nick/output/output_carboxamide','/home/gpu/Nick/output/comb_database_carboxamide.csv',['AA','iFG_iRes_sc_N_O_valence'],['ASP',0],vdm_origin, vdm_x, vdm_inplane, iFG_sel_lists, len(iFG_sel_lists[0])) 
#'''
