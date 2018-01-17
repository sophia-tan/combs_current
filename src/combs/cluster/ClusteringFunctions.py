__all__ = ['hierarchical_cluster', 'do_dbscan_for_bb_and_sc']

from prody import *
import numpy as np
import sys
import os, pickle as pkl
import shutil
from ..apps.constants import interactamer_atoms, residue_sc_names
import scipy
from scipy.cluster import hierarchy as hier
from scipy.cluster.hierarchy import fcluster
import scipy, scipy.spatial
from sklearn.cluster import DBSCAN
import math

def concat_hier_pca_inputs(list_iFGcoords, list_pdbs):
    '''Concatenates inputs for hierarchical clustering and PCA. Use if you're analyzing >1 VDM AA 
    (ex: ASN+GLN or ASP+GLU)'''
    multiple_iFGcoords = np.vstack(list_iFGcoords)
    if len(list_pdbs)==1:
        multiple_pdbs = list_pdbs[0]
    elif len(list_pdbs)==2:
        multiple_pdbs=list_pdbs[0]+list_pdbs[1]
    else:
        multiple_pdbs = []
        for ix in range(len(list_pdbs)):
            multiple_pdbs += list_pdbs[ix]
    return multiple_iFGcoords, multiple_pdbs

def getpdbpath(pkldata,ifg,pkldir,BBorSC):
    '''helper function. gets pdbpath from pickle file'''
    ls = []
    for item in pkldata:
        ifgcount,vdmcount = item[0],item[1]
        resn = item[9]
        resngroup = list(interactamer_atoms[resn].keys())[0]
        path = 'iFG_%s_vdM_%s_%s_oriented.pdb.gz'\
            % (ifgcount,vdmcount,ifg)
        if BBorSC == 'SC':
            pdbdir = pkldir+'pdbs/%s/%s/'%(resn,resngroup)
        elif BBorSC == 'N_CA' or BBorSC=='C_O':
            pdbdir = pkldir+'pdbs/%s/'%(resn) # doesn't have resngroup folder
        path = pdbdir+path
        ls.append([[resn, ifgcount, vdmcount], path])
    return ls

def do_dbscan_for_bb_and_sc(combsoutputdir, ifg, vdmresidues, print_clus, typ):
    '''print_clus needs to be booleani. typ is BBorSC'''
    
    pkldir = combsoutputdir+ifg+'/clusters/%s/'%typ
    if typ == 'SC':
        vdmresidues = vdmresidues.split(',')
    
    # pool all residues if bb
    elif typ=='N_CA' or typ=='C_O':
        vdmresidues = []
        for res,val in residue_sc_names.items():
            try:
                pkl.load(open(pkldir+'pickle/%s_rel_vdms.pickle'%res,'rb')) 
                vdmresidues.append(res)
            except:
                pass

    # make outputdir 
    outdir = combsoutputdir+ifg+'/clusters/%s/clustered/'%typ
    try:
        os.mkdir(outdir)
    except:
        pass
    # load ifgcoords and ifgcount/vdmcount info
    vdmnames = '+'.join(vdmresidues) # for naming pickle file
    pkls = [pkl.load(open(pkldir+'pickle/%s_rel_vdms.pickle'%vdm,'rb')) for vdm in vdmresidues]
    ifgcoords = [x[:,6] for x in pkls]
    ifgcoords = [np.array([y.flatten() for y in x]) for x in ifgcoords]
    pdbpaths = [getpdbpath(x,ifg,pkldir,typ) for x in pkls]
    ifgcoords, pdbs = concat_hier_pca_inputs(ifgcoords,pdbpaths)
    rms = (len(ifgcoords[0])/3-1)*0.25 # rmsd cutoff is 0.75 for len(4), 0.5 for len(3)
    # adjust rms for # of samples

    rms = rms - math.log(len(ifgcoords)/10000,10)*.25
    print(rms, 'rms', len(ifgcoords), ifg, vdmresidues, typ)
    if print_clus==False:
        parsed_clusters, centroids = dbcluster(ifgcoords, pdbs, rmsd=rms, min_core_num=1,BBorSC=typ,print_clusters=False,outputdir=None,vdmnames=None)
    elif print_clus==True:
        parsed_clusters, centroids = dbcluster(ifgcoords, pdbs, rmsd=rms, min_core_num=1,BBorSC=typ,print_clusters=True, outputdir=outdir, vdmnames=vdmnames)
        if typ == 'SC':
            pklfile = open(outdir+'clusterinfo_%s_%s_%s.pkl'%(ifg,vdmnames,typ), 'wb')
        else:
            pklfile = open(outdir+'clusterinfo_%s_%s.pkl'%(ifg,typ), 'wb')
        pkl.dump([parsed_clusters,centroids], pklfile)
        # format of parsed_clusters is list of clusters. in each cluster, each member has [labels, filename, ifgcoords]
        # format of centroids is list of centroids. each centroid is [labels, filename]

def make_rmsd_matrix(iFGcoords_array):
    pair_dist = scipy.spatial.distance.pdist(iFGcoords_array)
    pair_dist_matrix = scipy.spatial.distance.squareform(pair_dist)
    # pair_dist_matrix is the rmsd matrix
    return pair_dist_matrix
    
def hierarchical_cluster(iFGcoords, pdbs_of_iFG_coords, max_d=0.6, print_clusters = False, outputdir = None, vdmnames=None,BBorSC=None):
    '''Performs hierarchical clustering on a set of iFGcoords array
    Inputs: max_d (hierarchical clustering cutoff), ifgcoordsarray, pdbs. Option to print the 
            clusters to a directory to view what the clusters look like.
    Returns: sorted_clusters (list where each element is a list of all the pdbs in that cluster)
    '''

    rmsd_mat = make_rmsd_matrix(iFGcoords)
    '''
    Z = hier.linkage(iFGcoords,method='single')
    #Z = hier.linkage(scipy.spatial.distance.squareform(rmsd_mat),method='single')
    clusters = fcluster(Z, max_d, criterion='distance') # retrieve clusters
    #print(clusters)
    n_clus = max(clusters)
    unsorted_clusters = []
    for n in range(1,n_clus+1):
        members = [np.where(clusters == n)][0][0]
        members = [pdbs_of_iFG_coords[m] for m in members]
        unsorted_clusters.append(members)
    clustersizes = np.array([np.size(k,0) for k in unsorted_clusters])
    # sort by cluster sizes (largest cluster to smallest cluster)
    sorted_clusters = sorted(unsorted_clusters,key=len,reverse=True)
    
    if print_clusters == True:
        for i, members in enumerate(sorted_clusters):
            clus_num = str(i+1)
            for r in members:
                orig = r[1]
                name = '_'+r[1].split('/')[-1] 
                if BBorSC=='SC':
                    outfile = outputdir+'%s_cluster%s'%(vdmnames,clus_num+name)
                elif BBorSC=='N_CA' or BBorSC=='C_O':
                    outfile = outputdir+'cluster%s'%(clus_num+name)
                shutil.copy(orig, outfile)
    return sorted_clusters
    '''


def dbcluster(iFGcoords_array,pdbs,BBorSC,print_clusters,vdmnames,outputdir,rmsd,min_core_num):
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
            parsed_clusters[j].append([pdbs[each][0], pdbs[each][1], iFGcoords_array[each]])
        j+=1
    cluster_centroids = []
    try:
        pair_dist = scipy.spatial.distance.pdist(iFGcoords_array)
        pair_dist_matrix = scipy.spatial.distance.squareform(pair_dist)
        # pair_dist_matrix is the rmsd matrix
        for k in range(1,n_clusters+1):
            ind = np.argmin(np.mean(pair_dist_matrix[sorted_labels==k,:][:,sorted_labels==k],0))
            cluster_centroids.append([parsed_clusters[k-1][ind][0],'centroid_cluster'+str(k)+'_'+parsed_clusters[k-1][ind][1]])
    except:
        pass # clusters too large to compute pdist to get centroids

    ##To print clusters to file:
    if print_clusters:
        for i, members in enumerate(parsed_clusters):
            clus_num = str(i+1)
            for r in members:
                orig = r[1]
                name = '_'+r[1].split('/')[-1] 
                if r[0] in [x[0] for x in cluster_centroids]:
                    head = clus_num+'_centroid'
                else:
                    head = clus_num
                resn = r[0][0]
                outfile = outputdir+'%s_cluster%s'%(resn,head+name)
                shutil.copy(orig, outfile)
    # format of parsed_clusters is list of clusters. in each cluster, each member has [labels, filename, ifgcoords]
    # format of cluster_centroids is list of centroids. each centroid is [labels, filename]
    return parsed_clusters, cluster_centroids

#def do_per_cluster_pca_analysis(parsed_clusters,iFGatomlist,originatomlist,vec1atomlist,vec2atomlist,stdcond):
#    param_arrays = []
#    pca_analysis =[]
#    for cluster in parsed_clusters:
#        param_array = make_6param_iFG_array(cluster,iFGatomlist,originatomlist,vec1atomlist,vec2atomlist)
#        pca_components,pca_explained_var_ratio,pca_covar_matrix = do_pca(param_array,std_cond=stdcond)
#        param_arrays.append(param_array)
#        pca_analysis.append([pca_components,pca_explained_var_ratio,pca_covar_matrix])
#    return param_arrays, pca_analysis
#
#def do_db_clustering_and_pca(iFGcoords_array,pdbs,iFGatomlist,originatomlist,vec1atomlist,vec2atomlist,rmsdcutoff=1.1, mincorenum=6, print_clusters = False,outdir = None,stdcond=False):
#    n_clusters,sorted_clustersize_db,sorted_labels, parsed_clusters, cluster_centroids = dbcluster(iFGcoords_array,pdbs,rmsd=rmsdcutoff, min_core_num=mincorenum,print_clusters_to_file = print_clusters,outputdir = outdir)
#    
#    param_arrays, pca_analysis = do_per_cluster_pca_analysis(parsed_clusters,iFGatomlist,originatomlist,vec1atomlist,vec2atomlist,stdcond)
#    return n_clusters,sorted_clustersize_db,sorted_labels, parsed_clusters, cluster_centroids, param_arrays, pca_analysis 
#    
## you can only run this after you generate the pkl file!
#def count_vdm(AA, two_vdm=False, second_AA=None):
#    ''' Find out how many vdms of this AA are in the whole dataset, and how many are clustered.
#    Input: string of 3-letter AA code (ex: 'ASN')
#    Returns: # vdms of that AA in dataset, # clustered, and ratio'''
#    
#    vdm_count = 0
#    clustered_vdm_count = 0
#    with open('/home/gpu/Nick/output/comb_database_carboxamide.csv') as infile:
#        for line in infile:
#            line = line.strip()
#            split = line.split(',')
#            if split[2] == AA:
#                vdm_count += 1
#    all_data = pic.load(open('data_%s.pkl' % AA,'rb'))
#    param_arrays = all_data[7]
#    for cluster in param_arrays:
#        clustered_vdm_count += len(cluster)
#    ratio = clustered_vdm_count / vdm_count
#
#    return clustered_vdm_count, vdm_count, ratio
#
#def make_rmsd_matrix(iFGcoords_array):
#    pair_dist = scipy.spatial.distance.pdist(iFGcoords_array)
#    pair_dist_matrix = scipy.spatial.distance.squareform(pair_dist)
#    # pair_dist_matrix is the rmsd matrix
#    return pair_dist_matrix
#
#def plot_pc_space(param_arrays,title):
#    #### to plot boring 2D scatterplot in PC space###
#    fig, axarr = plt.subplots(5,10,sharex=True, sharey=True)
#    fig_row = 0
#    fig_column = 0
#    for i in range(len(param_arrays)): # i is cluster number
#        sklearn_pca = sklearnPCA(n_components=2)
#        sklearn_transf = sklearn_pca.fit_transform(param_arrays[i])  
#        x_values = sklearn_transf[:,0]
#        y_values = sklearn_transf[:,1]
#        axarr[fig_row,fig_column].scatter(x_values, y_values, s=0.5)
#        axes = plt.gca()
#        if fig_column < 9:
#            fig_column += 1
#        else:
#            fig_column = 0
#            fig_row += 1
#    
#    fig.text(0.5, 0.01, 'PC1', ha='center')
#    fig.text(0.04, 0.5, 'PC2', va='center', rotation='vertical')
#    plt.suptitle('PCA for carboxamide interactomers with ASP+GLU vandermers')
#    plt.show()
#
#def plot_original_parameters_space(param_arrays, pca_analysis):
#    ## limit param_arrays
#    param_arrays = param_arrays[:20]
#    pca_analysis = pca_analysis[:20]
#    radians = []
#    for i, clus in enumerate(param_arrays):
#        radians_clus = np.array([])
#        for ix, c in enumerate(clus):
#            p1 = c[3]/1.33
#            p2 = c[4]/1.33
#            p3 = c[5]/1.33
#            array = np.array([c[0], c[1], c[2], p1, p2, p3])
#            if ix == 0:
#                radians_clus = array
#            else:
#                radians_clus = np.vstack((radians_clus, array))
#
#        radians.append(radians_clus)
#    param_arrays = np.copy(radians)
#
#    sns.set_context("notebook", font_scale=0.8)
#    params_dict = {}
#    params_dict[0] = 'r'
#    params_dict[1] = 'phi'
#    params_dict[2] = 'theta'
#    params_dict[3] = 'x'
#    params_dict[4] = 'y'
#    params_dict[5] = 'z'
#    f, axarr = plt.subplots(4,7,figsize=(12,6.5))
#    sns.despine(left=True)
#    fig_row = 0
#    fig_column = 0
#    for clusnum, cluster in enumerate(pca_analysis):
#        status = 0 # when status changes to 1, it means only top 1 pc vector. if status is 2, top 2 pca
#        exp_var = 0 # will be updated
#        explained_var = cluster[1]
#        important_parameters = [] # integers that represent parameters 
#        if explained_var[0] > 0.7: # only look at pc1
#            weights = cluster[0][0] # for pc1 only
#            indices = np.where(abs(weights)>0.6)[0]
#            for ix in indices:
#                important_parameters.append(ix)
#            if len(indices) == 0: 
#                indices = np.where(abs(weight)>0.5)[0]
#                for ix in indices:
#                    important_parameters.append(ix)
#            status = 1
#            exp_var = explained_var[0]
#        else: # look at top 2 pc
#            count = 0
#            for weights in (cluster[0][0], cluster[0][1]):
#                indices = np.where(abs(weights)>0.6)[0]
#                for ix in indices:
#                    important_parameters.append(ix)
#                    count += 1
#            if count==0:
#                for weights in (cluster[0][0],cluster[0][1]):
#                    indices = np.where(abs(weights)>0.5)[0]
#                    for ix in indices:
#                        important_parameters.append(ix)
#            status = 2 
#            exp_var = explained_var[0] + explained_var[1]
#        important_parameters = list(set(important_parameters))
#        if len(important_parameters) == 1:
#            sns.distplot(param_arrays[clusnum][:,important_parameters[0]], ax=axarr[fig_row,fig_column])
#            axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[0]]), verticalalignment='top', horizontalalignment='right',
#                            transform=axarr[fig_row,fig_column].transAxes,
#                            fontsize=10, weight='bold')
#            if fig_column < 6:
#                fig_column += 1
#            else:
#                fig_column = 0
#                fig_row += 1
#        if len(important_parameters) == 2:
#            axarr[fig_row,fig_column].scatter(param_arrays[clusnum][:,important_parameters[0]], param_arrays[clusnum][:,important_parameters[1]], s=1.3)
#            axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s vs %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[0]],params_dict[important_parameters[1]]), verticalalignment='top', horizontalalignment='right',
#                    transform=axarr[fig_row,fig_column].transAxes,
#                    fontsize=10, weight='bold')
#    
#            if fig_column < 6:
#                fig_column += 1
#            else:
#                fig_column = 0
#                fig_row += 1
#        if len(important_parameters) == 3:
#            axarr[fig_row,fig_column].scatter(param_arrays[clusnum][:,important_parameters[0]], param_arrays[clusnum][:,important_parameters[1]], s=1.3)
#            axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s vs %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[0]],params_dict[important_parameters[1]]), verticalalignment='top', horizontalalignment='right',
#                    transform=axarr[fig_row,fig_column].transAxes,
#                    fontsize=10, weight='bold')
#            if fig_column < 6:
#                fig_column += 1
#            else:
#                fig_column = 0
#                fig_row += 1
#            axarr[fig_row,fig_column].scatter(param_arrays[clusnum][:,important_parameters[1]], param_arrays[clusnum][:,important_parameters[2]], s=1.3)
#            axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s vs %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[1]],params_dict[important_parameters[2]]), verticalalignment='top', horizontalalignment='right',
#                    transform=axarr[fig_row,fig_column].transAxes,
#                    fontsize=10, weight='bold')
#            if fig_column < 6:
#                fig_column += 1
#            else:
#                fig_column = 0
#                fig_row += 1
#            axarr[fig_row,fig_column].scatter(param_arrays[clusnum][:,important_parameters[0]], param_arrays[clusnum][:,important_parameters[2]], s=1.3)
#            axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s vs %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[0]],params_dict[important_parameters[2]]), verticalalignment='top', horizontalalignment='right',
#                    transform=axarr[fig_row,fig_column].transAxes,
#                    fontsize=10, weight='bold')
#            if fig_column < 6:
#                fig_column += 1
#            else:
#                fig_column = 0
#                fig_row += 1
#    plt.show()
#
#
#def plot_contours_original_parameters_space(param_arrays, pca_analysis):
#    ## limit param_arrays
#    param_arrays = param_arrays[:2]
#    pca_analysis = pca_analysis[:2]
#    f, axarr = plt.subplots(1,2,figsize=(12,6.5))
#    for clusnum, cluster in enumerate(pca_analysis):
#        explained_var = cluster[1]
#        important_parameters = [] # integers that represent parameters 
#        if explained_var[0] > 0.7: # only look at pc1
#            weights = cluster[0][0] # for pc1 only
#            indices = np.where(abs(weights)>0.6)[0]
#            for ix in indices:
#                important_parameters.append(ix)
#            if len(indices) == 0: 
#                indices = np.where(abs(weight)>0.5)[0]
#                for ix in indices:
#                    important_parameters.append(ix)
#        else: # look at top 2 pc
#            count = 0
#            for weights in (cluster[0][0], cluster[0][1]):
#                indices = np.where(abs(weights)>0.6)[0]
#                for ix in indices:
#                    important_parameters.append(ix)
#                    count += 1
#            if count==0:
#                for weights in (cluster[0][0],cluster[0][1]):
#                    indices = np.where(abs(weights)>0.5)[0]
#                    for ix in indices:
#                        important_parameters.append(ix)
#        important_parameters = list(set(important_parameters))
#        
#        if len(important_parameters) == 1:
#            pass
#        if len(important_parameters) == 2:
#            counts,xbins,ybins=np.histogram2d(param_arrays[clusnum][:,important_parameters[0]],param_arrays[clusnum][:,important_parameters[1]])
#            # get the avg of each xbin and ybin 
#            x_vals = []
#            y_vals = []
#            for i in range(len(xbins)-1):
#                x_vals.append((xbins[i]+xbins[i+1])/2)
#                y_vals.append((ybins[i]+ybins[i+1])/2)
#            
#            CS = axarr[clusnum].contour(x_vals,y_vals,counts)
#            #axarr[clusnum].clabel(CS, inline=1, fontsize=10)
#            axarr[clusnum].scatter(param_arrays[clusnum][:,important_parameters[0]], param_arrays[clusnum][:,important_parameters[1]], s=1.3)
#            #plt.clabel(CS, inline=1, fontsize=10)
#            #plt.title('ASN/GLN cluster 1')
#            #axarr[fig_row,fig_column].scatter(param_arrays[clusnum][:,important_parameters[0]], param_arrays[clusnum][:,important_parameters[1]], s=1.3)
#            #axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s vs %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[0]],params_dict[important_parameters[1]]), verticalalignment='top', horizontalalignment='right',
#            #        transform=axarr[fig_row,fig_column].transAxes,
#            #        fontsize=10, weight='bold')
#    plt.show()
#'''
#        if len(important_parameters) == 3:
#            axarr[fig_row,fig_column].scatter(param_arrays[clusnum][:,important_parameters[0]], param_arrays[clusnum][:,important_parameters[1]], s=1.3)
#            axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s vs %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[0]],params_dict[important_parameters[1]]), verticalalignment='top', horizontalalignment='right',
#                    transform=axarr[fig_row,fig_column].transAxes,
#                    fontsize=10, weight='bold')
#            if fig_column < 6:
#                fig_column += 1
#            else:
#                fig_column = 0
#                fig_row += 1
#            axarr[fig_row,fig_column].scatter(param_arrays[clusnum][:,important_parameters[1]], param_arrays[clusnum][:,important_parameters[2]], s=1.3)
#            axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s vs %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[1]],params_dict[important_parameters[2]]), verticalalignment='top', horizontalalignment='right',
#                    transform=axarr[fig_row,fig_column].transAxes,
#                    fontsize=10, weight='bold')
#            if fig_column < 6:
#                fig_column += 1
#            else:
#                fig_column = 0
#                fig_row += 1
#            axarr[fig_row,fig_column].scatter(param_arrays[clusnum][:,important_parameters[0]], param_arrays[clusnum][:,important_parameters[2]], s=1.3)
#            axarr[fig_row,fig_column].text(0.95, 0.95,'%s, %s, %s, %s vs %s' % (clusnum+1, status,np.around(exp_var, decimals=1),params_dict[important_parameters[0]],params_dict[important_parameters[2]]), verticalalignment='top', horizontalalignment='right',
#                    transform=axarr[fig_row,fig_column].transAxes,
#                    fontsize=10, weight='bold')
#            if fig_column < 6:
#                fig_column += 1
#            else:
#                fig_column = 0
#                fig_row += 1
#
#
#
#
#
#
#
#    
#            plt.tight_layout()
#            plt.show()
#
#'''
#
#
