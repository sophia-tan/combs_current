import pickle as pkl
import matplotlib.pyplot as plt
from cluster_for_sophia import *

D=pkl.load(open('./output_data/A755_755PHE_606ASP_pairwisematrix_planar_group.pkl','rb'))
adj_mat = make_adj_mat(D,.5)
clusters = greedy(adj_mat)
clus_sizes = [len(i) for i in clusters[0]]
print(np.mean(clus_sizes))
print(np.median(clus_sizes))
no_sing = [i for i in clus_sizes if i > 1]
print(np.mean(no_sing))
print(np.median(no_sing))
plt.hist(clus_sizes,bins=50)
plt.show()
