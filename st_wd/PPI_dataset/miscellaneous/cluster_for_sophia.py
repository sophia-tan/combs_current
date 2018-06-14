import numpy as np
from scipy.sparse import csr_matrix
from numba import jit


#Order of operations:
# 1. First make coordinate array of all coords *X* to be clustered.
#       Do this using prody.  Make sure atoms are printed in same order.
# 2. Coordinate array X is input to make_pairwise_rmsd_mat().  This gives
#       the lower triangular of the rmsd matrix *D*.
# 3. Make the rmsd matrix square using make_square()
# 4. Convert the rmsd matrix D into an adjacency matrix *adj_mat*.
#       This requires an rmsd cutoff threshold. Choose rmsd_cutoff=0.5 for now.
# 5. Cluster the adj matrix with greedy(adj_mat)


@jit("f4[:,:](f4[:,:,:])", nopython=True, cache=True)
def make_pairwise_rmsd_mat(X):
    """This function takes a numpy array of coordinates.
    Each row of the array is a 1d vector containing the
    coordoinates of the atoms to be superimposed and rmsd
    calculated.  The function returns the lower triangular
    the RMSD matrix."""

    M = X.shape[0]
    N = X.shape[1]
    O = X.shape[2]
    D = np.zeros((M, M), dtype=np.float32)
    m_com = np.zeros(O, dtype=np.float32)
    t_com = np.zeros(O, dtype=np.float32)
    m = np.zeros((N, O), dtype=np.float32)
    mtrans = np.zeros((O, N), dtype=np.float32)
    mtr = np.zeros((N, O), dtype=np.float32)
    t = np.zeros((N, O), dtype=np.float32)
    c = np.zeros((O, O), dtype=np.float32)
    U = np.zeros((O, O), dtype=np.float32)
    S = np.zeros(O, dtype=np.float32)
    Wt = np.zeros((O, O), dtype=np.float32)
    R = np.zeros((O, O), dtype=np.float32)
    mtr_re = np.zeros(N * O, dtype=np.float32)
    t_re = np.zeros(N * O, dtype=np.float32)
    sub = np.zeros(N * O, dtype=np.float32)
    for i in range(M):
        for j in range(i + 1, M):
            for k in range(O):
                m_com[k] = np.mean(X[i, :, k])
                t_com[k] = np.mean(X[j, :, k])
            m = np.subtract(X[i, :, :], m_com)
            for a in range(N):
                for b in range(O):
                    mtrans[b, a] = m[a, b]
            t = np.subtract(X[j, :, :], t_com)
            c = np.dot(mtrans, t)
            U, S, Wt = np.linalg.svd(c)
            R = np.dot(U, Wt)
            if np.linalg.det(R) < 0.0:
                Wt[-1, :] *= -1.0
                R = np.dot(U, Wt)
            mtr = np.add(np.dot(m, R), t_com)
            q = 0
            for a in range(N):
                for b in range(O):
                    mtr_re[q] = mtr[a, b]
                    t_re[q] = X[j, :, :][a, b]
                    q += 1
            sub = np.subtract(mtr_re, t_re)
            D[i, j] = np.sqrt(1.0 / N * np.dot(sub, sub))
    return D


def make_square(D):
    return D.T + D


def make_adj_mat(D, rmsd_cutoff):
    """Makes an adjacency matrix from the RMSD matrix"""
    adj_mat = np.zeros(D.shape)
    adj_mat[D <= rmsd_cutoff] = 1
    return csr_matrix(adj_mat)


def greedy(adj_mat):
    """Takes an adjacency matrix as input.
        All values of adj_mat are 1 or 0:  1 if <= to cutoff, 0 if > cutoff.
        Can generate adj_mat from data in column format with:
        sklearn.neighbors.NearestNeighbors(metric='euclidean',
        radius=cutoff).fit(data).radius_neighbors_graph(data)

        Returns list of cluster members *all_mems* and centroids *cents*
        (indices of the original coordinate array X) in order of
        largest to smallest.
        """

    if not isinstance(adj_mat, csr_matrix):
        try:
            adj_mat = csr_matrix(adj_mat)
        except:
            print('adj_mat distance matrix must be scipy csr_matrix '
                  '(or able to convert to one)')
            return

    assert adj_mat.shape[0] == adj_mat.shape[1], 'Distance matrix is not square.'

    all_mems = []
    cents = []
    indices = np.arange(adj_mat.shape[0])

    while adj_mat.shape[0] > 0:
        cent = adj_mat.sum(axis=1).argmax()
        cents.append(indices[cent])
        row = adj_mat.getrow(cent)
        tf = ~row.toarray().astype(bool)[0]
        mems = indices[~tf]
        all_mems.append(mems)
        indices = indices[tf]
        adj_mat = adj_mat[tf][:, tf]

    return all_mems, cents

