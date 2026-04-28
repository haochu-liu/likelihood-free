import numpy as np
import scipy.stats as stats
from ClonalOrigin_ARG import ARG
from add_mutation import add_mutation
from G3_test import G3_test
from LD_r import LD_r
from homoplasy_index import homoplasy_index


def ClonalOrigin_seq_sim(tree, rho_site, theta_site, L, delta):
    """
    Simulate approximated ARG using ClonalOrigin with sequential method.
    
    Simulate pairs of sites in the ClonalOrigin approximation from a given clonal tree.
    Provide summary statistics for pairs in three different distances.
    
    Parameters
    ----------
    tree : ClonalTree
        The clonal genealogy.
    rho_site : float
        The recombination parameter per site.
    theta_site : float
        The mutation rate per site.
    L : int
        The number of sites.
    delta : float
        The mean of recombinant segment length.
    
    Returns
    -------
    np.ndarray
        A 7-dimensional vector as the summary statistics of simulations.
    """
    if not isinstance(L, int):
        raise ValueError("`L` must be a single integer!")

    s_vec = np.full(7, np.nan)
    tree_width = tree.n

    ARG_sim = ARG(tree, rho_site, L, delta, L, "seq")
    node_site = add_mutation(ARG_sim, theta_site)
    mat = node_site[:tree_width, :]

    # Identify segregating sites
    has_true = mat.any(axis=0)
    has_false = ~mat.all(axis=0)
    idx_seg = np.where(has_true & has_false)[0]

    # Summary statistics LD r and G3 test
    ld_near, ld_far, g3_near, g3_far = 0, 0, 0, 0
    if idx_seg.size >= 2:
        for i in range(idx_seg.size - 1):
            for j in range(i + 1, idx_seg.size):
                dist_ij = idx_seg[j] - idx_seg[i]
                idx_pair = [idx_seg[i], idx_seg[j]]
                if dist_ij < L/2:
                    ld_near += LD_r(mat[:, idx_pair])
                    g3_near += G3_test(mat[:, idx_pair])
                else:
                    ld_far += LD_r(mat[:, idx_pair])
                    g3_far += G3_test(mat[:, idx_pair])
        
        s_far = (int(L/2) + 1) * (int(L/2)) / 2
        s_near = L * (L - 1) / 2 - s_far
        s_vec[0] = ld_near / s_near
        s_vec[1] = ld_far / s_far
        s_vec[2] = g3_near / s_near
        s_vec[3] = g3_far / s_far
    else:
        s_vec[:4] = 0
    
    # Summary statistic homoplasy index
    s_vec[4] = homoplasy_index(ARG_sim, node_site)

    # Summary statistic proportion of segregating sites
    count_S = idx_seg.size
    s_vec[5] = count_S / L

    # Add the length of sequence as a summary statistic
    s_vec[6] = L
    
    return s_vec
