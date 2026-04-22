import numpy as np
import scipy.stats as stats
from ClonalOrigin_ARG import ARG
from add_mutation_truncated import add_mutation_truncated
from G3_test import G3_test
from LD_r import LD_r
from homoplasy_index import homoplasy_index


def ClonalOrigin_seq_sim(tree, rho_site, theta_site, L, delta, k_vec=[50, 200, 2000]):
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
    k_vec : list of int in ascending order, optional
        An ascending list of three distance values between two sites. Default is [50, 200, 2000].
    
    Returns
    -------
    np.ndarray
        A 8-dimensional vector as the summary statistics of simulations.
    """
    if max(k_vec) >= L:
        raise ValueError("Site distance cannot be greater or equal to the number of sites!")
    if not np.all(np.sort(k_vec) == k_vec):
        raise ValueError("`k_vec` must be in ascending order!")
    
    s_size = int(2 * len(k_vec) + 3)
    s_vec = np.full(s_size, np.nan)
    tree_width = tree.n

    ARG_sim = ARG(tree, rho_site, L, delta, L, "seq")
    node_site = add_mutation_truncated(ARG_sim, theta_site)

    # Summary statistics LD r and G3 test
    for i in range(len(k_vec)):
        mat_num = node_site.shape[1] - k_vec[i]
        v_r = np.full(mat_num, np.nan)
        v_g3 = np.full(mat_num, np.nan)
        for j in range(mat_num):
            mat = node_site[:tree_width, [j, j+k_vec[0]]]
            v_r[j] = LD_r(mat)
            v_g3[j] = G3_test(mat)
        s_vec[i] = np.mean(v_r)
        s_vec[i + len(k_vec)] = np.mean(v_g3)
    
    # Summary statistic homoplasy index
    s_vec[2 * len(k_vec)] = homoplasy_index(ARG_sim, node_site)

    # Summary statistic proportion of segregating sites
    log_weight = 0
    for i in range(ARG_sim.edge_mat.shape[1]):
        arg_i_length = np.sum(ARG_sim.edge[ARG_sim.edge_mat[:, i] == 1, 2])
        log_weight += stats.expon.logcdf(arg_i_length, scale=1/(theta_site/2))
    s_vec[2 * len(k_vec) + 1] = np.exp(log_weight) / L

    s_vec[2 * len(k_vec) + 2] = L
    
    return s_vec
