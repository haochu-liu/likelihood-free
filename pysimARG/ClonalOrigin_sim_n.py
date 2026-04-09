import numpy as np
import scipy.stats as stats
from ClonalOrigin_pair import ARG
from add_mutation_truncated import add_mutation_truncated
from G3_test import G3_test
from LD_r import LD_r
from homoplasy_index import homoplasy_index


def ClonalOrigin_sim_n(tree, rho_site, theta_site, L, delta, N,
                       k_vec=[50, 200, 2000]):
    """
    Simulate approximated ARG using ClonalOrigin with pair method.
    
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
    N : list of int
        A list of integers representing the number of pairs to simulate at each distance.
    k_vec : list of int, optional
        A list of three distance values between two sites. Default is [50, 200, 2000].
    
    Returns
    -------
    np.ndarray
        A len(N) * 9 size matrix as the summary statistics of simulations.
    """
    if max(k_vec) > L:
        raise ValueError("Site distance cannot be greater than the number of sites!")
    if len(k_vec) != 3:
        raise ValueError("`k_vec` must be a list of three integer values!")
    
    desc_indices = np.argsort(N)[::-1]
    N_max = int(N[desc_indices[0]])
    
    s_mat = np.full((len(N), 9), np.nan)
    s_mat[:, 8] = N
    tree_width = tree.n
    v_h = np.full((3, N_max), np.nan)
    v_s = np.full((3, N_max), np.nan)
    v_r = np.full((3, N_max), np.nan)
    v_g3 = np.full((3, N_max), np.nan)

    for j in range(3):
        arg_list = [None] * N_max
        log_weight_list = [None] * N_max
        for i in range(N_max):
            ARG_i = ARG(tree, rho_site, L, delta, k_vec[j])
            arg_list[i] = ARG_i
            arg_1_length = np.sum(ARG_i.edge[ARG_i.edge_mat[:, 0] == 1, 2])
            arg_2_length = np.sum(ARG_i.edge[ARG_i.edge_mat[:, 1] == 1, 2])
            log_weight_list[i] = stats.expon.logcdf(arg_1_length, scale=1/(theta_site/2)) + stats.expon.logcdf(arg_2_length, scale=1/(theta_site/2))
        
        log_weight_list = np.array(log_weight_list)
        shifted_log_p = log_weight_list - np.max(log_weight_list)
        weights = np.exp(shifted_log_p)
        probs = weights / np.sum(weights)

        sampled_values = np.random.choice(N_max, size=N_max, replace=True, p=probs)
        
        for i in range(N_max):
            node_site = add_mutation_truncated(arg_list[sampled_values[i]], theta_site)
            mat = node_site[:tree_width, :]
            
            v_r[j, i] = LD_r(mat)
            v_g3[j, i] = G3_test(mat)
            v_h[j, i] = homoplasy_index(arg_list[sampled_values[i]], node_site)
            v_s[j, i] = weights[sampled_values[i]]
    
    for i in range(len(N)):
        idx = desc_indices[i]
        N_val = N[idx]
        if i > 0:
            choose_idx = np.random.randint(0, N_max, size=N_val)
        else:
            choose_idx = np.arange(N_val)

        s_mat[idx, 0] = np.mean(v_r[0, choose_idx])
        s_mat[idx, 1] = np.mean(v_r[1, choose_idx])
        s_mat[idx, 2] = np.mean(v_r[2, choose_idx])
        s_mat[idx, 3] = np.mean(v_g3[0, choose_idx])
        s_mat[idx, 4] = np.mean(v_g3[1, choose_idx])
        s_mat[idx, 5] = np.mean(v_g3[2, choose_idx])
        s_mat[idx, 6] = np.mean(v_h[:, choose_idx])
        s_mat[idx, 7] = np.mean(v_s[:, choose_idx])
    
    return s_mat
