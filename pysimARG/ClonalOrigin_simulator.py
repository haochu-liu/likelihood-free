import numpy as np
from ClonalOrigin_pair import ARG
from G3_test import G3_test
from LD_r import LD_r


def ClonalOrigin_simulator(tree, rho_site, theta_site, L, delta, N,
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
    N : int
        Number of pairs to simulate.
    k_vec : list of int, optional
        A list of three distance values between two sites. Default is [50, 200, 2000].
    
    Returns
    -------
    np.ndarray
        A 7-dimensional vector as the summary statistics of simulations.
    """
    if max(k_vec) > L:
        raise ValueError("Site distance cannot be greater than the number of sites!")
    if len(k_vec) != 3:
        raise ValueError("`k_vec` must be a list of three integer values!")
    
    s_vec = np.full(7, np.nan)
    tree_width = tree.n
    v_s = np.full(N * 3, np.nan)
    
    for j in range(3):
        v_r = np.full(N, np.nan)
        v_g3 = np.full(N, np.nan)
        
        for i in range(N):
            arg = ARG(tree, rho_site, L, delta, k_vec[j])
            arg.add_mutation(theta_site)
            mat = arg.node_site[:tree_width, :]
            
            v_r[i] = LD_r(mat)
            v_g3[i] = G3_test(mat)
            v_s[i + j * N] = np.any(mat[:, 0].astype(bool))
        
        s_vec[j] = np.mean(v_r)
        s_vec[j + 3] = np.mean(v_g3)
    
    s_vec[6] = np.mean(v_s)
    
    return s_vec
