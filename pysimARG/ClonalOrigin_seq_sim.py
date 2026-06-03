import numpy as np
from ClonalOrigin_ARG import ARG
from add_mutation import add_mutation
from segment_summary_stats import segment_summary_stats


def ClonalOrigin_seq_sim(tree, rho_site, theta_site, L, delta):
    """
    Simulate approximated ARG using ClonalOrigin with sequential method.
    
    Simulate pairs of sites in the ClonalOrigin approximation from a given clonal tree.
    Provide summary statistics for pairs in three different distances.
    Kelly's Zns and Hudson's Rm computation are integrated.
    
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
        A 46-dimensional vector as the summary statistics of simulations.
    """
    if not isinstance(L, int):
        raise ValueError("`L` must be a single integer!")

    tree_width = tree.n

    # Simulate ARG and add mutations
    ARG_sim = ARG(tree, rho_site, L, delta, L, "seq")
    node_site = add_mutation(ARG_sim, theta_site)
    # mat = node_site[:tree_width, :]
    mat = node_site[:tree_width, :] != node_site[0]

    # Compute summary statistics
    s_vec = segment_summary_stats(tree, mat)
    
    return s_vec
