import numpy as np


def clade_homoplasy(tree, node_site):
    """
    Calculate the clade homoplasy index for an ARG with mutations.

    Given a clonal genealogy, it is divided into its two largest clades and for biallelic sites
    if both alleles are present in both clades, we say that the site is clade homoplasic.

    Parameters
    ----------
    tree : ClonalTree
        A clonal tree object representing the clonal genealogy for the ARG.
    node_site : np.ndarray
        Node-site incidence matrix (boolean), typically from add_mutation().
        Shape: (n_nodes, n_sites)

    Returns
    -------
    float
        The clade homoplasy index, ranging from 0 (no homoplasy) to 1 (maximum homoplasy).
    """
    n_leaf = tree.n
    n_site = node_site.shape[1]
    hi_vec = np.full(n_site, np.nan)

    clade1 = []
    clade2 = []

    for site_loc in range(n_site):
        leaf_states = node_site[:n_leaf, site_loc]
        partition1 = leaf_states[clade1]
        partition2 = leaf_states[clade2]

        len1 = len(np.unique(partition1))
        len2 = len(np.unique(partition2))

        if len1 > 1 and len2 > 1:
            hi_vec[site_loc] = 1.0
        else:
            hi_vec[site_loc] = 0.0

    return np.mean(hi_vec)
