import numpy as np
from scipy.stats import poisson

from pysimARG.ClonalOrigin_pair import ARG


def truncated_poisson(lam, size=1):
    """
    Samples from a Zero-Truncated Poisson distribution.
    """
    p0 = np.exp(-lam)

    eps = np.finfo(float).eps
    u = np.random.uniform(low=eps+p0, high=1.0, size=size)

    return int(poisson.ppf(u, lam).astype(int)[0])


def add_mutation_truncated(arg, theta_site):
    """
    Add mutations uniformly onto the edges of ARG with infinite site assumption.
    
    Parameters
    ----------
    arg : ARG
        An ARG object to add mutations to.
    theta_site : float
        The mutation rate per site.
    
    Returns
    -------
    node_site : np.ndarray incidence matrix
    """
    num_site = arg.node_mat.shape[1]
    n_list = np.zeros(num_site, dtype=int)
    mutate_edge = list()
    mutate_site = list()

    # Initialize node_site matrix (boolean)
    node_site = np.zeros(
        (arg.node_mat.shape[0], arg.node_mat.shape[1]), dtype=bool
    )

    # Simulate mutations for each site
    for i in range(num_site):
        # Length of local tree without reduction
        local_edge = np.where(arg.edge_mat[:, i])[0]
        local_length = arg.edge[local_edge, 2]
        # Truncated Poisson distribution
        local_n = truncated_poisson(theta_site * np.sum(local_length) / 2)
        n_list[i] = local_n
        mutate_site.extend([i] * local_n)
        # Simulate edges
        probs = local_length / np.sum(local_length)
        mutate_edge.extend(np.random.choice(local_edge, local_n, replace=True, p=probs).tolist())
    
    # Simulate the mutations at every node
    # Process edges from last to first (bottom-up in the tree)
    for i in range(len(arg.edge) - 1, -1, -1):
        edge_mutation = mutate_site[mutate_edge == i]

        parent_idx = int(arg.edge[i, 0]) - 1
        child_idx = int(arg.edge[i, 1]) - 1

        parent_seq = node_site[parent_idx, :].copy()
        flip_counts = np.bincount(edge_mutation, minlength=len(parent_seq))
        parent_seq = np.logical_xor(parent_seq, flip_counts % 2 == 1)

        material_range = arg.edge_mat[i, :]
        node_site[child_idx, material_range==True] = parent_seq[material_range==True]
    
    return node_site
