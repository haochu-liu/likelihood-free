import numpy as np
import sys
sys.path.append('pysimARG')
from clonal_genealogy import ClonalTree
from ClonalOrigin_pair import ARG


@profile
def add_mutation(arg, theta_site):
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
    theta = theta_site * arg.node_mat.shape[1]
    n_mutations = np.random.poisson(theta * arg.length / 2)
    
    # Initialize node_site matrix (boolean)
    node_site = np.zeros(
        (arg.node_mat.shape[0], arg.node_mat.shape[1]), dtype=bool
    )

    # If there is no mutation
    if n_mutations == 0:
        return node_site
    
    # If there are mutations
    # Sample edges with probability proportional to edge length
    edge_probs = arg.edge[:, 2] / arg.length
    mutate_edge = np.random.choice(len(arg.edge), n_mutations, replace=True, p=edge_probs)
    # Sample sites uniformly (0-indexed)
    mutate_site = np.random.choice(arg.node_mat.shape[1], n_mutations, replace=True)

    # Ignore mutations not in the edge material
    keep_mutation = []
    for i in range(n_mutations):
        if arg.edge_mat[mutate_edge[i], mutate_site[i]]:
            keep_mutation.append(i)

    mutate_edge = mutate_edge[keep_mutation]
    mutate_site = mutate_site[keep_mutation]

    if len(keep_mutation) == 0:
        return node_site
    
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


if __name__ == "__main__":
    tree = ClonalTree(n=15)
    rho_site = 0.2
    L = 100000
    delta = 300
    k = 2000
    theta_site = 0.2
    ARG = ARG(tree, rho_site, L, delta, k)
    node_site = add_mutation(ARG, theta_site)
