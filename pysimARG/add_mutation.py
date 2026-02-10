import numpy as np


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
    arg
        The ARG object with added mutation attributes:
        - node_site: np.ndarray (if binary=True) incidence matrix
    """
    theta = theta_site * arg.L
    n_mutations = np.random.poisson(theta * arg.length / 2)
    
    # Initialize node_site matrix (boolean)
    arg.node_site = np.zeros(
        (arg.node_mat.shape[0], arg.node_mat.shape[1]), dtype=bool
    )

    # If there is no mutation
    if n_mutations == 0:
        return arg
    
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
        return arg
    
    # Simulate the mutations at every node
    # Process edges from last to first (bottom-up in the tree)
    for i in range(len(arg.edge) - 1, -1, -1):
        # Get mutations on this edge (sites that have mutations on edge i)
        edge_mutation = mutate_site[mutate_edge == i]

        parent_idx = int(arg.edge[i, 0]) - 1
        child_idx = int(arg.edge[i, 1]) - 1

        # Get parent sequence
        parent_seq = arg.node_site[parent_idx, :].copy()

        # Count mutations at each site (equivalent to R's tabulate)
        flip_counts = np.bincount(edge_mutation, minlength=len(parent_seq))

        # XOR with whether flip count is odd
        parent_seq = np.logical_xor(parent_seq, flip_counts % 2 == 1)

        # Get material range for this edge
        material_range = arg.edge_mat[i, :]

        # Update child node sequence only where there's material
        arg.node_site[child_idx, material_range==True] = parent_seq[material_range==True]
    
    return arg
