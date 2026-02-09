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
        - mutation: np.ndarray with columns [edge_index, pos, site]
        - node_site: np.ndarray (if binary=True) incidence matrix
        - node_gene: list of lists containing mutation sites at each node
    """
    theta = theta_site * arg.L
    n_mutations = np.random.poisson(theta * arg.length / 2)
    
    # Initialize mutation matrix
    arg.mutation = np.full((n_mutations, 3), np.nan)
    # Columns: edge_index, pos, site
    
    num_nodes = arg.node_mat.shape[0]
    
    # Initialize node gene information
    arg.node_gene = [[] for _ in range(num_nodes)]
    
    # If there is no mutation
    if n_mutations == 0:
        return arg
    
    # If there are mutations
    # Sample edges with probability proportional to edge length
    edge_probs = arg.edge[:, 2] / arg.length
    mutate_edge = np.random.choice(
        len(arg.edge), n_mutations, replace=True, p=edge_probs
    )
    # Sample sites uniformly (0-indexed)
    mutate_site = np.random.choice(arg.L, n_mutations, replace=True)
    
    # Ignore mutations not in the edge material
    keep_mutation = []
    for i in range(n_mutations):
        if arg.edge_mat[mutate_edge[i], mutate_site[i]]:
            keep_mutation.append(i)
    
    if len(keep_mutation) == 0:
        arg.mutation = np.full((0, 3), np.nan)
        return arg
    
    mutate_edge = mutate_edge[keep_mutation]
    mutate_site = mutate_site[keep_mutation]
    arg.mutation = arg.mutation[keep_mutation, :]
    
    # Store mutation information
    arg.mutation[:, 0] = mutate_edge
    arg.mutation[:, 2] = mutate_site
    for i in range(len(mutate_edge)):
        arg.mutation[i, 1] = np.random.uniform(0, arg.edge[mutate_edge[i], 2])
    
    # Simulate the mutations at every node
    # Process edges from last to first (bottom-up in the tree)
    for i in range(len(arg.edge) - 1, -1, -1):
        # Get mutations on this edge
        edge_mutations = arg.mutation[arg.mutation[:, 0] == i, 2].astype(int).tolist()
        
        parent_idx = int(arg.edge[i, 0])
        child_idx = int(arg.edge[i, 1])
        
        # Get parent sequence and add edge mutations
        parent_seq = arg.node_gene[parent_idx].copy()
        parent_seq.extend(edge_mutations)
        
        # Remove mutations for sites not in edge material
        for j in range(arg.L):
            if arg.edge_mat[i, j] == 0:
                parent_seq = [m for m in parent_seq if m != j]
        
        # Combine with existing child sequence and sort
        arg.node_gene[child_idx] = sorted(arg.node_gene[child_idx] + parent_seq)
    
    return arg
