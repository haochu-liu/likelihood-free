import numpy as np
from localtree import LocalTree


def homoplasy_index(arg, node_site):
    """
    Calculate the homoplasy index for an ARG with mutations.

    The homoplasy index measures the amount of homoplasy (parallel or convergent
    evolution) in the data. It is calculated as 1 - m/s, where m is the minimum
    number of changes required, and s is the actual number of changes observed.

    Parameters
    ----------
    arg : ARG
        An ARG object (from ClonalOrigin_pair).
    node_site : np.ndarray
        Node-site incidence matrix (boolean), typically from add_mutation().
        Shape: (n_nodes, n_sites)

    Returns
    -------
    float
        The homoplasy index, ranging from 0 (no homoplasy) to 1 (maximum homoplasy).
    """
    n_leaf = arg.n
    n_site = arg.node_mat.shape[1]
    m_vec = np.zeros(n_site, dtype=int)
    s_vec = np.zeros(n_site, dtype=int)

    for site_loc in range(n_site):
        # Get local tree at this site
        arg_site = LocalTree(arg, site_loc)
        node_vec = arg_site.node_index
        node_site_vec = node_site[node_vec - 1, site_loc]  # Convert to 0-indexed

        # Compute minimum possible changes
        leaf_states = node_site_vec[:n_leaf]
        if np.any(leaf_states) and not np.all(leaf_states):
            m_vec[site_loc] = 1

        # Compute actual changes using Fitch algorithm
        s_site = 0
        site_dict = {}

        # Initialize leaf nodes
        for i in range(n_leaf):
            node_id = node_vec[i]
            site_dict[node_id] = {int(node_site_vec[i])}

        # Process internal nodes
        for i in range(n_leaf, len(node_vec)):
            parent_node = node_vec[i]
            node_indices = np.where(arg_site.edge[:, 0] == parent_node)[0]
            if len(node_indices) == 2:
                # Coalescent structure (two children)
                children_node = arg_site.edge[node_indices, 1].astype(int)
                child_1_states = site_dict[children_node[0]]
                child_2_states = site_dict[children_node[1]]

                intersec = child_1_states & child_2_states
                if len(intersec) == 0:
                    # Intersection is empty -> mutation occurred
                    site_dict[parent_node] = child_1_states | child_2_states
                    s_site += 1
                else:
                    # Intersection is not empty -> no mutation
                    site_dict[parent_node] = intersec

            elif len(node_indices) == 1:
                # Recombination structure (single child)
                children_node = arg_site.edge[node_indices, 1].astype(int)
                site_dict[parent_node] = site_dict[children_node]

        s_vec[site_loc] = s_site

    # Calculate homoplasy index
    m_sum = np.sum(m_vec)
    s_sum = np.sum(s_vec)

    if m_sum == 0 and s_sum == 0:
        return 0.0

    return 1 - m_sum / s_sum
