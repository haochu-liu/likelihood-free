import numpy as np
from localtree import LocalTree
from localtree_simbac import get_tree_at_position
from newick_to_tree import newick_to_tree


def homoplasy_index(blocks, node_site, start_pos, L):
    """
    Calculate the homoplasy index for an ARG with mutations.

    The homoplasy index measures the amount of homoplasy (parallel or convergent
    evolution) in the data. It is calculated as 1 - m/s, where m is the minimum
    number of changes required, and s is the actual number of changes observed.

    Parameters
    ----------
    blocks : list
        List of (span, tree) tuples representing the blocks in the ARG.
    node_site : np.ndarray
        Node-site incidence matrix (boolean), typically from add_mutation().
        Shape: (n_nodes, n_sites)
    start_pos : int
        Starting position of the gene region.
    L : int
        Length of the gene region.


    Returns
    -------
    float
        The homoplasy index, ranging from 0 (no homoplasy) to 1 (maximum homoplasy).
    """
    n_leaf = node_site.shape[0]
    m_vec = np.zeros(L, dtype=int)
    s_vec = np.zeros(L, dtype=int)

    # Initialize local tree
    local_tree, start, end = get_tree_at_position(blocks, start_pos)
    local_edge, local_node_height = newick_to_tree(local_tree)

    local_pos = start_pos

    for i in range(L):
        # Get local tree at this site
        if not start <= local_pos <= end:
            local_tree, start, end = get_tree_at_position(blocks, local_pos)
            local_edge, local_node_height = newick_to_tree(local_tree)

        node_vec = np.arange(1, 2*n_leaf)

        # Compute minimum possible changes
        leaf_states = node_site[:, local_pos - 1]  # Convert to 0-indexed
        if np.any(leaf_states) and not np.all(leaf_states):
            m_vec[i] = 1

        # Compute actual changes using Fitch algorithm
        s_site = 0
        site_dict = {}

        # Initialize leaf nodes
        site_dict = {node: {int(state)} for node, state in zip(node_vec[:n_leaf], leaf_states)}

        # Process internal nodes
        for i in range(n_leaf, len(node_vec)):
            parent_node = node_vec[i]
            node_indices = np.where(local_edge[:, 0] == parent_node)[0]

            # Coalescent structure (two children)
            children_node = local_edge[node_indices, 1].astype(int)
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

        s_vec[i] = s_site
        local_pos += 1

    # Calculate homoplasy index
    m_sum = np.sum(m_vec)
    s_sum = np.sum(s_vec)

    if m_sum == 0 and s_sum == 0:
        return 0.0

    return 1 - m_sum / s_sum
