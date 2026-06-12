import numpy as np
from Bio import Phylo


def newick_to_tree(tree):
    '''Convert a Newick tree to a ClonalTree object.'''
    counter = 1
    for leaf in tree.get_terminals():
        if leaf.name and leaf.name.isdigit():
            leaf.node_index = int(leaf.name)
        else:
            leaf.node_index = counter
        counter += 1

    terminal_indices = [leaf.node_index for leaf in tree.get_terminals()]
    next_index = max(terminal_indices) + 1 if terminal_indices else 1

    for clade in tree.find_clades(terminal=False, order='postorder'):
        clade.node_index = next_index
        next_index += 1

    edge = []
    node_height = np.full(next_index, np.nan)
    node_height[:(counter-1)] = 0

    for clade in tree.find_clades(order='preorder'):
        for child in clade.clades:
            edge.append([clade.node_index, child.node_index, child.branch_length])

    sort_idx = np.argsort(node_height)
    index_mapping = np.empty_like(sort_idx)
    index_mapping[sort_idx] = np.arange(len(node_height))

    edge = np.array(edge)
    edge[:, :2] = index_mapping[np.int64(edge[:, :2])]
    index_difference = 1 - np.min(edge[:, :2])
    edge[:, :2] += index_difference
    edge = edge[edge[:, 0].argsort()]

    height_to_root = np.zeros(counter - 1)
    for leaf in range(1, counter):
        height_to_leaf = 0
        start_node = leaf
        while start_node < (counter - 1)*2 - 1:
            mat_indices = np.where(edge[:, 1] == start_node)[0]
            parent_node = edge[mat_indices, 0][0]
            height_to_leaf += edge[mat_indices, 2][0]
            start_node = parent_node
        height_to_root[leaf - 1] = height_to_leaf
    
    min_height = np.min(height_to_root)
    diff_height = height_to_root - min_height
    for leaf in range(1, counter):
        mat_indices = np.where(edge[:, 1] == leaf)[0]
        edge[mat_indices, 2] -= diff_height[leaf - 1]
    
    for node in range(counter, (counter - 1)*2):
        mat_indices = np.where(edge[:, 0] == node)[0][0]
        leaf = edge[mat_indices, 1]
        node_height[node-1] = node_height[int(leaf)-1] + edge[mat_indices, 2]

    return edge, node_height[sort_idx]
