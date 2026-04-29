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
        
        leaf_distances = []
        descendant_leaves = clade.get_terminals()
        for leaf in descendant_leaves:
            dist = tree.distance(clade, leaf)
            leaf_distances.append(dist)
        node_height[clade.node_index] = np.median(leaf_distances)
    
    sort_idx = np.argsort(node_height)
    index_mapping = np.empty_like(sort_idx)
    index_mapping[sort_idx] = np.arange(len(node_height))

    edge = np.array(edge)
    edge[:, :2] = index_mapping[np.int64(edge[:, :2])]
    edge[:, :2] += 1

    return edge, node_height[sort_idx]
