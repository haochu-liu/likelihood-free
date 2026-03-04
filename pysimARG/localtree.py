import numpy as np
from copy import deepcopy


class LocalTree:
    """
    Local tree at a chosen site, constructed from an ARG object.

    Attributes
    ----------
    edge : np.ndarray
        Filtered edge matrix with columns [node1, node2, length].
    edge_mat : np.ndarray
        Filtered edge material matrix.
    edge_index : np.ndarray
        Indices of the kept edges from the original ARG.
    node_height : np.ndarray
        Node heights (inherited from ARG).
    node_mat : np.ndarray
        Node material matrix (inherited from ARG).
    node_clonal : np.ndarray
        Boolean array indicating clonal nodes (inherited from ARG).
    n : int
        Number of leaf lineages.
    """

    def __init__(self, ARG, location):
        """
        Construct a local tree from an ARG object at a given site.

        Parameters
        ----------
        ARG : ARG
            An ARG object (from ClonalOrigin_pair).
        location : int
            The site location for the local tree (0-indexed column of edge_mat).
        """
        # Copy ARG attributes so we don't mutate the original
        edge = deepcopy(ARG.edge)
        edge_mat = deepcopy(ARG.edge_mat)

        # Keep edges whose material covers the given site
        keep_edge = np.where(edge_mat[:, location].astype(bool))[0]

        edge = edge[keep_edge, :]
        edge_mat = edge_mat[keep_edge, :]
        edge_index = keep_edge

        # Delete long root edges:
        # Sort by parent node (column 0)
        sort_order = np.argsort(edge[:, 0])
        edge = edge[sort_order, :]
        edge_mat = edge_mat[sort_order, :]
        edge_index = edge_index[sort_order]

        # Find duplicated parent nodes (nodes that appear more than once as parent)
        parents = edge[:, 0]
        duplicated_forward = np.zeros(len(parents), dtype=bool)
        duplicated_backward = np.zeros(len(parents), dtype=bool)

        seen = {}
        for i, val in enumerate(parents):
            if val in seen:
                duplicated_forward[i] = True
                duplicated_forward[seen[val]] = True
            seen[val] = i

        # Combine: an element is "duplicated" if it appears more than once
        duplicated_edge = duplicated_forward

        # Find the last index that is duplicated
        dup_indices = np.where(duplicated_edge)[0]
        if len(dup_indices) > 0:
            last_duplicated = dup_indices[-1]
            # Keep edges up to and including the last duplicated entry (0-indexed, so +1)
            edge = edge[:last_duplicated + 1, :]
            edge_mat = edge_mat[:last_duplicated + 1, :]
            edge_index = edge_index[:last_duplicated + 1]

        # Store attributes
        self.edge = edge
        self.edge_mat = edge_mat
        self.edge_index = edge_index
        self.node_height = deepcopy(ARG.node_height)
        self.node_mat = deepcopy(ARG.node_mat)
        self.node_clonal = deepcopy(ARG.node_clonal)
        self.n = ARG.n
        self.location = location

    def __repr__(self):
        return f"LocalTree(n={self.n}, location={self.location})"
