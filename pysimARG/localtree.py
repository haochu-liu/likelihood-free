import numpy as np
from copy import deepcopy


class LocalTree:
    """
    Local tree at a chosen site, constructed from an ARG object.

    Attributes
    ----------
    edge : np.ndarray
        Filtered edge matrix with columns [node1, node2, length].
    edge_index : np.ndarray
        Indices of the kept edges from the original ARG.
    node_index : np.ndarray
        The nodes in the local tree.
    n : int
        Number of leaf lineages.
    location : int
        The site location for the local tree.
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
        edge = deepcopy(ARG.edge)

        # Keep edges whose material covers the given site
        keep_edge = np.where(ARG.edge_mat[:, location])[0]
        edge = edge[keep_edge, :]

        # 1. Sort edge by the first column (index 0)
        sort_order = np.argsort(edge[:, 0])
        edge = edge[sort_order]

        # 2. Find duplicates in the first column
        col_0 = edge[:, 0]
        unique_vals, counts = np.unique(col_0, return_counts=True)
        duplicate_vals = unique_vals[counts > 1]
        duplicated_edge = np.isin(col_0, duplicate_vals)

        # 3. Find the index of the last True value in the boolean array
        last_duplicated = np.where(duplicated_edge)[0][-1]

        # 4. Slice the arrays (Python slicing is exclusive, so we add 1)
        edge = edge[:last_duplicated + 1]
        keep_edge = keep_edge[:last_duplicated + 1]

        # Store attributes
        self.edge = edge
        self.edge_index = keep_edge
        self.node_index = np.unique(edge[:, :2])
        self.n = ARG.n
        self.location = location

    def __repr__(self):
        return f"LocalTree(n={self.n}, location={self.location})"
