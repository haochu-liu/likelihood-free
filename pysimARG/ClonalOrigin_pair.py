import numpy as np
from tree import tree
from clonal_genealogy import ClonalTree
from ClonalOrigin_nodes import ClonalOrigin_nodes
from add_mutation import add_mutation as _add_mutation


class ARG(tree):
    """
    An approximated ancestral recombination graph (ARG) using ClonalOrigin algorithm.
    
    Simulate coalescent and recombination events by ClonalOrigin algorithm at two sites.
    """
    
    def __init__(self, tree_obj, rho_site, L, delta, k):
        """
        Initialize an ARG by simulating recombination on a clonal tree.
        
        Parameters
        ----------
        tree_obj : ClonalTree
            The clonal genealogy.
        rho_site : float
            The recombination parameter per site.
        L : int
            The number of sites.
        delta : float
            The mean of recombinant segment length (must be > 0).
        k : int
            The distance between two sites.
        
        Attributes
        ----------
        edge : np.ndarray
            Edge matrix with columns [node1, node2, length].
        edge_mat : np.ndarray
            Edge material matrix.
        node_height : np.ndarray
            Node heights.
        node_mat : np.ndarray
            Node material matrix.
        node_clonal : np.ndarray
            Boolean array indicating clonal nodes.
        rho : float
            Total recombination rate (L * rho_site).
        L : int
            Number of sites.
        delta : float
            Mean recombinant segment length.
        height : float
            Height of the ARG.
        length : float
            Total length of the ARG.
        """
        # Input validation
        if not isinstance(tree_obj, ClonalTree):
            raise TypeError("Object must be of class 'ClonalTree'")
        if not isinstance(L, int):
            raise ValueError("`L` must be a single integer!")
        if not isinstance(k, int):
            raise ValueError("`k` must be a single integer!")
        if delta <= 0:
            raise ValueError("`delta` must be greater than zero!")
        
        n = tree_obj.n
        super().__init__(n)
        
        self.rho = L * rho_site
        self.L = L
        self.delta = delta
        
        self._simulate(tree_obj, rho_site, L, delta, k)
    
    def _simulate(self, tree_obj, rho_site, L, delta, k):
        """Simulate the ARG using ClonalOrigin algorithm."""
        n = self.n
        clonal_edge = tree_obj.edge.copy()
        clonal_node_height = tree_obj.node_height.copy()
        
        # Tree length
        tree_length = np.sum(clonal_edge[:, 2])
        
        # Initialize recombination edges
        nrow_max = 1000
        # Columns: b_edge, b_height, a_edge, a_height, x, y
        recomb_edge = np.full((nrow_max, 6), np.nan)
        
        # Add recombination sequentially
        n_recomb = 0
        remain_index = np.array([], dtype=int)
        
        for i in range(1, 3):  # i = 1, 2
            if i == 1:
                R_new = np.random.poisson(rho_site * delta * tree_length / 2)
                R_old = 0
            else:  # i == 2
                survive_index = np.where(recomb_edge[:n_recomb, 5] == 1)[0] if n_recomb > 0 else np.array([], dtype=int)
                delta2 = np.sum((1 - 1/delta) ** np.arange(k))
                R_new = np.random.poisson(rho_site * delta2 * tree_length / 2)

                R_old = np.random.binomial(len(survive_index), (1 - 1/delta) ** k)
                remain_index = np.random.choice(survive_index, R_old, replace=False)
            
            if R_new > 0:
                # Expand matrix if needed
                if n_recomb + R_new >= nrow_max:
                    recomb_edge = np.vstack([recomb_edge, np.full((nrow_max, 6), np.nan)])
                    nrow_max = 2 * nrow_max
                
                # Set x and y columns
                recomb_edge[n_recomb:n_recomb + R_new, 4] = i  # x
                recomb_edge[n_recomb:n_recomb + R_new, 5] = i  # y
                
                a_rexp = np.random.exponential(1.0, size=R_new)
                
                # Simulate b_edge (similar to mutation)
                # Sample edges with probability proportional to edge length
                edge_probs = clonal_edge[:, 2] / np.sum(clonal_edge[:, 2])
                recomb_edge[n_recomb:n_recomb + R_new, 0] = np.random.choice(
                    range(1, (2*n-1)), R_new, replace=True, p=edge_probs
                )
                
                for j in range(R_new):
                    idx = n_recomb + j
                    b_edge_idx = int(recomb_edge[idx, 0]) - 1
                    
                    # Simulate b_height
                    recomb_edge[idx, 1] = (
                        np.random.uniform(0, clonal_edge[b_edge_idx, 2]) +
                        clonal_node_height[int(clonal_edge[b_edge_idx, 1])-1]
                    )
                    
                    # Identify a_height
                    # t_above_b: heights of internal nodes minus b_height
                    t_above_b = clonal_node_height[n:2*n-1] - recomb_edge[idx, 1]
                    
                    # Get positive values (nodes above b)
                    positive_mask = t_above_b >= 0
                    positive_t = t_above_b[positive_mask]
                    
                    # i_above_b with 0 prepended
                    i_above_b_full = np.concatenate([[0], np.sort(positive_t)])
                    i_above_b = np.diff(i_above_b_full)
                    
                    # Calculate cumulative values
                    num_intervals = len(i_above_b)
                    lineage_counts = np.arange(num_intervals + 1, 1, -1)
                    cuml_above_b = np.cumsum(i_above_b * lineage_counts)
                    
                    # Determine number of lineages at coalescence time
                    num_lineage = (num_intervals + 1) - np.sum(a_rexp[j] > cuml_above_b)
                    
                    if num_lineage == (num_intervals + 1):
                        recomb_edge[idx, 3] = a_rexp[j] / num_lineage + recomb_edge[idx, 1]
                    else:
                        idx_cuml = num_intervals - num_lineage
                        recomb_edge[idx, 3] = (
                            (a_rexp[j] - cuml_above_b[idx_cuml]) / num_lineage +
                            np.sum(i_above_b[:idx_cuml + 1]) +
                            recomb_edge[idx, 1]
                        )
                    
                    # Simulate a_edge
                    if num_lineage > 1:
                        a_height = recomb_edge[idx, 3]
                        # Find edges that span the a_height
                        pool_edge = np.where(
                            (clonal_node_height[clonal_edge[:, 0].astype(int)-1] >= a_height) &
                            (clonal_node_height[clonal_edge[:, 1].astype(int)-1] <  a_height)
                        )[0] + 1
                        recomb_edge[idx, 2] = np.random.choice(pool_edge)
                    else:
                        # Root edge (using edge index 2*n-2 for 0-indexed)
                        recomb_edge[idx, 2] = 2 * n - 1
            
            if R_old > 0 and len(remain_index) > 0:
                recomb_edge[remain_index, 5] = i
            
            n_recomb = n_recomb + R_new
        
        # Handle case with no recombination
        if n_recomb == 0:
            self.edge = clonal_edge
            self.edge_mat = np.full((2 * (n - 1), 2), True)
            self.node_height = clonal_node_height
            self.node_mat = np.full((2 * n - 1, 2), True)
            self.node_clonal = np.full(2 * n - 1, True)
            self.height = np.max(clonal_node_height)
            self.length = np.sum(clonal_edge[:, 2])
            return
        
        # Trim recomb_edge to actual size
        recomb_edge = recomb_edge[:n_recomb, :]
        
        # Build full ARG
        node_max = 2 * n - 1 + 3 * n_recomb
        edge_max = 2 * (n - 1) + 4 * n_recomb

        edge_matrix = np.full((edge_max, 3), np.nan)
        edge_mat_index = np.full(edge_max, np.nan)
        node_mat = np.full((node_max, 2), np.nan)
        node_info = np.full((node_max, 4), np.nan)
        # Columns: index, height, recomb, clonal

        node_mat[:n, :] = True
        node_info[:, 0] = np.arange(1, node_max + 1)

        # Set clonal node info
        node_info[:2*n-1, 1] = clonal_node_height
        node_info[:2*n-1, 2] = 0
        node_info[:2*n-1, 3] = 1  # True -> 1

        # Set recombination out nodes (b nodes)
        # Interleave: for each recomb event, add two nodes at b_height
        for r in range(n_recomb):
            base_idx = 2 * n - 1 + 2 * r
            node_info[base_idx, 1] = recomb_edge[r, 1]      # b_height
            node_info[base_idx, 2] = -(r + 1)               # negative recomb index (1-indexed)
            node_info[base_idx, 3] = 1                      # clonal = True
            
            node_info[base_idx + 1, 1] = recomb_edge[r, 1]  # same b_height
            node_info[base_idx + 1, 2] = np.nan             # NA
            node_info[base_idx + 1, 3] = 0                  # clonal = False

        # Set recombination in nodes (a nodes)
        for r in range(n_recomb):
            idx = 2 * n - 1 + 2 * n_recomb + r
            node_info[idx, 1] = recomb_edge[r, 3]           # a_height
            node_info[idx, 2] = r + 1                       # positive recomb index (1-indexed)
            node_info[idx, 3] = 1                           # clonal = True
        
        # Sort by height
        sort_order = np.lexsort((node_info[:, 0], node_info[:, 1]))
        node_info = node_info[sort_order, :]

        # Recombination nodes on every edge
        # Use 1-indexed edge indices to match R behavior
        recomb_node = []
        for edge_idx in range(1, 2 * n):
            # Convert to 1-indexed for ClonalOrigin_nodes
            nodes = ClonalOrigin_nodes(recomb_edge, edge_idx)
            recomb_node.append(nodes)
        
        # Build ARG edges and track ancestral material
        i = n  # Start after leaf nodes (0-indexed)
        edge_index = 0

        while i < node_max:
            recomb_val = node_info[i, 2]
            
            if recomb_val == 0:
                # Clonal tree node
                node_index = int(node_info[i, 0])
                leaf_edge = np.where(clonal_edge[:, 0] == node_index)[0]
                leaf_index = np.full(2, np.nan)
                leaf_node = np.full(2, np.nan)

                for le_idx in range(2):
                    le = leaf_edge[le_idx]
                    if len(recomb_node[le]) > 0:
                        # Target node is last element
                        tar_node = recomb_node[le][-1]
                        leaf_index[le_idx] = np.where(node_info[:, 2] == tar_node)[0][0]
                        leaf_node[le_idx] = node_info[int(leaf_index[le_idx]), 0]
                    else:
                        leaf_node[le_idx] = clonal_edge[le, 1]
                        leaf_index[le_idx] = np.where(node_info[:, 0] == leaf_node[le_idx])[0][0]

                # Append edges
                edge_matrix[edge_index:edge_index+2, 0] = i + 1
                edge_matrix[edge_index:edge_index+2, 1] = leaf_index + 1
                edge_matrix[edge_index:edge_index+2, 2] = node_info[i, 1] - node_info[leaf_index.astype(int), 1]
                edge_mat_index[edge_index:edge_index+2] = leaf_index + 1

                # Append root node material
                li0, li1 = int(leaf_index[0]), int(leaf_index[1])
                node_mat[i, :] = np.logical_or(
                    np.nan_to_num(node_mat[li0, :], nan=0).astype(bool),
                    np.nan_to_num(node_mat[li1, :], nan=0).astype(bool)
                )

                edge_index += 2
                i += 1
                
            elif recomb_val < 0:
                # Recombination edge out node
                node_index = node_info[i:i+2, 0].astype(int)
                recomb_idx = int(abs(node_info[i, 2])) - 1
                leaf_edge = int(recomb_edge[recomb_idx, 0]) - 1

                tar_node = np.where(recomb_node[leaf_edge] == node_info[i, 2])[0]
                if tar_node == 0:
                    leaf_node = clonal_edge[leaf_edge, 1]
                else:
                    leaf_node = node_info[np.where(recomb_node[leaf_edge][tar_node-1] == node_info[:, 2])[0][0], 0]

                leaf_index = int(np.where(node_info[:, 0] == leaf_node)[0][0])

                # Append edges
                edge_matrix[edge_index:edge_index+2, 0] = [i + 1, i + 2]
                edge_matrix[edge_index:edge_index+2, 1] = leaf_index + 1
                edge_matrix[edge_index:edge_index+2, 2] = node_info[i, 1] - node_info[leaf_index, 1]
                edge_mat_index[edge_index:edge_index+2] = [i + 1, i + 2]

                x = int(recomb_edge[recomb_idx, 4])
                y = int(recomb_edge[recomb_idx, 5])

                # Append root node material
                node_mat[i:i+2, :] = False
                # x:y in R is inclusive, in Python x-1:y is equivalent for 1-indexed to 0-indexed
                node_mat[i + 1, x-1:y] = np.nan_to_num(node_mat[leaf_index, x-1:y], nan=0).astype(bool)

                # For the complement (-(x:y) in R means all except x:y)
                mask = np.ones(2, dtype=bool)
                mask[x-1:y] = False
                node_mat[i, mask] = np.nan_to_num(node_mat[leaf_index, mask], nan=0).astype(bool)

                edge_index += 2
                i += 2
                
            elif recomb_val > 0:
                # Recombination edge in node
                node_index = int(node_info[i, 0])
                recomb_idx = int(node_info[i, 2]) - 1
                leaf_edge = int(recomb_edge[recomb_idx, 2]) - 1

                tar_node = np.where(recomb_node[leaf_edge] == node_info[i, 2])[0]
                if tar_node == 0:
                    if leaf_edge== 2*n - 2:
                        leaf_node = 2*n - 1
                    else:
                        leaf_node = clonal_edge[leaf_edge, 1]
                else:
                    leaf_node = node_info[np.where(recomb_node[leaf_edge][tar_node-1] == node_info[:, 2])[0][0], 0]

                leaf_index = np.full(2, np.nan)
                leaf_index[0] = int(np.where(node_info[:, 0] == leaf_node)[0][0])
                leaf_index[1] = int(np.where(node_info[:, 2] == (-node_info[i, 2]))[0][0]) + 1

                # Append edges
                edge_matrix[edge_index:edge_index+2, 0] = i + 1
                edge_matrix[edge_index:edge_index+2, 1] = leaf_index + 1
                edge_matrix[edge_index:edge_index+2, 2] = node_info[i, 1] - node_info[leaf_index.astype(int), 1]
                edge_mat_index[edge_index:edge_index+2] = leaf_index + 1

                # Append root node material
                li0, li1 = int(leaf_index[0]), int(leaf_index[1])
                node_mat[i, :] = np.logical_or(
                    np.nan_to_num(node_mat[li0, :], nan=0).astype(bool),
                    np.nan_to_num(node_mat[li1, :], nan=0).astype(bool)
                )

                edge_index += 2
                i += 1
            else:
                # NaN case - skip
                i += 1
        
        # Store results
        self.edge = edge_matrix
        self.edge_mat = node_mat[edge_mat_index.astype(int) - 1, :]
        self.node_height = node_info[:, 1]
        self.node_mat = node_mat
        self.node_clonal = node_info[:, 3].astype(bool)
        self.height = np.max(self.node_height)
        self.length = np.sum(self.edge[:, 2])
    
    def __repr__(self):
        return f"ARG(n={self.n}, rho={self.rho}, L={self.L}, delta={self.delta})"
    
    def add_mutation(self, theta_site):
        return _add_mutation(self, theta_site)
    
    def to_dict(self):
        """Returns the ARG as a dictionary, similar to the R list output."""
        return {
            "edge": self.edge,
            "edge_mat": self.edge_mat,
            "node_height": self.node_height,
            "node_mat": self.node_mat,
            "node_clonal": self.node_clonal,
            "n": self.n,
            "rho": self.rho,
            "L": self.L,
            "delta": self.delta,
            "height": self.height,
            "length": self.length
        }
