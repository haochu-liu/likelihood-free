import numpy as np
from tree import tree


class clonal_tree(tree):
    """Subclass representing a clonal genealogy tree."""
    
    def __init__(self, n):
        """
        Initialize a clonal tree by simulating coalescent backwards in time.
        
        Parameters
        ----------
        n : int
            The number of leaf lineages.
        
        Attributes
        ----------
        edge : np.ndarray
            Array of shape (2*(n-1), 3) containing [node1, node2, length] for each edge.
        node_height : np.ndarray
            Array of length (2*n-1) containing node heights to recent time.
        """
        if not isinstance(n, int) or n < 2:
            raise ValueError("`n` must be a single integer >= 2!")
        
        super().__init__(n)
        self._simulate_clonal_tree()
    
    def _simulate_clonal_tree(self):
        """Simulate coalescent backwards in time to construct the clonal tree."""
        n = self.n
        k = n
        t_sum = 0.0
        
        # Initialize variables for clonal tree
        # edge columns: node1, node2, length
        edge = np.full((2 * (n - 1), 3), np.nan)
        node_height = np.full(2 * n - 1, np.nan)
        node_height[:n] = 0.0  # initialize first n nodes
        
        # Initialize variables
        edge_index = 0
        node_index = n  # 0-indexed, so internal nodes start at n
        pool = list(range(n))  # 0-indexed leaf nodes
        
        # Clonal tree by coalescent only
        while k > 1:
            # Sample a new event time
            rate = k * (k - 1) / 2
            event_time = np.random.exponential(1.0 / rate)
            t_sum += event_time
            
            # Coalescent event: sample two lineages to coalesce
            leaf_indices = np.random.choice(len(pool), size=2, replace=False)
            leaf_node = [pool[leaf_indices[0]], pool[leaf_indices[1]]]
            
            # Append edges
            edge[edge_index, 0] = node_index
            edge[edge_index, 1] = leaf_node[0]
            edge[edge_index, 2] = t_sum - node_height[leaf_node[0]]
            
            edge[edge_index + 1, 0] = node_index
            edge[edge_index + 1, 1] = leaf_node[1]
            edge[edge_index + 1, 2] = t_sum - node_height[leaf_node[1]]
            
            # Append root node height
            node_height[node_index] = t_sum
            
            # Updates for iteration
            pool = [p for p in pool if p not in leaf_node] + [node_index]
            edge_index += 2
            node_index += 1
            k -= 1
        
        self.edge = edge
        self.node_height = node_height
        self.height = np.max(node_height)
        self.length = np.sum(edge[:, 2])
    
    def __repr__(self):
        return f"{self.__class__.__name__}(n={self.n})"
