import numpy as np
from tree import tree


class ClonalTree(tree):
    """
    Clonal Genealogy simulation.
    Simulate coalescent backwards in time to construct a clonal tree.
    """

    def __init__(self, n: int):
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
        if not isinstance(n, int) or n < 1:
            raise ValueError("`n` must be a single positive integer!")
        
        super().__init__(n)
        self._simulate()

    def _simulate(self):
        n = self.n
        k = n
        t_sum = 0.0

        # Initialize variables for clonal tree
        # Cols: node1 (parent), node2 (child), length
        self.edge = np.zeros((2 * (n - 1), 3))
        self.node_height = np.zeros(2 * n - 1)
        self.node_height[:n] = 0.0

        # Initialize variables and vector
        edge_index = 0
        node_index = n + 1
        pool = list(range(1, n+1))

        # Clonal tree by coalescent only
        while k > 1:
            # Sample a new event time
            # Rate is k*(k-1)/2. Python uses scale = 1/rate
            rate = k * (k - 1) / 2.0
            event_time = np.random.exponential(scale=1.0/rate)
            t_sum += event_time

            # Coalescent event: sample 2 distinct nodes from the pool
            # Using random.choice from numpy without replacement
            leaf_nodes = np.random.choice(pool, size=2, replace=False)

            # Append edges
            # We add two edges: Parent -> Child 1 and Parent -> Child 2
            
            # 1. Parent index (node1)
            self.edge[edge_index : edge_index+2, 0] = node_index
            
            # 2. Child indices (node2)
            self.edge[edge_index : edge_index+2, 1] = leaf_nodes
            
            # 3. Edge lengths (current time - height of child)
            child_heights = self.node_height[leaf_nodes]
            self.edge[edge_index : edge_index+2, 2] = t_sum - child_heights

            # Append root node height
            self.node_height[node_index-1] = t_sum

            # Updates for iteration
            pool = [x for x in pool if x not in leaf_nodes]
            pool.append(node_index)

            edge_index += 2
            node_index += 1
            k -= 1

        self.height = np.max(self.node_height)
        self.length = np.sum(self.edge[:, 2])

    def __repr__(self):
        return f"{self.__class__.__name__}(n={self.n})"

    def to_dict(self):
        """Returns the tree as a dictionary, similar to the R list output."""
        return {
            "edge": self.edge,
            "node_height": self.node_height,
            "height": self.height,
            "length": self.length,
            "n": self.n
        }
