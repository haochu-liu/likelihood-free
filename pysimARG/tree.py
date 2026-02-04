import numpy as np


class tree:
    """Base class for representing a phylogenetic tree structure."""
    
    def __init__(self, n):
        """
        Initialize a tree.
        n: an integer for the number of leaf lineages.
        """
        self.n = n
        self.edge = None
        self.node_height = None
    
    def __repr__(self):
        return f"{self.__class__.__name__}(n={self.n})"
