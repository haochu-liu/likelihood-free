class tree:
    """Base class for representing a phylogenetic tree structure."""
    
    def __init__(self, n):
        """
        Initialize a tree.
        n: an integer for the number of leaf lineages.
        """
        self.n = n
    
    def __repr__(self):
        return f"{self.__class__.__name__}(n={self.n})"


class clonal_tree(tree):
    """Subclass representing a clonal genealogy tree."""
    
    def __init__(self, n):
        """
        Initialize a clonal tree.
        n: an integer for the number of leaf lineages.
        """
        super().__init__(n)
    
    def __repr__(self):
        return f"{self.__class__.__name__}(n={self.n})"
    