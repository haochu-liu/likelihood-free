import numpy as np
from tree import tree
from clonal_genealogy import ClonalTree
from pair_simulator import pair_simulator
from seq_simulator import seq_simulator


class ARG(tree):
    """
    An approximated ancestral recombination graph (ARG) using ClonalOrigin algorithm.
    
    Simulate coalescent and recombination events by ClonalOrigin algorithm at two sites.
    """
    
    def __init__(self, tree_obj, rho_site, L, delta, k, type):
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
            The distance between two sites in pair simulator,
            and the sequence length in sequence simulator.
        type: str
            The type of simulator used (pair or seq)
        
        Attributes
        ----------
        type : str
            The type of simulator used (pair or seq)
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
        if not (type == "seq" or type == "pair"):
            raise ValueError("`type` must be a string equal to 'pair' or 'seq'.")
        
        n = tree_obj.n
        super().__init__(n)
        
        self.rho = L * rho_site
        self.L = L
        self.delta = delta
        self.type = type
        
        self._simulate(tree_obj, rho_site, L, delta, k, type)
    
    def _simulate(self, tree_obj, rho_site, L, delta, k, type):
        if type == "pair":
            pair_simulator(self, tree_obj, rho_site, L, delta, k)
        elif type == "seq":
            seq_simulator(self, tree_obj, rho_site, L, delta, k)

    def __repr__(self):
        return f"ARG(n={self.n}, rho={self.rho}, L={self.L}, delta={self.delta})"
    
    def to_dict(self):
        """Returns the ARG as a dictionary, similar to the R list output."""
        return {
            "type": self.type,
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
