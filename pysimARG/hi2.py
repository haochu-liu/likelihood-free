import numpy as np

def homoplasy_index(ARG):
    # 1. Class validation
    # Assuming FSM_ARG is a custom class. In Python, we usually check type directly.
    if type(ARG).__name__ != "FSM_ARG":
        raise ValueError("Object must be of class 'FSM_ARG'")

    # Initialize dimensions based on the NumPy matrices
    n_node = ARG.node_mat.shape[0]
    n_leaf = ARG.n
    n_site = ARG.node_mat.shape[1]
    
    # Initialize vectors
    m_vec = np.zeros(n_site)
    s_vec = np.full(n_site, np.nan)

    for site_loc in range(n_site):
        # Assuming local_tree_FSM is defined elsewhere in your Python script
        ARG_site = local_tree_FSM(ARG, site_loc)
        
        # Flatten edge[:, 0:2], get unique values, and sort
        # (Equivalent to sort(unique(as.vector(ARG_site$edge[, 1:2]))))
        node_vec = np.sort(np.unique(ARG_site.edge[:, 0:2].flatten()))
        
        # Note: We assume node_vec contains 0-based integer indices suitable for slicing
        node_site_vec = ARG_site.node_site[node_vec, site_loc]

        # Compute minimum possible changes
        leaf_states = node_site_vec[:n_leaf]
        if np.any(leaf_states) and not np.all(leaf_states):
            m_vec[site_loc] = 1

        # Compute actual changes
        s_site = 0
        
        # Use a dictionary instead of R's named list. 
        # We store the leaf states as Python 'sets' to make intersections/unions easy.
        site_dict = {
            node: {int(state)} for node, state in zip(node_vec[:n_leaf], leaf_states)
        }

        # Iterate through the remaining internal nodes
        # Python range is exclusive at the end, equivalent to (n_leaf+1):length(node_vec)
        for i in range(n_leaf, len(node_vec)):
            parent_node = node_vec[i]
            
            # Find indices where the first column (parent column) equals parent_node
            node_index = np.where(ARG_site.edge[:, 0] == parent_node)[0]
            
            if len(node_index) == 2:
                # Coalescent structure
                children_node = ARG_site.edge[node_index, 1]
                child_1 = site_dict[children_node[0]]
                child_2 = site_dict[children_node[1]]
                
                # Check set intersection
                intersec = child_1.intersection(child_2)
                
                if len(intersec) == 0:
                    # intersection is empty -> mutation
                    site_dict[parent_node] = child_1.union(child_2)
                    s_site += 1
                else:
                    # intersection is not empty -> no mutation
                    site_dict[parent_node] = intersec
                    
            elif len(node_index) == 1:
                # Recombination structure
                child_node = ARG_site.edge[node_index[0], 1]
                # Copy the set to avoid accidental reference mutation
                site_dict[parent_node] = set(site_dict[child_node])
                
        s_vec[site_loc] = s_site

    # Final calculation
    sum_m = np.sum(m_vec)
    sum_s = np.sum(s_vec)
    
    if sum_m == 0 and sum_s == 0:
        hi = 0.0
    else:
        hi = 1.0 - (sum_m / sum_s)
        
    return hi
