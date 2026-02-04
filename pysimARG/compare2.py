import numpy as np


def clonal_origin_pair_seq(tree, rho_site, L, delta, k):
    """
    An approximated ancestral recombination graph (ARG) using ClonalOrigin algorithm.
    Simulates coalescent and recombination events at two sites.

    Args:
        tree (ClonalTree): The clonal genealogy object.
        rho_site (float): Recombination parameter per site.
        L (int): Number of sites.
        delta (float): Mean of recombinant segment length.
        k (int): Distance between two sites.

    Returns:
        dict: A dictionary containing the ARG structure (edges, nodes, etc.)
    """
    
    # --- Input Validation ---
    # Assuming 'tree' is an instance of the ClonalTree class defined previously
    if not hasattr(tree, 'edge') or not hasattr(tree, 'node_height'):
        raise ValueError("Object must be a valid ClonalTree instance.")
    if not isinstance(L, int):
        raise ValueError("`L` must be a single integer!")
    if not isinstance(k, int):
        raise ValueError("`k` must be a single integer!")
    if delta <= 0:
        raise ValueError("`delta` must be greater than zero!")

    # --- Initialization ---
    rho = L * rho_site
    clonal_edge = tree.edge
    clonal_node_height = tree.node_height
    n = tree.n

    # Number of recombination edges
    # clonal_edge column 2 is length
    tree_length = np.sum(clonal_edge[:, 2])

    # Initialize recombination edges matrix
    # Columns: 0:b_edge, 1:b_height, 2:a_edge, 3:a_height, 4:x, 5:y
    nrow_max = 1000
    recomb_edge = np.full((nrow_max, 6), np.nan) 

    n_recomb = 0

    # --- Add Recombination Sequentially (Sites 1 and 2) ---
    for i in range(1, 3): # Loops i = 1, 2
        R_new = 0
        R_old = 0
        remain_index = []

        if i == 1:
            # Poisson sample
            rate = rho_site * delta * tree_length / 2.0
            R_new = np.random.poisson(rate)
        else: # i == 2
            # Find indices where site y (col 5) == 1. 
            # Note: We store sites as 1 and 2 to match R logic for now.
            survive_index = np.where(recomb_edge[:n_recomb, 5] == 1)[0]
            
            # Sum calculation: sum((1 - 1/delta)**k for k in 0..k-1)
            # This is a geometric series sum
            exponents = np.arange(k)
            delta2 = np.sum((1.0 - 1.0/delta)**exponents)
            
            rate = rho_site * delta2 * tree_length / 2.0
            R_new = np.random.poisson(rate)

            if len(survive_index) > 0:
                prob = (1.0 - 1.0/delta)**k
                R_old = np.random.binomial(len(survive_index), prob)
                
                if len(survive_index) == 1:
                    if R_old == 1:
                        remain_index = survive_index
                    else:
                        remain_index = []
                else:
                    # Sample R_old elements from survive_index
                    remain_index = np.random.choice(survive_index, size=R_old, replace=False)
            else:
                R_old = 0

        # Process New Recombinations
        if R_new > 0:
            # Expand matrix if needed
            if n_recomb + R_new >= nrow_max:
                extension = np.full((nrow_max, 6), np.nan)
                recomb_edge = np.vstack((recomb_edge, extension))
                nrow_max *= 2

            # Set site x and y (cols 4, 5) to i
            start_idx = n_recomb
            end_idx = n_recomb + R_new
            recomb_edge[start_idx:end_idx, 4] = i
            recomb_edge[start_idx:end_idx, 5] = i

            # Sample exponential waiting times for 'above' nodes
            a_rexp = np.random.exponential(scale=1.0, size=R_new)

            # Simulate b_edge (similar to mutation)
            # Sample edge indices based on length probability
            edge_indices = np.arange(2 * (n - 1))
            probs = clonal_edge[:, 2] / np.sum(clonal_edge[:, 2])
            selected_b_edges = np.random.choice(edge_indices, size=R_new, replace=True, p=probs)
            recomb_edge[start_idx:end_idx, 0] = selected_b_edges

            # Loop through each new event to calculate heights
            for j in range(R_new):
                curr_idx = start_idx + j
                b_edge_idx = int(recomb_edge[curr_idx, 0])
                
                # Simulate b_height: Uniform along the edge
                # Height = uniform(0, len) + child_node_height
                edge_len = clonal_edge[b_edge_idx, 2]
                child_node = int(clonal_edge[b_edge_idx, 1])
                b_height = np.random.uniform(0, edge_len) + clonal_node_height[child_node]
                recomb_edge[curr_idx, 1] = b_height

                # Identify a_height (Coalescent simulation upwards)
                # Internal node heights: indices n to 2n-2
                internal_heights = clonal_node_height[n : 2*n-1]
                t_above_b = internal_heights - b_height
                
                # Filter positive time differences and sort (implicitly sorted by node index usually)
                valid_t = np.sort(t_above_b[t_above_b >= 0])
                
                # Create intervals
                i_above_b = np.concatenate(([0], valid_t))
                i_above_b = np.diff(i_above_b) # differences between adjacent elements

                # Weights logic: (1 + num_intervals) down to 2
                weights = np.arange(len(i_above_b) + 1, 1, -1)
                cuml_above_b = np.cumsum(i_above_b * weights)

                # Determine number of lineages
                # Count how many cuml steps are less than random value a_rexp[j]
                survived_steps = np.where(a_rexp[j] > cuml_above_b)[0]
                num_lineage = (1 + len(i_above_b)) - len(survived_steps)

                if num_lineage == (1 + len(i_above_b)):
                     recomb_edge[curr_idx, 3] = (a_rexp[j] / num_lineage) + b_height
                else:
                    # Logic: (remaining_time / num_lineages) + time_offset + base_height
                    idx_cutoff = 1 + len(i_above_b) - num_lineage
                    # Python 0-based adjustment for indexing cuml_above_b
                    val_cuml = cuml_above_b[idx_cutoff - 1] 
                    sum_interval = np.sum(i_above_b[:idx_cutoff])
                    
                    recomb_edge[curr_idx, 3] = ((a_rexp[j] - val_cuml) / num_lineage) + sum_interval + b_height

                # Simulate a_edge
                if num_lineage > 1:
                    # Find edges that span the calculated a_height
                    # parent_height >= a_height AND child_height < a_height
                    parent_indices = clonal_edge[:, 0].astype(int)
                    child_indices = clonal_edge[:, 1].astype(int)
                    
                    mask = (clonal_node_height[parent_indices] >= recomb_edge[curr_idx, 3]) & \
                           (clonal_node_height[child_indices] < recomb_edge[curr_idx, 3])
                    
                    pool_edge = np.where(mask)[0]
                    
                    if len(pool_edge) > 0:
                        recomb_edge[curr_idx, 2] = np.random.choice(pool_edge)
                    else:
                        # Fallback (should theoretically not happen if logic holds)
                        recomb_edge[curr_idx, 2] = 2*n - 2 # Root edge approx
                else:
                    # Root
                    recomb_edge[curr_idx, 2] = 2*n - 2 # Since Python max index is 2n-2

        # Update Old Recombinations (Extend them to site i)
        if R_old > 0 and len(remain_index) > 0:
            recomb_edge[remain_index.astype(int), 5] = i

        n_recomb += R_new

    # --- Handling No Recombination ---
    if n_recomb == 0:
        node_mat = np.ones((2*n - 1, 2), dtype=bool)
        edge_mat = np.ones((2*(n - 1), 2), dtype=bool)
        return {
            "edge": clonal_edge,
            "edge_mat": edge_mat,
            "node_height": clonal_node_height,
            "node_mat": node_mat,
            "node_clonal": np.full(2*n - 1, True),
            "sum_time": np.max(clonal_node_height),
            "n": n, "rho": rho, "L": L, "delta": delta,
            "class": "FSM_ARG"
        }
    
    # Trim recomb_edge to actual size
    recomb_edge = recomb_edge[:n_recomb, :]

    # --- Recombination Segment and Ancestral Material ---
    node_max = (2*n - 1) + 3*n_recomb
    edge_max = 2*(n - 1) + 4*n_recomb
    
    # Initialize Main Output Matrices
    edge_matrix = np.full((edge_max, 3), np.nan)
    # edge_mat_index stores indices pointing to leaf nodes in node_info
    edge_mat_index = np.full(edge_max, -1, dtype=int) 
    
    node_mat = np.full((node_max, 2), False, dtype=bool)
    
    # node_info cols: 0:index, 1:height, 2:recomb, 3:clonal (0=False, 1=True)
    # Note: Using float for everything, but cols 0, 2, 3 act as ints/bools
    node_info = np.full((node_max, 4), np.nan)
    node_info[:, 0] = np.arange(node_max)

    # 1. Clonal Nodes (0 to 2n-2)
    node_mat[:n, :] = True # Leaves have all material
    node_info[:2*n-1, 1] = clonal_node_height
    node_info[:2*n-1, 2] = 0 # 0 indicates no active recomb event transition
    node_info[:2*n-1, 3] = 1 # True (Clonal)

    # 2. Recombination Nodes (Start/End points)
    # Range indices
    idx_start = 2*n - 1
    idx_mid = idx_start + 2*n_recomb
    
    # "Out" nodes (bottom of recomb edge) and "In" nodes? 
    # R: mapply with recomb_edge[, 2] (b_height)
    # We alternate two entries per recombination event
    b_heights = recomb_edge[:, 1]
    
    # Interleave heights: [h1, h1, h2, h2, ...]
    interleaved_heights = np.repeat(b_heights, 2)
    node_info[idx_start:idx_mid, 1] = interleaved_heights
    
    # Interleave recomb IDs: [-1, NA, -2, NA, ...]
    # R uses 1-based index for IDs. -1, -2 etc.
    # Python: We will use -1, -2... corresponding to recomb index + 1
    recomb_ids = -(np.arange(n_recomb) + 1)
    interleaved_ids = np.empty(2 * n_recomb)
    interleaved_ids[0::2] = recomb_ids
    interleaved_ids[1::2] = np.nan # NA in R
    node_info[idx_start:idx_mid, 2] = interleaved_ids
    
    # Interleave clonal flag: [True, False, True, False...]
    interleaved_clonal = np.empty(2 * n_recomb)
    interleaved_clonal[0::2] = 1 # True
    interleaved_clonal[1::2] = 0 # False
    node_info[idx_start:idx_mid, 3] = interleaved_clonal

    # 3. Recombination "Top" Nodes (a_height)
    idx_end = node_max
    node_info[idx_mid:idx_end, 1] = recomb_edge[:, 3] # a_height
    node_info[idx_mid:idx_end, 2] = np.arange(n_recomb) + 1 # Positive IDs 1..N
    node_info[idx_mid:idx_end, 3] = 1 # True

    # Sort node_info by height
    # np.argsort is stable enough usually, but height is float.
    order = np.argsort(node_info[:, 1])
    node_info = node_info[order]
    
    # Remap indices after sort because we need to find them later
    # We rely on 'index' column (col 0) to track original identities if needed,
    # but the loop below iterates strictly by the sorted order.

    # --- Recombination Nodes on Every Edge ---
    # We need the helper function 'clonal_origin_nodes' here.
    # It returns signed indices (1-based) matching the logic we used for 'recomb' col.
    recomb_node_map = []
    for edge_idx in range(2 * n - 1): # For every edge in ClonalTree
        # Pass raw recomb_edge matrix and edge index
        nodes = clonal_origin_nodes(recomb_edge, edge_idx)
        recomb_node_map.append(nodes)

    # --- Construct Full ARG Edges ---
    i = n # Start processing from first internal node (skipping leaves 0..n-1 in sorted list? No, R starts i = n+1, which is n in 0-based)
    # Wait, R: i = n + 1. But R indices are 1..N. n+1 is the first non-leaf if sorted purely by height (since leaves are 0 height).
    # In Python, leaves are 0..n-1. So we start at n.
    
    edge_idx_counter = 0

    while i < node_max:
        # Get info for current node in the sorted list
        curr_row = node_info[i]
        curr_node_original_idx = int(curr_row[0])
        curr_recomb_val = curr_row[2] # 0, <0, >0, or NaN
        curr_height = curr_row[1]

        # Case 1: Clonal Tree Node (recomb == 0)
        if curr_recomb_val == 0:
            node_index = curr_node_original_idx
            # Find rows in clonal_edge where parent is node_index
            leaf_edges = np.where(clonal_edge[:, 0] == node_index)[0]
            
            leaf_indices = []
            
            for edge_id in leaf_edges:
                # Check if there are recomb nodes on this edge
                r_nodes = recomb_node_map[edge_id]
                
                if len(r_nodes) > 0:
                    # target node is the last one in r_nodes (tail)
                    tar_node = r_nodes[-1]
                    # Find where this target node is in node_info's recomb column (col 2)
                    # Note: we need to find the index in the CURRENT sorted node_info
                    match_idx = np.where(node_info[:, 2] == tar_node)[0]
                    # R takes the first match usually? should be unique.
                    if len(match_idx) > 0:
                        leaf_indices.append(match_idx[0])
                    else:
                        raise ValueError(f"Target node {tar_node} not found in node_info")
                else:
                    # No recomb node, connects to original child
                    original_child = clonal_edge[edge_id, 1]
                    # Find original child in node_info (col 0)
                    match_idx = np.where(node_info[:, 0] == original_child)[0]
                    leaf_indices.append(match_idx[0])

            # Append edges
            # Set node1 (Parent) to i
            edge_matrix[edge_idx_counter : edge_idx_counter+2, 0] = i
            # Set node2 (Children)
            edge_matrix[edge_idx_counter : edge_idx_counter+2, 1] = leaf_indices
            # Set Length
            child_heights = node_info[leaf_indices, 1]
            edge_matrix[edge_idx_counter : edge_idx_counter+2, 2] = curr_height - child_heights
            
            # Store material index
            edge_mat_index[edge_idx_counter : edge_idx_counter+2] = leaf_indices

            # Update genetic material (Bitwise OR of children)
            node_mat[i, :] = node_mat[leaf_indices[0], :] | node_mat[leaf_indices[1], :]

            edge_idx_counter += 2
            i += 1

        # Case 2: Recombination Edge OUT Node (recomb < 0)
        # R check: node_info[i, 3] < 0. In numpy, NaN comparisons are false, so safe.
        elif not np.isnan(curr_recomb_val) and curr_recomb_val < 0:
            # This logic handles i and i+1 in R. 
            # In our sorted list, i and i+1 should be the paired nodes we created earlier
            # (interleaved -ID and NaN).
            
            recomb_id_abs = int(abs(curr_recomb_val)) # 1-based index
            recomb_row_idx = recomb_id_abs - 1 # 0-based index for recomb_edge
            
            # Identify leaf edge from recomb_edge (col 1 is b_edge in R -> col 0 in Py)
            leaf_edge_id = int(recomb_edge[recomb_row_idx, 0])
            
            # Find target node in the chain
            r_nodes = recomb_node_map[leaf_edge_id]
            # tar_node index in r_nodes array
            # which(recomb_node[[...]] == node_info[i, 3])
            tar_node_loc = np.where(r_nodes == curr_recomb_val)[0]
            
            leaf_node_original = -1
            if len(tar_node_loc) > 0 and tar_node_loc[0] == 0:
                # Corresponds to tar_node == 1 in R (1-based index 1 is 0)
                leaf_node_original = clonal_edge[leaf_edge_id, 1]
            elif len(tar_node_loc) > 0:
                # Look at previous node in chain
                prev_recomb_id = r_nodes[tar_node_loc[0] - 1]
                # Find that ID in node_info to get original index? 
                # R: node_info[which(...), 1] -> gets index
                match_row = np.where(node_info[:, 2] == prev_recomb_id)[0]
                leaf_node_original = node_info[match_row[0], 0]

            # Find this leaf node in current node_info list
            leaf_index = np.where(node_info[:, 0] == leaf_node_original)[0][0]

            # Append edges for both split nodes (i and i+1)
            # Both point to the same leaf_index
            edge_matrix[edge_idx_counter : edge_idx_counter+2, 0] = [i, i+1]
            edge_matrix[edge_idx_counter : edge_idx_counter+2, 1] = leaf_index
            edge_matrix[edge_idx_counter : edge_idx_counter+2, 2] = curr_height - node_info[leaf_index, 1]
            
            edge_mat_index[edge_idx_counter : edge_idx_counter+2] = [i, i+1]

            # Recombination sites x and y
            # x is col 4, y is col 5
            x = int(recomb_edge[recomb_row_idx, 4])
            y = int(recomb_edge[recomb_row_idx, 5])
            
            # R logic: 
            # node_mat[c(i, i+1), ] <- FALSE
            node_mat[i, :] = False
            node_mat[i+1, :] = False
            
            # R: x:y range. 
            # In R x:y is inclusive. In Python slice x-1 : y
            # However, code uses 1 and 2 for sites.
            # If x=1, y=1 -> slice 0:1. If x=2, y=2 -> slice 1:2.
            # Adjust to 0-based: x-1, y
            
            slice_start = x - 1
            slice_end = y 
            
            # Node i+1 inherits material IN the range
            node_mat[i+1, slice_start:slice_end] = node_mat[leaf_index, slice_start:slice_end]
            
            # Node i inherits material OUTSIDE the range
            # Python boolean indexing or manual mask
            mask = np.ones(2, dtype=bool)
            mask[slice_start:slice_end] = False
            node_mat[i, mask] = node_mat[leaf_index, mask]

            edge_idx_counter += 2
            i += 2 # Skip the pair

        # Case 3: Recombination Edge IN Node (recomb > 0)
        elif not np.isnan(curr_recomb_val) and curr_recomb_val > 0:
            recomb_id_abs = int(curr_recomb_val)
            recomb_row_idx = recomb_id_abs - 1
            
            leaf_edge_id = int(recomb_edge[recomb_row_idx, 2]) # col 3 in R (a_edge) -> col 2 in Py
            
            r_nodes = recomb_node_map[leaf_edge_id]
            tar_node_loc = np.where(r_nodes == curr_recomb_val)[0]
            
            leaf_node_original = -1
            if len(tar_node_loc) > 0 and tar_node_loc[0] == 0:
                 # tar_node == 1 in R
                 if leaf_edge_id == (2*n - 2): # Root edge logic (2n-1 in R)
                     leaf_node_original = 2*n - 2
                 else:
                     leaf_node_original = clonal_edge[leaf_edge_id, 1]
            elif len(tar_node_loc) > 0:
                prev_recomb_id = r_nodes[tar_node_loc[0] - 1]
                match_row = np.where(node_info[:, 2] == prev_recomb_id)[0]
                leaf_node_original = node_info[match_row[0], 0]

            leaf_index_1 = np.where(node_info[:, 0] == leaf_node_original)[0][0]
            
            # leaf_index_2: node_info[, 3] == -curr_recomb_val + 1 (the 'next' node in the pair)
            # Wait, logic in R: which(node_info[, 3] == -val) + 1.
            # In our sorted list, -val is at index X. We want X+1 (the blank partner).
            neg_match = np.where(node_info[:, 2] == -curr_recomb_val)[0]
            if len(neg_match) == 0:
                 raise ValueError(f"Could not find pair for {-curr_recomb_val}")
            leaf_index_2 = neg_match[0] + 1
            
            leaf_indices = [leaf_index_1, leaf_index_2]

            # Append edges
            edge_matrix[edge_idx_counter : edge_idx_counter+2, 0] = i
            edge_matrix[edge_idx_counter : edge_idx_counter+2, 1] = leaf_indices
            edge_matrix[edge_idx_counter : edge_idx_counter+2, 2] = curr_height - node_info[leaf_indices, 1]
            
            edge_mat_index[edge_idx_counter : edge_idx_counter+2] = leaf_indices

            # Update genetic material (Union)
            node_mat[i, :] = node_mat[leaf_indices[0], :] | node_mat[leaf_indices[1], :]

            edge_idx_counter += 2
            i += 1
        
        else:
            # Fallback for NA or weird states, simply increment
            i += 1

    # --- Final Output Construction ---
    # Slice edge_mat to valid entries based on index
    # We computed edge_mat_index to point to rows in node_mat
    valid_indices = edge_mat_index[:edge_idx_counter]
    final_edge_mat = node_mat[valid_indices, :]

    # Truncate edge_matrix
    edge_matrix = edge_matrix[:edge_idx_counter, :]

    arg_result = {
        "edge": edge_matrix,
        "edge_mat": final_edge_mat,
        "node_height": node_info[:, 1],
        "node_mat": node_mat,
        "node_clonal": node_info[:, 3] == 1,
        "sum_time": node_info[-1, 1], # Last node height (root)
        "n": n,
        "rho": rho,
        "L": L,
        "delta": delta,
        "class": "FSM_ARG"
    }
    
    return arg_result
