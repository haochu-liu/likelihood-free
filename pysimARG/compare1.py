import numpy as np

def clonal_origin_nodes(matrix_data, n):
    """
    Find the row index for recombination nodes on a given edge.

    Args:
        matrix_data (np.ndarray): A recomb_edge matrix.
        n (int): The edge index.

    Returns:
        np.ndarray: A vector of signed indices indicating where the nodes are.
                    Positive indicates Type A match, Negative indicates Type B match.
                    (Returns 1-based indices to preserve sign distinction for index 0).
    """
    # Ensure input is a numpy array
    matrix_data = np.array(matrix_data)

    # --- Column Mapping ---
    # R uses 1-based column indices; Python uses 0-based.
    # R col 1 -> Python col 0
    # R col 2 -> Python col 1
    # R col 3 -> Python col 2
    # R col 4 -> Python col 3

    # Find rows for Type A: where matrix_data[, 3] == n (Python col 2)
    rows_type_a = np.where(matrix_data[:, 2] == n)[0]

    # Find rows for Type B: where matrix_data[, 1] == n (Python col 0)
    rows_type_b = np.where(matrix_data[:, 0] == n)[0]

    # Return empty if no nodes found
    if len(rows_type_a) == 0 and len(rows_type_b) == 0:
        return np.array([])

    found_nodes_list = []

    # --- Process Type A ---
    # R Logic: found_nodes[, 1] <- rows_type_a
    #          found_nodes[, 2] <- matrix_data[rows_type_a, 4]
    if len(rows_type_a) > 0:
        # Values from Python col 3 (R col 4)
        vals_a = matrix_data[rows_type_a, 3]
        for idx, val in zip(rows_type_a, vals_a):
            # Using idx + 1 to make it 1-based
            found_nodes_list.append([idx + 1, val])

    # --- Process Type B ---
    # R Logic: found_nodes[, 1] <- -rows_type_b
    #          found_nodes[, 2] <- matrix_data[rows_type_b, 2]
    if len(rows_type_b) > 0:
        # Values from Python col 1 (R col 2)
        vals_b = matrix_data[rows_type_b, 1]
        for idx, val in zip(rows_type_b, vals_b):
            # Using -(idx + 1) to make it 1-based and negative
            found_nodes_list.append([-(idx + 1), val])

    # Convert to numpy array for operations
    found_nodes = np.array(found_nodes_list)

    # --- Sort ---
    # Sort by the second column (value/position) if there is more than 1 row
    if len(found_nodes) > 1:
        # argsort returns the indices that would sort the array based on column 1
        sort_order = np.argsort(found_nodes[:, 1])
        found_nodes = found_nodes[sort_order]

    # Return the first column (the signed indices)
    return found_nodes[:, 0]
