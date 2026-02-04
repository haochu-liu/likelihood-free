import numpy as np


def ClonalOrigin_nodes(matrix_data, n):
    """
    Find the row index for recombination nodes on a given edge.
    
    Parameters
    ----------
    matrix_data : np.ndarray
        A recomb_edge matrix.
    n : int
        The edge index.
    
    Returns
    -------
    np.ndarray
        Vector of where the nodes are. Positive values indicate rows where
        column 3 matches n (type a), negative values indicate rows where
        column 1 matches n (type b).
    """
    # Find rows where column 3 (0-indexed: 2) equals n
    rows_type_a = np.where(matrix_data[:, 2] == n)[0]
    # Find rows where column 1 (0-indexed: 0) equals n
    rows_type_b = np.where(matrix_data[:, 0] == n)[0]
    
    total_rows = len(rows_type_a) + len(rows_type_b)
    
    if total_rows == 0:
        return np.array([], dtype=float)
    
    found_nodes = np.full((total_rows, 2), np.nan)
    
    if len(rows_type_a) > 0:
        # R is 1-indexed, Python is 0-indexed
        # Store 1-indexed row numbers to match R behavior
        found_nodes[:len(rows_type_a), 0] = rows_type_a + 1
        found_nodes[:len(rows_type_a), 1] = matrix_data[rows_type_a, 3]
    
    if len(rows_type_b) > 0:
        start_idx = len(rows_type_a)
        end_idx = start_idx + len(rows_type_b)
        # Negative 1-indexed row numbers to match R behavior
        found_nodes[start_idx:end_idx, 0] = -(rows_type_b + 1)
        found_nodes[start_idx:end_idx, 1] = matrix_data[rows_type_b, 1]
    
    if found_nodes.shape[0] > 1:
        # Sort by second column
        found_nodes = found_nodes[found_nodes[:, 1].argsort()]
    
    return found_nodes[:, 0]
