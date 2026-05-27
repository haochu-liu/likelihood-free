import numpy as np


def Wall_BQ(seg_mat):
    """
    Wall's BQ estimator.
    
    Compute Wall's B and Q estimators for a given incidence matrix.
    
    Parameters
    ----------
    seg_mat : np.ndarray
        An incidence matrix where each row represents a sequence.
        And each column represents a segregating site.

    Returns
    -------
    dict
        Wall's B and Q estimators in a dictionary.
    """
    seg_mat = seg_mat.astype(bool)
    n = seg_mat.shape[0]
    S = seg_mat.shape[1]

    if S < 2:
        return {'B': 0, 'Q': 0}

    B = 0
    A = set()
    for i in range(S - 1):
        if seg_mat[0, i]:
            col1 = np.logical_not(seg_mat[:, i])
        else:
            col1 = seg_mat[:, i]
        
        if seg_mat[0, i + 1]:
            col2 = np.logical_not(seg_mat[:, i + 1])
        else:
            col2 = seg_mat[:, i + 1]
        
        if np.array_equal(col1, col2):
            B += 1
            A.add(tuple(col1))
    
    B_norm = B / (S - 1)
    Q = (B + len(A)) / S
    return {'B': B_norm, 'Q': Q}
