import numpy as np
from LD import LD


def Kelly_Z(mat):
    """
    Kelly's Z estimator.
    
    Compute Kelly's Z estimator for a given incidence matrix.
    
    Parameters
    ----------
    mat : np.ndarray
        An incidence matrix where each row represents a sequence.

    Returns
    -------
    float
        Kelly's Z estimator.
    """
    has_true = mat.any(axis=0)
    has_false = ~mat.all(axis=0)
    idx_seg = np.where(has_true & has_false)[0]

    r_squares = []
    if idx_seg.size >= 2:
        for i in range(idx_seg.size - 1):
            for j in range(i + 1, idx_seg.size):
                idx_pair = [idx_seg[i], idx_seg[j]]
                r_sq = LD(mat[:, idx_pair])['r_square']
                r_squares.append(r_sq)
    
    if len(r_squares) == 0:
        return 0.0
    
    return np.mean(r_squares)
