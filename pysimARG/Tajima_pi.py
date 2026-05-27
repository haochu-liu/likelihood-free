import numpy as np


def Tajima_pi(mat):
    """
    Tajima's pi estimator.
    
    Compute Tajima's pi estimator for a given incidence matrix.
    
    Parameters
    ----------
    mat : np.ndarray
        An incidence matrix where each row represents a sequence.
    
    Returns
    -------
    float
        Tajima's pi estimator.
    """
    n = mat.shape[0]

    pi_sum = 0
    for i in range(n):
        for j in range(i + 1, n):
            pi_sum += np.sum(mat[i] != mat[j])
    
    pi = 2 * pi_sum / (n * (n - 1))
    return pi
