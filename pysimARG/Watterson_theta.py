import numpy as np


def Watterson_theta(mat, S):
    """
    Watterson's theta estimator.
    
    Compute Watterson's theta estimator for a given incidence matrix.
    
    Parameters
    ----------
    mat : np.ndarray
        An incidence matrix where each row represents a sequence.
    
    S : int
        The number of segregating sites.

    Returns
    -------
    float
        Watterson's theta estimator.
    """
    n = mat.shape[0]
    
    a_n = np.sum(1 / np.arange(1, n))
    
    theta = S / a_n
    
    return theta
