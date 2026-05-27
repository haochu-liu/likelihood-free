import numpy as np


def Watterson_theta(mat):
    """
    Watterson's theta estimator.
    
    Compute Watterson's theta estimator for a given incidence matrix.
    
    Parameters
    ----------
    mat : np.ndarray
        An incidence matrix where each row represents a sequence.
    
    Returns
    -------
    float
        Watterson's theta estimator.
    """
    n = mat.shape[0]
    
    has_true = mat.any(axis=0)
    has_false = ~mat.all(axis=0)
        
    S = np.sum(has_true & has_false)
    
    a_n = np.sum(1 / np.arange(1, n))
    
    theta = S / a_n
    
    return theta
