import numpy as np


def LD(mat):
    """
    The square of correlation coefficient for LD.
    
    Compute the square of Pearson's correlation coefficient for linkage disequilibrium.
    
    Parameters
    ----------
    mat : np.ndarray
        An incidence matrix with two columns.
    
    Returns
    -------
    dict
        A dictionary containing the following keys:
        - 'D': The linkage disequilibrium coefficient.
        - 'D_prime': The normalized linkage disequilibrium coefficient.'
        - 'r_square': The square of Pearson's correlation coefficient for LD.
    """
    if mat.shape[1] != 2:
        raise ValueError("`mat` must have two columns.")
    
    n = mat.shape[0]
    
    p_A = np.sum(mat[:, 0]) / n
    p_a = 1 - p_A
    
    p_B = np.sum(mat[:, 1]) / n
    p_b = 1 - p_B
    
    p_AB = np.sum(mat[:, 0] & mat[:, 1]) / n
    
    D_AB = p_AB - (p_A * p_B)
    r_square = D_AB ** 2 / (p_A * p_a * p_B * p_b + np.finfo(float).eps)

    if D_AB > 0:
        D_prime = D_AB / (min(p_A * p_b, p_a * p_B) + np.finfo(float).eps)
    else:
        D_prime = D_AB / (-min(p_A * p_B, p_a * p_b) + np.finfo(float).eps)
    
    return {'D': D_AB, 'D_prime': D_prime, 'r_square': r_square}
