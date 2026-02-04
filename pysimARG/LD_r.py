import numpy as np


def LD_r(mat):
    """
    The square of correlation coefficient for LD.
    
    Compute the square of Pearson's correlation coefficient for linkage disequilibrium.
    
    Parameters
    ----------
    mat : np.ndarray
        An incidence matrix with two columns.
    
    Returns
    -------
    float
        r square value.
    """
    if mat.shape[1] != 2:
        raise ValueError("`mat` must have two columns.")
    
    n = mat.shape[0]
    
    n_A = np.sum(mat[:, 0])
    n_a = n - n_A
    
    n_B = np.sum(mat[:, 1])
    n_b = n - n_B
    
    n_AB = np.sum(mat[:, 0] & mat[:, 1])
    
    D_AB = n_AB / n - (n_A * n_B) / (n ** 2)
    r_square = D_AB ** 2 * n ** 4 / (n_A * n_a * n_B * n_b + np.finfo(float).eps)
    
    return r_square
