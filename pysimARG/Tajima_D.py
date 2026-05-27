import numpy as np


def Tajima_D(mat, pi_hat, theta_hat, S):
    """
    Tajima's D estimator.
    
    Compute Tajima's D estimator for a given incidence matrix.
    
    Parameters
    ----------
    mat : np.ndarray
        An incidence matrix where each row represents a sequence.
    
    pi_hat : float
        Tajima's pi estimator.
    
    theta_hat : float
        Watterson's theta estimator.

    S : int
        The number of segregating sites.
    
    Returns
    -------
    float
        Tajima's D estimator.
    """
    n = mat.shape[0]
    an = np.sum(1 / np.arange(1, n))
    bn = np.sum(1 / np.arange(1, n) ** 2)

    e1 = (n + 1) / (3 * an * (n - 1)) - 1 / an**2
    e2 = (2 * (n**2 + n + 3) / (9 * n * (n - 1)) - (n + 2) / (n * an) + bn / an**2) / (an**2 + bn)

    D = (pi_hat - theta_hat) / np.sqrt(e1 * S + e2 * S * (S - 1))
    return D
