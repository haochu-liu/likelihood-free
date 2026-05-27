import numpy as np


def Tajima_pi(mat, Wakeley=False):
    """
    Tajima's pi estimator.
    
    Compute Tajima's pi estimator for a given incidence matrix.
    And optionally, compute Wakeley's pi squared estimator if `Wakeley` is set to True.
    
    Parameters
    ----------
    mat : np.ndarray
        An incidence matrix where each row represents a sequence.
    
    Wakeley : bool, optional
        Whether to compute Wakeley's pi squared estimator, by default False

    Returns
    -------
    dict
        Tajima's pi estimator and optionally Wakeley's pi squared estimator.
    """
    n = mat.shape[0]

    pi_list = []
    for i in range(n-1):
        for j in range(i + 1, n):
            pi_list.append(np.sum(mat[i] != mat[j]))

    pi = 2 * np.sum(pi_list) / (n * (n - 1))
    result = {'pi': pi}

    if Wakeley:
        pi2 = 2 * np.sum([(x - pi)**2 for x in pi_list]) / (n * (n - 1))
        result['pi2'] = pi2

    return result
