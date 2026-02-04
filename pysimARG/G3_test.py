import numpy as np


def G3_test(mat):
    """
    Three-gamete test.
    
    Compute the number of combinations that have patterns 01, 10, 11.
    
    Parameters
    ----------
    mat : np.ndarray
        An incidence matrix with two columns.
    
    Returns
    -------
    bool
        True if all three patterns (01, 10, 11) are present, False otherwise.
    """
    if mat.shape[1] != 2:
        raise ValueError("`mat` must have two columns.")
    if mat.shape[0] < 3:
        raise ValueError("`mat` must have at least 3 rows.")
    
    # Count pattern 01: first column False, second column True
    num01 = np.sum(~mat[:, 0] & mat[:, 1])
    # Count pattern 10: first column True, second column False
    num10 = np.sum(mat[:, 0] & ~mat[:, 1])
    # Count pattern 11: both columns True
    num11 = np.sum(mat[:, 0] & mat[:, 1])
    
    # Check if all patterns are present (non-zero counts)
    if num01 > 0 and num10 > 0 and num11 > 0:
        return True
    else:
        return False
