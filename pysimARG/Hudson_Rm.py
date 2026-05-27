import numpy as np
from G4_test import G4_test


def Hudson_Rm(mat):
    """
    Hudson's Rm statistic.
    
    Compute Hudson's Rm statistic for a given incidence matrix.
    
    Parameters
    ----------
    mat : np.ndarray
        An incidence matrix where each row represents a sequence.

    Returns
    -------
    int
        Hudson's Rm statistic.
    """
    n_sites = mat.shape[1]
    
    incompatible_intervals = []
    
    for i in range(n_sites - 1):
        for j in range(i + 1, n_sites):
            
            if G4_test(mat[:, [i, j]]):
                incompatible_intervals.append((i, j))

    if not incompatible_intervals:
        return 0

    incompatible_intervals.sort(key=lambda x: x[1])
    
    rm_count = 0
    last_cut = -1
    
    for start, end in incompatible_intervals:
        if start > last_cut:
            last_cut = end - 1
            rm_count += 1
            
    return rm_count
