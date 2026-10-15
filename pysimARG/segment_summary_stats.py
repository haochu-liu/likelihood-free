import numpy as np
import pandas as pd
from G4_test import G4_test
from LD import LD
from homoplasy_index import homoplasy_index
from Watterson_theta import Watterson_theta
from Tajima_pi import Tajima_pi
from Tajima_D import Tajima_D
from Wall_BQ import Wall_BQ
from exp_regression import exp_regression
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)


def segment_summary_stats(tree, seg_mat):
    """
    Compute summary statistics for a given clonal tree and genetic matrix.
    
    Parameters
    ----------
    tree : ClonalTree
        The clonal genealogy.
    seg_mat : np.ndarray
        The genetic matrix.
    
    Returns
    -------
    np.ndarray
        A 30-dimensional vector as the summary statistics of simulations.
    """
    tree_width = tree.n
    if tree_width != seg_mat.shape[0]:
        raise ValueError("Inconsistent number of sequences between tree and genetic matrix!")
    L = seg_mat.shape[1]

    s_vec = np.full(30, np.nan).astype(float)

    # Identify segregating sites
    has_true = seg_mat.any(axis=0)
    has_false = ~seg_mat.all(axis=0)
    idx_seg = np.where(has_true & has_false)[0]

    # Summary statistics LD r and G4 test
    seg_near, seg_far = 0, 0
    D_near, D_far, D_prime_near, D_prime_far, r2_near, r2_far = 0, 0, 0, 0, 0, 0
    g4_near, g4_far = 0, 0

    r_squares = []
    distances = []
    incompatible_intervals = []

    if idx_seg.size >= 2:
        for i in range(idx_seg.size - 1):
            for j in range(i + 1, idx_seg.size):
                dist_ij = idx_seg[j] - idx_seg[i]
                idx_pair = [idx_seg[i], idx_seg[j]]

                # Compute LD statistics
                LD_result = LD(seg_mat[:, idx_pair])
                r_sq = LD_result['r_square']
                r_squares.append(r_sq)
                distances.append(dist_ij)

                # Compute G4 test
                g4_result = G4_test(seg_mat[:, idx_pair])
                if g4_result:
                    incompatible_intervals.append((idx_seg[i], idx_seg[j]))

                if dist_ij < L/2:
                    D_near += abs(LD_result['D'])
                    D_prime_near += LD_result['D_prime']
                    r2_near += LD_result['r_square']
                    g4_near += g4_result
                    seg_near += 1
                else:
                    D_far += abs(LD_result['D'])
                    D_prime_far += LD_result['D_prime']
                    r2_far += LD_result['r_square']
                    g4_far += g4_result
                    seg_far += 1
        
        s_vec[0] = D_near
        s_vec[1] = D_far
        s_vec[2] = D_prime_near
        s_vec[3] = D_prime_far
        s_vec[4] = r2_near
        s_vec[5] = r2_far

        s_vec[6] = g4_near
        s_vec[7] = g4_far

        s_vec[8] = D_near / seg_near if seg_near > 0 else 0
        s_vec[9] = D_far / seg_far if seg_far > 0 else 0
        s_vec[10] = D_prime_near / seg_near if seg_near > 0 else 0
        s_vec[11] = D_prime_far / seg_far if seg_far > 0 else 0
        s_vec[12] = r2_near / seg_near if seg_near > 0 else 0
        s_vec[13] = r2_far / seg_far if seg_far > 0 else 0

        s_vec[14] = g4_near / seg_near if seg_near > 0 else 0
        s_vec[15] = g4_far / seg_far if seg_far > 0 else 0

        # Kelly's Zns estimator
        s_vec[25] = np.mean(r_squares)
    else:
        s_vec[:16] = 0
        s_vec[25] = 0 # Kelly's Zns estimator
    
    # Summary statistic homoplasy index
    s_vec[16] = homoplasy_index(tree, seg_mat)

    # Summary statistic proportion of segregating sites
    count_S = idx_seg.size
    s_vec[17] = count_S / L

    # Watterson's theta estimator
    s_vec[18] = Watterson_theta(seg_mat, count_S)

    # Tajima's pi estimators
    tajima_dict = Tajima_pi(seg_mat, Wakeley=True)
    s_vec[19] = tajima_dict['pi']
    s_vec[20] = tajima_dict['pi2']

    # Tajima's D statistic
    s_vec[21] = Tajima_D(seg_mat, s_vec[19], s_vec[18], count_S)

    # Wall's B and Q statistics
    wall_dict = Wall_BQ(seg_mat[:, idx_seg])
    s_vec[22] = wall_dict['B']
    s_vec[23] = wall_dict['Q']

    # Hudson's Rm estimator
    if not incompatible_intervals:
        s_vec[24] = 0

    incompatible_intervals.sort(key=lambda x: x[1])
    
    rm_count = 0
    last_cut = -1
    
    for start, end in incompatible_intervals:
        if start > last_cut:
            last_cut = end - 1
            rm_count += 1

    s_vec[24] = rm_count

    # Regression coefficient of r^2 on distance
    if len(r_squares) >= 4: # Ensure enough data points for regression
        df = pd.DataFrame({'x': distances, 'y': r_squares})
        mean_df = df.groupby('x')['y'].mean().reset_index()
        try:
            coeff = exp_regression(np.array(mean_df['x']), np.array(mean_df['y']))
            s_vec[26:29] = coeff
        except RuntimeError as e:
            s_vec[26:29] = 0
    else:
        s_vec[26:29] = 0

    # Add the length of sequence as a summary statistic
    s_vec[29] = L
    
    return s_vec
