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
        A 46-dimensional vector as the summary statistics of simulations.
    """
    tree_width = tree.n
    if tree_width != seg_mat.shape[0]:
        raise ValueError("Inconsistent number of sequences between tree and genetic matrix!")
    L = seg_mat.shape[1]

    s_vec = np.full(46, np.nan).astype(float)

    # Identify segregating sites
    has_true = seg_mat.any(axis=0)
    has_false = ~seg_mat.all(axis=0)
    idx_seg = np.where(has_true & has_false)[0]

    # Summary statistics LD r and G4 test
    seg_near, seg_far = 0, 0
    seg_20_50, seg_50_80 = 0, 0
    D_near, D_far, D_prime_near, D_prime_far, r2_near, r2_far = 0, 0, 0, 0, 0, 0
    g4_near, g4_far = 0, 0
    D_20_50, D_50_80, D_prime_20_50, D_prime_50_80, r2_20_50, r2_50_80 = 0, 0, 0, 0, 0, 0
    g4_20_50, g4_50_80 = 0, 0

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
                    D_near += LD_result['D']
                    D_prime_near += LD_result['D_prime']
                    r2_near += LD_result['r_square']
                    g4_near += g4_result
                    seg_near += 1
                else:
                    D_far += LD_result['D']
                    D_prime_far += LD_result['D_prime']
                    r2_far += LD_result['r_square']
                    g4_far += g4_result
                    seg_far += 1
                if 20 <= dist_ij < 50:
                    D_20_50 += LD_result['D']
                    D_prime_20_50 += LD_result['D_prime']
                    r2_20_50 += LD_result['r_square']
                    g4_20_50 += g4_result
                    seg_20_50 += 1
                if 50 <= dist_ij <= 80:
                    D_50_80 += LD_result['D']
                    D_prime_50_80 += LD_result['D_prime']
                    r2_50_80 += LD_result['r_square']
                    g4_50_80 += g4_result
                    seg_50_80 += 1
        
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

        s_vec[16] = D_20_50
        s_vec[17] = D_50_80
        s_vec[18] = D_prime_20_50
        s_vec[19] = D_prime_50_80
        s_vec[20] = r2_20_50
        s_vec[21] = r2_50_80

        s_vec[22] = g4_20_50
        s_vec[23] = g4_50_80

        s_vec[24] = D_20_50 / seg_20_50 if seg_20_50 > 0 else 0
        s_vec[25] = D_50_80 / seg_50_80 if seg_50_80 > 0 else 0
        s_vec[26] = D_prime_20_50 / seg_20_50 if seg_20_50 > 0 else 0
        s_vec[27] = D_prime_50_80 / seg_50_80 if seg_50_80 > 0 else 0
        s_vec[28] = r2_20_50 / seg_20_50 if seg_20_50 > 0 else 0
        s_vec[29] = r2_50_80 / seg_50_80 if seg_50_80 > 0 else 0

        s_vec[30] = g4_20_50 / seg_20_50 if seg_20_50 > 0 else 0
        s_vec[31] = g4_50_80 / seg_50_80 if seg_50_80 > 0 else 0
    else:
        s_vec[:32] = 0
        s_vec[42] = 0 # Kelly's Zns estimator
    
    # Summary statistic homoplasy index
    s_vec[32] = homoplasy_index(tree, seg_mat)

    # Summary statistic clade homoplasy index
    # s_vec[k] = clade_homoplasy(tree, node_site)

    # Summary statistic proportion of segregating sites
    count_S = idx_seg.size
    s_vec[33] = count_S / L

    # Watterson's theta estimator
    s_vec[34] = Watterson_theta(seg_mat, count_S)

    # Tajima's pi estimators
    tajima_dict = Tajima_pi(seg_mat, Wakeley=True)
    s_vec[35] = tajima_dict['pi']
    s_vec[36] = tajima_dict['pi2']

    # Tajima's D statistic
    s_vec[37] = Tajima_D(seg_mat, s_vec[35], s_vec[34], count_S)

    # Wall's B and Q statistics
    wall_dict = Wall_BQ(seg_mat[:, idx_seg])
    s_vec[38] = wall_dict['B']
    s_vec[39] = wall_dict['Q']

    # Hudson's Rm estimator
    if not incompatible_intervals:
        s_vec[40] = 0

    incompatible_intervals.sort(key=lambda x: x[1])
    
    rm_count = 0
    last_cut = -1
    
    for start, end in incompatible_intervals:
        if start > last_cut:
            last_cut = end - 1
            rm_count += 1

    s_vec[40] = rm_count

    # Kelly's Zns estimator
    s_vec[41] = np.mean(r_squares)

    # Regression coefficient of r^2 on distance
    if len(r_squares) >= 2:
        df = pd.DataFrame({'x': distances, 'y': r_squares})
        mean_df = df.groupby('x')['y'].mean().reset_index()
        try:
            coeff = exp_regression(np.array(mean_df['x']), np.array(mean_df['y']))
            s_vec[42:45] = coeff
        except RuntimeError as e:
            s_vec[42:45] = 0
    else:
        s_vec[42:45] = 0

    # Add the length of sequence as a summary statistic
    s_vec[45] = L
    
    return s_vec
