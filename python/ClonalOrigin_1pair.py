import sys
import multiprocessing
import concurrent.futures
import numpy as np
import scipy.stats as stats
from pathlib import Path

# --- Path Setup & Imports ---
current_path = Path(__file__).resolve()
project_root = current_path.parent.parent
pysimARG_path = project_root / 'pysimARG'
data_path = project_root / 'data'

sys.path.append(str(pysimARG_path))
from clonal_genealogy import ClonalTree
from ClonalOrigin_pair import ARG
from add_mutation_truncated import add_mutation_truncated
from G3_test import G3_test
from LD_r import LD_r
from homoplasy_index import homoplasy_index

# --- Worker Functions ---
def generate_arg(seed, tree, rho_site, L, delta, k, theta_site):
    """Phase 1: Generates the ARG and computes its log weight."""
    # Crucial: Set a unique seed per worker so they don't generate identical random streams
    np.random.seed(seed)
    
    ARG_i = ARG(tree, rho_site, L, delta, k)
    arg_1_length = np.sum(ARG_i.edge[ARG_i.edge_mat[:, 0] == 1, 2])
    arg_2_length = np.sum(ARG_i.edge[ARG_i.edge_mat[:, 1] == 1, 2])
    log_weight = stats.expon.logcdf(arg_1_length, scale=1/(theta_site/2)) + \
                 stats.expon.logcdf(arg_2_length, scale=1/(theta_site/2))
                 
    return ARG_i, log_weight

def compute_stats(seed, arg_obj, theta_site, tree_width):
    """Phase 2: Adds mutations and calculates statistics."""
    np.random.seed(seed)
    
    node_site = add_mutation_truncated(arg_obj, theta_site)
    mat = node_site[:tree_width, :]
    
    r_val = LD_r(mat)
    g3_val = G3_test(mat)
    h_val = homoplasy_index(arg_obj, node_site)
    s_val = np.mean(np.any(mat, axis=0) & ~np.all(mat, axis=0))
    
    return r_val, g3_val, h_val, s_val

# --- Main Execution Block ---
if __name__ == '__main__':
    rho_site = 0.02
    theta_site = 0.05
    L = 1000000
    delta = 300
    N = 500
    k_vec = [50, 200, 2000]

    np.random.seed(100)
    tree = ClonalTree(n=15)
    tree_width = tree.n

    # Initialize result arrays
    v_h = np.full(N * 3, np.nan)
    v_s = np.full(N * 3, np.nan)
    v_r = np.zeros((N, 3), dtype=np.float64)
    v_g3 = np.zeros((N, 3), dtype=np.float64)

    # Determine core count
    num_cores = 10
    print(f"Starting parallel execution on {num_cores} cores...")

    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
        for j in range(3):
            print(f"--- Processing k_vec {k_vec[j]} (Iteration {j+1}/3) ---")
            
            # Phase 1: Submit ARG generation tasks
            futures_arg = []
            for i in range(N):
                # Create a unique, reproducible seed for each task
                task_seed = 1000 + (j * N) + i 
                futures_arg.append(
                    executor.submit(generate_arg, task_seed, tree, rho_site, L, delta, k_vec[j], theta_site)
                )
            
            # Gather Phase 1 results
            arg_list = []
            log_weight_list = []
            for f in futures_arg:
                arg, lw = f.result()
                arg_list.append(arg)
                log_weight_list.append(lw)
            
            # Synchronization Step: Compute probabilities and sample
            log_weight_list = np.array(log_weight_list)
            shifted_log_p = log_weight_list - np.max(log_weight_list)
            weights = np.exp(shifted_log_p)
            probs = weights / np.sum(weights)

            sampled_values = np.random.choice(N, size=N, replace=True, p=probs)
            
            # Phase 2: Submit Statistics tasks using the sampled ARGs
            futures_stats = []
            for i in range(N):
                task_seed = 5000 + (j * N) + i
                sampled_arg = arg_list[sampled_values[i]]
                futures_stats.append(
                    executor.submit(compute_stats, task_seed, sampled_arg, theta_site, tree_width)
                )
            
            # Gather Phase 2 results and populate arrays
            for i, f in enumerate(futures_stats):
                r_val, g3_val, h_val, s_val = f.result()
                
                v_r[i, j] = r_val
                v_g3[i, j] = g3_val
                v_h[i + j * N] = h_val
                v_s[i + j * N] = s_val

    # --- Save Data ---
    print("Saving files to data/ClonalOrigin directory...")
    output_file = data_path / 'ClonalOrigin'
    output_file.mkdir(parents=True, exist_ok=True) # Ensure directory exists

    np.savetxt(output_file / 'v_r.csv', v_r, delimiter=",")
    np.savetxt(output_file / 'v_g3.csv', v_g3, delimiter=",")
    np.savetxt(output_file / 'v_h.csv', v_h, delimiter=",")
    np.savetxt(output_file / 'v_s.csv', v_s, delimiter=",")
    
    print("Simulation complete!")
