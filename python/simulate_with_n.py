import numpy as np
import torch
import concurrent.futures
import os
from tqdm import tqdm
import sys
from pathlib import Path
current_path = Path(__file__).resolve()
project_root = current_path.parent.parent
pysimARG_path = project_root / 'pysimARG'
data_path = project_root / 'data'

sys.path.append(str(pysimARG_path))
from clonal_genealogy import ClonalTree
from ClonalOrigin_simulator import ClonalOrigin_simulator
from ClonalOrigin_sim_n import ClonalOrigin_sim_n

torch_device = "cpu"


np.random.seed(100)
tree = ClonalTree(n=15)

rho_site = 0.02
theta_site = 0.05
L = 1000000
delta = 300

x_o = ClonalOrigin_simulator(tree, rho_site, theta_site, L, delta, N=2000, k_vec=[50, 200, 2000])
x_o = torch.tensor(x_o, device=torch_device)
x_o = x_o.flatten()
x_o_numpy = x_o.cpu().numpy()

def prior_simulator():
    theta = [None] * 2
    theta[0] = np.random.uniform(0.0, 0.2)
    theta[1] = np.random.uniform(0.0, 0.2)
    N = np.random.randint(100, 2001, size=10)
    return theta, N

def run_single_simulation(i, seed, tree, L, delta):
    np.random.seed(seed)

    theta, N_vec = prior_simulator()
    s_mat = ClonalOrigin_sim_n(tree, theta[0], theta[1], L, delta, N=N_vec, k_vec=[50, 200, 2000])

    return i, theta, N_vec, s_mat


if __name__ == "__main__":
    # Run simulations
    simulation_budget = 100
    save_interval = 20

    theta_mat = np.zeros((simulation_budget * 10, 3))
    x_mat = np.zeros((simulation_budget * 10, 9))

    seeds = np.random.randint(0, 2**31 - 1, size=simulation_budget)

    output_file = data_path
    
    num_cores = 10
    print(f"Starting parallel execution on {num_cores} cores...")
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = []
        for i in range(simulation_budget):
            future = executor.submit(run_single_simulation, i, seeds[i], tree, L, delta)
            futures.append(future)

        completed_count = 0
        for future in tqdm(concurrent.futures.as_completed(futures), total=simulation_budget, desc="Simulation Progress"):
            i, theta, N_vec, s_mat = future.result()

            theta_mat[i*10:(i+1)*10, :2] = theta
            theta_mat[i*10:(i+1)*10, 2] = N_vec
            x_mat[i*10:(i+1)*10, :] = s_mat
            
            completed_count += 1

            if completed_count % save_interval == 0:
                np.savez(output_file / 'simulation_checkpoint_latest.npz', 
                         theta=theta_mat, 
                         x=x_mat, 
                         completed=completed_count)

    np.savez(output_file / 'simulation_final_results.npz', theta=theta_mat, x=x_mat)
    print("All simulations finished and saved successfully!")
