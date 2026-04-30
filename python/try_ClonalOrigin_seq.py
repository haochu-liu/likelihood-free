import numpy as np
from Bio import Phylo
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
import sys
import os
from pathlib import Path
current_path = Path(__file__).resolve()
project_root = current_path.parent.parent
pysimARG_path = project_root / 'pysimARG'
data_path = project_root / 'data'

sys.path.append(str(pysimARG_path))
from clonal_genealogy import ClonalTree
from ClonalOrigin_seq_sim import ClonalOrigin_seq_sim
from newick_to_tree import newick_to_tree

torch_device = "cpu"


np.random.seed(100)
tree = ClonalTree(n=15)

# Load phylo tree and convert to ClonalTree format
phylo_tree = Phylo.read(data_path / "SimBac/clonal_frame.nwk", "newick")

edge, node_height = newick_to_tree(phylo_tree)
tree.edge = edge
tree.node_height = node_height
tree.height = np.max(node_height)
tree.length = np.sum(edge[:, 2])

rho_site = 0.02
theta_site = 0.05
L = 200
delta = 30

x_o = ClonalOrigin_seq_sim(tree, rho_site, theta_site, L, delta)

def prior():
    theta = []
    theta.append(np.random.uniform(0.0, 0.2))    # rho_site
    theta.append(np.random.uniform(1.0, 100.0))  # delta
    theta.append(np.random.uniform(0.0, 0.2))    # theta_site
    theta.append(np.random.randint(100, 500))  # L
    theta = np.array(theta)
    return theta

def simulator(theta):
    summary_stats = ClonalOrigin_seq_sim(tree,
                                         theta[0].item(),
                                         theta[2].item(),
                                         int(theta[3].item()),
                                         theta[1].item())
    return summary_stats

def run_task(seed_value):
    """
    This function isolates the workload for a single CPU core.
    It sets the seed, generates the prior, and runs the simulation.
    """
    np.random.seed(seed_value)
    theta_val = prior()
    stats = simulator(theta_val)

    return seed_value, theta_val, stats

if __name__ == "__main__":
    Phylo.draw_ascii(phylo_tree)
    
    simulation_budget = 50
    
    # Matrices to store results
    theta_matrix = np.zeros((simulation_budget, 4))
    x_matrix = np.zeros((simulation_budget, 7))
    
    # If running on an HPC with SLURM, use the allocated CPUs
    if "SLURM_CPUS_PER_TASK" in os.environ:
        num_cores = int(os.environ["SLURM_CPUS_PER_TASK"])
    else:
        # Local laptop (Windows/Linux)
        num_cores = max(1, multiprocessing.cpu_count() - 1) 

    print(f"Running on {num_cores} cores...")

    # Execute in parallel
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        # map() automatically distributes the range of seeds to the available cores
        seeds = range(simulation_budget)
        results = executor.map(run_task, seeds)

    # Gather results as they finish and place them in the matrices
    for seed, theta_val, stats in results:
        theta_matrix[seed, :] = theta_val
        x_matrix[seed, :] = stats
        print(f"Completed Random seed: {seed}")
        
    print("All simulations finished.")

    # Reset directory and save simulation results
    output_file = data_path
    np.savetxt(output_file / 'x_o.csv', x_o, delimiter=",")
    np.savetxt(output_file / 'theta.csv', theta_matrix, delimiter=",")
    np.savetxt(output_file / 'x.csv', x_matrix, delimiter=",")
