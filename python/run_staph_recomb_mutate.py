import numpy as np
from Bio import Phylo
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from tqdm import tqdm
import sys
from pathlib import Path
current_path = Path(__file__).resolve()
project_root = current_path.parent.parent
pysimARG_path = project_root / 'pysimARG'
data_path = project_root / 'data'

sys.path.append(str(pysimARG_path))
from clonal_genealogy import ClonalTree
from ClonalOrigin_seq_sim import ClonalOrigin_seq_sim


torch_device = "cpu"
output_file = data_path


def run_sim(params, clonal_tree):
    rho, theta, L = params
    return ClonalOrigin_seq_sim(clonal_tree, rho, theta, int(L), 453.6649321544094)


if __name__ == "__main__":
    np.random.seed(100)
    clonal_tree = ClonalTree(n=110)

    # Load phylo tree and convert to ClonalTree format
    phylo_tree = Phylo.read(data_path / "staph/saureus_clonal.nwk", "newick")
    Phylo.draw_ascii(phylo_tree)

    clonal_edge = np.loadtxt(data_path / "staph/clonal_edge.csv", delimiter=",", dtype=float)
    clonal_node_height = np.loadtxt(data_path / "staph/clonal_node_height.csv", delimiter=",", dtype=float)
    clonal_tree.edge = clonal_edge
    clonal_tree.node_height = clonal_node_height
    clonal_tree.height = np.max(clonal_node_height)
    clonal_tree.length = np.sum(clonal_edge[:, 2])

    seed = 1
    np.random.seed(seed)
    prior_rho = np.random.uniform(0, 0.1, size=100)
    prior_theta = np.random.uniform(0, 0.05, size=100)
    prior_L = np.random.randint(20, 10000, size=100)
    prior_param = np.column_stack((prior_rho, prior_theta, prior_L))

    np.savetxt(output_file / 'theta_test.csv', prior_param, delimiter=",")

    tasks = list(zip(prior_rho, prior_theta, prior_L))

    run_sim_with_tree = partial(run_sim, clonal_tree=clonal_tree)
    
    chunk_size = 10
    total_sims = len(tasks)
    
    with ProcessPoolExecutor(max_workers=8) as executor, tqdm(total=total_sims, desc="Running Simulations") as pbar:
        for chunk_start in range(0, total_sims, chunk_size):
            chunk_tasks = tasks[chunk_start : chunk_start + chunk_size]
            chunk_results = []
            for result in executor.map(run_sim_with_tree, chunk_tasks):
                chunk_results.append(result)
                pbar.update(1)

            # Using the memory-efficient "append" mode so we don't hold 10,000 results in RAM
            chunk_array = np.array(chunk_results)
            
            with open(output_file / 'x_test.csv', "a") as f:
                np.savetxt(f, chunk_array, delimiter=",")
