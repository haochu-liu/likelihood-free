import sys
import concurrent.futures
import numpy as np
from Bio import Phylo
from pathlib import Path


current_path = Path(__file__).resolve()
project_root = current_path.parent.parent
pysimARG_path = project_root / 'pysimARG'
data_path = project_root / 'data'

sys.path.append(str(pysimARG_path))
from clonal_genealogy import ClonalTree
from ClonalOrigin_pair_sim import ClonalOrigin_pair_sim
from newick_to_tree import newick_to_tree


def sim_rho(x_val, tree, theta, L, delta):
    return ClonalOrigin_pair_sim(tree, x_val, theta, L, delta)

def sim_delta(x_val, tree, rho, theta, L):
    return ClonalOrigin_pair_sim(tree, rho, theta, L, x_val)

def sim_theta(x_val, tree, rho, L, delta):
    return ClonalOrigin_pair_sim(tree, rho, x_val, L, delta)

def sim_L(x_val, tree, rho, theta, delta):
    return ClonalOrigin_pair_sim(tree, rho, theta, int(x_val), delta)


if __name__ == '__main__':
    
    np.random.seed(100)
    clonal_tree = ClonalTree(n=10)

    # Load phylo tree and convert to ClonalTree format
    phylo_tree = Phylo.read(data_path / "SimBac/clonal_frame.nwk", "newick")
    Phylo.draw_ascii(phylo_tree)

    edge, node_height = newick_to_tree(phylo_tree)
    clonal_tree.edge = edge
    clonal_tree.node_height = node_height
    clonal_tree.height = np.max(node_height)
    clonal_tree.length = np.sum(edge[:, 2])

    rho_site = 0.02
    theta_site = 0.05
    L = 2000
    delta = 300

    x_mat = np.zeros((20, 4))
    x_mat[:, 0] = np.random.uniform(0, 0.2, size=20)
    x_mat[:, 1] = np.random.uniform(1, 100, size=20)
    x_mat[:, 2] = np.random.uniform(0, 0.2, size=20)
    x_mat[:, 3] = np.random.randint(100, 501, size=20)

    num_cores = 8
    print(f"Starting parallel execution on {num_cores} cores...")

    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
        print("Simulating rho_site_mat...")
        rho_futures = [executor.submit(sim_rho, x, clonal_tree, theta_site, L, delta) for x in x_mat[:, 0]]
        rho_site_mat = np.array([f.result() for f in rho_futures])

        print("Simulating delta_mat...")
        delta_futures = [executor.submit(sim_delta, x, clonal_tree, rho_site, theta_site, L) for x in x_mat[:, 1]]
        delta_mat = np.array([f.result() for f in delta_futures])

        print("Simulating theta_site_mat...")
        theta_futures = [executor.submit(sim_theta, x, clonal_tree, rho_site, L, delta) for x in x_mat[:, 2]]
        theta_site_mat = np.array([f.result() for f in theta_futures])

        print("Simulating L_mat...")
        L_futures = [executor.submit(sim_L, x, clonal_tree, rho_site, theta_site, delta) for x in x_mat[:, 3]]
        L_mat = np.array([f.result() for f in L_futures])

    print("Saving files to data directory...")
    data_path.mkdir(parents=True, exist_ok=True) 

    np.savetxt(data_path / 'x_pair.csv', x_mat, delimiter=",")
    np.savetxt(data_path / 'rho_site_pair.csv', rho_site_mat, delimiter=",")
    np.savetxt(data_path / 'delta_pair.csv', delta_mat, delimiter=",")
    np.savetxt(data_path / 'theta_site_pair.csv', theta_site_mat, delimiter=",")
    np.savetxt(data_path / 'L_pair.csv', L_mat, delimiter=",")
    
    print("Simulation and saves complete!")
