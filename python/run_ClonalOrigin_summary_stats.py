import sys
import concurrent.futures
import numpy as np
from pathlib import Path


current_path = Path(__file__).resolve()
project_root = current_path.parent.parent
pysimARG_path = project_root / 'pysimARG'
data_path = project_root / 'data'

sys.path.append(str(pysimARG_path))
from clonal_genealogy import ClonalTree
from ClonalOrigin_simulator import ClonalOrigin_simulator


def sim_rho(x_val, tree, theta, L, delta):
    return ClonalOrigin_simulator(tree, x_val, theta, L, delta, N=20)

def sim_delta(x_val, tree, rho, theta, L):
    return ClonalOrigin_simulator(tree, rho, theta, L, x_val, N=20)

def sim_theta(x_val, tree, rho, L, delta):
    return ClonalOrigin_simulator(tree, rho, x_val, L, delta, N=20)


if __name__ == '__main__':
    
    np.random.seed(100)
    tree = ClonalTree(n=15)

    rho_site = 0.02
    theta_site = 0.05
    L = 1000000
    delta = 300

    x_mat = np.zeros((200, 3))
    x_mat[:, 0] = np.random.uniform(0, 0.2, size=200)
    x_mat[:, 1] = np.random.uniform(1, 500, size=200)
    x_mat[:, 2] = np.random.uniform(0, 0.2, size=200)

    num_cores = 10
    print(f"Starting parallel execution on {num_cores} cores...")

    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
        print("Simulating rho_site_mat...")
        rho_futures = [executor.submit(sim_rho, x, tree, theta_site, L, delta) for x in x_mat[:, 0]]
        rho_site_mat = np.array([f.result() for f in rho_futures])

        print("Simulating delta_mat...")
        delta_futures = [executor.submit(sim_delta, x, tree, rho_site, theta_site, L) for x in x_mat[:, 1]]
        delta_mat = np.array([f.result() for f in delta_futures])

        print("Simulating theta_site_mat...")
        theta_futures = [executor.submit(sim_theta, x, tree, rho_site, L, delta) for x in x_mat[:, 2]]
        theta_site_mat = np.array([f.result() for f in theta_futures])

    print("Saving files to data directory...")
    data_path.mkdir(parents=True, exist_ok=True) 

    np.savetxt(data_path / 'x_mat_python.csv', x_mat, delimiter=",")
    np.savetxt(data_path / 'rho_site_python.csv', rho_site_mat, delimiter=",")
    np.savetxt(data_path / 'delta_python.csv', delta_mat, delimiter=",")
    np.savetxt(data_path / 'theta_site_python.csv', theta_site_mat, delimiter=",")
    
    print("Simulation and saves complete!")
