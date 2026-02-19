import numpy as np
import sys
from pathlib import Path
current_path = Path(__file__).resolve()
project_root = current_path.parent.parent
pysimARG_path = project_root / 'pysimARG'
data_path = project_root / 'data'

sys.path.append(str(pysimARG_path))
from clonal_genealogy import ClonalTree
from ClonalOrigin_simulator import ClonalOrigin_simulator


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

rho_site_mat = np.zeros((200, 7))
for i in range(200):
    rho_site_mat[i, :] = ClonalOrigin_simulator(tree, x_mat[i, 0], theta_site, L, delta, N=2000)

delta_mat = np.zeros((200, 7))
for i in range(200):
    delta_mat[i, :] = ClonalOrigin_simulator(tree, rho_site, theta_site, L, x_mat[i, 1], N=2000)

theta_site_mat = np.zeros((200, 7))
for i in range(200):
    theta_site_mat[i, :] = ClonalOrigin_simulator(tree, rho_site, x_mat[i, 2], L, delta, N=2000)


np.savetxt(data_path / 'x_mat_python.csv', x_mat, delimiter=",")
np.savetxt(data_path / 'rho_site_python.csv', rho_site_mat, delimiter=",")
np.savetxt(data_path / 'delta_python.csv', delta_mat, delimiter=",")
np.savetxt(data_path / 'theta_site_python.csv', theta_site_mat, delimiter=",")
