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


