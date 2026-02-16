import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from pysimARG.birth_death_sim import birth_death_sim
import numpy as np


birth_death_matrix = np.zeros((100, 2))
for i in range(1, 101):
    np.random.seed(i)
    result = birth_death_sim(100, 10)
    birth_death_matrix[i-1, :] = result
    print(f'Complete {i} iterations')

np.savetxt("data/birth_death_matrix.csv", birth_death_matrix, delimiter=",", fmt='%d')
