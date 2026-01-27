import birth_death_sim
import numpy as np
import pandas as pd
import time


time_vec = []
for i in range(1, 101):
    np.random.seed(i)
    start_time = time.time()
    result = birth_death_sim.birth_death_sim(100, 10)
    end_time = time.time()
    time_vec.append(end_time - start_time)
    print(result)
    print(f'Complete {i} iterations')

# print(time_vec)
python_bd_time = pd.DataFrame({'python': time_vec})
python_bd_time.to_csv('data/python_bd_time.csv', index=False)
