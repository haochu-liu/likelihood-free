import numpy as np


def birth_death_sim(n, rho):
    """
    Hitting time for a birth-death process

    Simulate a birth-death process and find its hitting time at boundary equals to 1.

    param:
    n: A single integer for initial value.
    rho: The birth rate.
    return:
    Hitting time and number of jumps.
    """
    if not isinstance(n, (int, np.integer)):
        raise ValueError(f"`n` must be a single integer!")

    t = 0.0
    k = n
    num = 0
    
    while k > 1:
        rate = k * (k - 1 + rho) / 2
        t += np.random.exponential(scale=1/rate)

        if np.random.random() < (k - 1) / (k - 1 + rho):
            k -= 1
            num += 1
        else:
            k += 1
            num += 2
            
    return {"t": t, "n": num}
