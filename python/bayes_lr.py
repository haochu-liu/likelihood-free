import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # Hides info and warning messages
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0' # Hides the oneDNN specific message
import numpy as np
from pathlib import Path
import keras
import bayesflow as bf


np.set_printoptions(suppress=True)

def likelihood(beta, sigma, N):
    # x: predictor variable
    x = np.random.normal(0, 1, size=N)
    # y: response variable
    y = np.random.normal(beta[0] + beta[1] * x, sigma, size=N)
    return dict(y=y, x=x)

data_draws = likelihood(beta = [2, 1], sigma = 1, N = 3)
print(data_draws["y"].shape)
print(data_draws["y"])

def prior():
    # beta: regression coefficients (intercept, slope)
    beta = np.random.normal([2, 0], [3, 1])
    # sigma: residual standard deviation
    sigma = np.random.gamma(1, 1)
    return dict(beta=beta, sigma=sigma)

prior_draws = prior()
print(prior_draws["beta"].shape)
print(prior_draws["beta"])

def meta():
    # N: number of observation in a dataset
    N = np.random.randint(5, 15)
    return dict(N=N)

meta_draws = meta()
print(meta_draws["N"])

simulator = bf.simulators.make_simulator([prior, likelihood], meta_fn=meta)

