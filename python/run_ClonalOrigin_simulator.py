import numpy as np
import torch
from sbi.utils.torchutils import BoxUniform
from sbi.inference import NPE_C, simulate_for_sbi
from sbi.utils.user_input_checks import (
    check_sbi_inputs,
    process_prior,
    process_simulator,
)
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from pysimARG.clonal_genealogy import ClonalTree
from pysimARG.ClonalOrigin_simulator import ClonalOrigin_simulator

torch_device = "cpu"


np.random.seed(100)
tree = ClonalTree(n=10)

rho_site = 0.02
theta_site = 0.05
L = 100000
delta = 300

x_o = ClonalOrigin_simulator(tree, rho_site, theta_site, L, delta, N=2000, k_vec=[50, 200, 2000])
x_o = torch.tensor(x_o, device=torch_device)
x_o = x_o.flatten()
x_o_numpy = x_o.cpu().numpy()

prior = BoxUniform(
    low=torch.tensor([0.0, 1.0, 0.0], device=torch_device),
    high=torch.tensor([0.2, 500.0, 0.2], device=torch_device), 
    device=torch_device
)

def simulator(theta):
    theta = theta.reshape(-1)
    summary_stats = ClonalOrigin_simulator(tree,
                                           theta[0].item(),
                                           theta[2].item(),
                                           L,
                                           theta[1].item(),
                                           N=100, k_vec=[50, 200, 2000])
    summary_stats = torch.tensor(summary_stats, device=torch_device)
    return summary_stats

prior, num_parameters, prior_returns_numpy = process_prior(prior)
simulator = process_simulator(simulator, prior, prior_returns_numpy)
check_sbi_inputs(simulator, prior)

simulation_budget = 50
seed = 100

inference = NPE_C(prior=prior, density_estimator="nsf", device=torch_device)
torch.manual_seed(seed)
np.random.seed(seed)

theta, x = simulate_for_sbi(
    simulator=simulator, proposal=prior, num_simulations=simulation_budget, num_workers=10
)

np.savetxt("data/ClonalOrigin/x_o.csv", x_o_numpy, delimiter=",", fmt='%d')
np.savetxt("data/ClonalOrigin/theta.csv", theta.cpu().numpy(), delimiter=",", fmt='%d')
np.savetxt("data/ClonalOrigin/x.csv", x.cpu().numpy(), delimiter=",", fmt='%d')
