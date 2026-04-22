import numpy as np
import torch
from torch.distributions import Uniform
from sbi.utils.user_input_checks import MultipleIndependent
from sbi.inference import simulate_for_sbi
from sbi.utils.user_input_checks import (
    check_sbi_inputs,
    process_prior,
    process_simulator,
)
import sys
from pathlib import Path
current_path = Path(__file__).resolve()
project_root = current_path.parent.parent
pysimARG_path = project_root / 'pysimARG'
data_path = project_root / 'data'

sys.path.append(str(pysimARG_path))
from clonal_genealogy import ClonalTree
from ClonalOrigin_seq_sim import ClonalOrigin_seq_sim
from discrete_uniform import DiscreteUniform

torch_device = "cpu"


np.random.seed(100)
tree = ClonalTree(n=15)

rho_site = 0.02
theta_site = 0.05
L = 2000
delta = 300

x_o = ClonalOrigin_seq_sim(tree, rho_site, theta_site, L, delta, k_vec=[20, 50, 90])
x_o = torch.tensor(x_o, device=torch_device)
x_o = x_o.flatten()
x_o_numpy = x_o.cpu().numpy()

prior_rho = Uniform(low=torch.tensor([0.0]), high=torch.tensor([0.2]))
prior_delta = Uniform(low=torch.tensor([1.0]), high=torch.tensor([500.0]))
prior_theta = Uniform(low=torch.tensor([0.0]), high=torch.tensor([0.2]))
prior_L = DiscreteUniform(low=torch.tensor([100.0]), high=torch.tensor([10000.0]))

prior = MultipleIndependent(
    dists=[prior_rho, prior_delta, prior_theta, prior_L],
    validate_args=False,
    device=torch_device
)

def simulator(theta):
    theta = theta.reshape(-1)
    summary_stats = ClonalOrigin_seq_sim(tree,
                                         theta[0].item(),
                                         theta[2].item(),
                                         theta[3].item(),
                                         theta[1].item(),
                                         k_vec=[20, 50, 90])
    summary_stats = torch.tensor(summary_stats, device=torch_device)
    return summary_stats


if __name__ == "__main__":
    # Run simulations
    prior, num_parameters, prior_returns_numpy = process_prior(prior)
    simulator = process_simulator(simulator, prior, prior_returns_numpy)
    check_sbi_inputs(simulator, prior)

    simulation_budget = 50
    seed = 100

    torch.manual_seed(seed)
    np.random.seed(seed)

    theta, x = simulate_for_sbi(
        simulator=simulator, proposal=prior, num_simulations=simulation_budget, num_workers=10
    )

    # Reset directory and save simulation results
    output_file = data_path
    np.savetxt(output_file / 'x_o.csv', x_o_numpy, delimiter=",")
    np.savetxt(output_file / 'theta.csv', theta.cpu().numpy(), delimiter=",")
    np.savetxt(output_file / 'x.csv', x.cpu().numpy(), delimiter=",")
