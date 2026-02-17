import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from pysimARG.birth_death_sim import birth_death_sim
import numpy as np
import torch
from sbi.utils.torchutils import BoxUniform
from sbi.inference import NPE_C, simulate_for_sbi
from sbi.utils.user_input_checks import (
    check_sbi_inputs,
    process_prior,
    process_simulator,
)

torch_device = "cpu"


x_o = birth_death_sim(20, 3.0) # True value is 3.0
x_o = torch.tensor(x_o, device=torch_device)
x_o = x_o.flatten()
x_o_numpy = x_o.cpu().numpy()

prior = BoxUniform(
    low=torch.tensor([0.0], device=torch_device),
    high=torch.tensor([10.0], device=torch_device), 
    device=torch_device
)

def simulator(theta):
    theta = theta.reshape(-1)
    x_o = birth_death_sim(20, theta[0].item())
    x_o = torch.tensor(x_o, device=torch_device)
    return x_o

prior, num_parameters, prior_returns_numpy = process_prior(prior)
simulator = process_simulator(simulator, prior, prior_returns_numpy)
check_sbi_inputs(simulator, prior)

simulation_budget = 10000
seed = 100
num_posterior_samples=1000
learning_rate = 0.0005

inference = NPE_C(prior=prior, density_estimator="nsf", device=torch_device)
torch.manual_seed(seed)
np.random.seed(seed)

theta, x = simulate_for_sbi(
    simulator=simulator, proposal=prior, num_simulations=simulation_budget, num_workers=10
)

density_estimator = inference.append_simulations(theta, x).train(
    max_num_epochs=100, learning_rate=learning_rate
)
posterior = inference.build_posterior(density_estimator).set_default_x(x_o)

theta_trained = posterior.sample((num_posterior_samples,), x=x_o)
theta_trained = theta_trained.reshape(-1)

np.savetxt("data/theta_trained.csv", theta_trained.cpu().numpy(), delimiter=",", fmt='%d')

# plt.hist(theta_trained.cpu().numpy(), bins=30, color='skyblue', edgecolor='black')
# plt.xlabel('Birth Rate')
# plt.ylabel('Frequency')
# plt.title('Histogram of posterior')
# plt.show()
