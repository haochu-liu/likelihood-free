import numpy as np
import pandas as pd
import torch
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import os
import warnings
from Bio import Phylo
from torch.distributions import Uniform
from sbi.utils.user_input_checks import MultipleIndependent
from sbi.inference import NPE_C, posterior_estimator_based_potential, MCMCPosterior
from sbi.analysis import conditional_potential
import sys
from pathlib import Path
current_path = Path(__file__).resolve()
project_root = current_path.parent.parent
pysimARG_path = project_root / 'pysimARG'
data_path = project_root / 'data'

sys.path.append(str(pysimARG_path))
from discrete_uniform import DiscreteUniform

torch_device = "cpu"

warnings.filterwarnings("ignore", category=UserWarning)


def process_mcmc(i):
    x_o = x_obs_torch[i:i+1]

    potential_fn, theta_transform = posterior_estimator_based_potential(
        density_estimator,
        prior=prior,
        x_o=x_o,
    )

    condition = torch.zeros(theta_dim)
    condition[fixed_dim] = fixed_value

    conditioned_potential_fn, restricted_tf, restricted_prior = conditional_potential(
        potential_fn=potential_fn,
        theta_transform=theta_transform,
        prior=prior,
        condition=condition,
        dims_to_sample=dims_to_sample,
    )

    mcmc_posterior = MCMCPosterior(
        potential_fn=conditioned_potential_fn,
        theta_transform=restricted_tf,
        proposal=restricted_prior,
        method="slice_np_vectorized",
        num_chains=10,
    )

    cond_samples = mcmc_posterior.sample(
        (num_posterior_samples,),
        show_progress_bars=False, 
    )

    return i, cond_samples.cpu().detach().numpy()


if __name__ == "__main__":
    # Plot the Bacillus clonal tree
    phylo_tree = Phylo.read(str(data_path / "bacillus" / "bacillus.nwk"), "newick")
    Phylo.draw_ascii(phylo_tree)

    # Load the observed data and simulation data
    x_obs_df = pd.read_csv(str(data_path / "bacillus" / "bacillus_block_summary_stats.csv"), header=None)
    x_obs_np = x_obs_df.to_numpy()
    x_obs_torch = torch.tensor(x_obs_np, device=torch_device)
    x_obs_torch = x_obs_torch.to(torch.float32)

    theta1 = np.loadtxt(str(data_path / "bacillus" / "ClonalOrigin_sim" / "change_theta" / "theta1.csv"), delimiter=",")
    x1 = np.loadtxt(str(data_path / "bacillus" / "ClonalOrigin_sim" / "change_theta" / "x1.csv"), delimiter=",")
    theta2 = np.loadtxt(str(data_path / "bacillus" / "ClonalOrigin_sim" / "change_theta" / "theta2.csv"), delimiter=",")
    x2 = np.loadtxt(str(data_path / "bacillus" / "ClonalOrigin_sim" / "change_theta" / "x2.csv"), delimiter=",")
    theta3 = np.loadtxt(str(data_path / "bacillus" / "ClonalOrigin_sim" / "change_theta" / "theta3.csv"), delimiter=",")
    x3 = np.loadtxt(str(data_path / "bacillus" / "ClonalOrigin_sim" / "change_theta" / "x3.csv"), delimiter=",")
    theta4 = np.loadtxt(str(data_path / "bacillus" / "ClonalOrigin_sim" / "change_theta" / "theta4.csv"), delimiter=",")
    x4 = np.loadtxt(str(data_path / "bacillus" / "ClonalOrigin_sim" / "change_theta" / "x4.csv"), delimiter=",")
    theta5 = np.loadtxt(str(data_path / "bacillus" / "ClonalOrigin_sim" / "change_theta" / "theta5.csv"), delimiter=",")
    x5 = np.loadtxt(str(data_path / "bacillus" / "ClonalOrigin_sim" / "change_theta" / "x5.csv"), delimiter=",")
    x = np.vstack([x1, x2, x3, x4, x5])
    theta = np.vstack([theta1, theta2, theta3, theta4, theta5])

    print("Begin training...")

    prior_rho = Uniform(low=torch.tensor([0.0]), high=torch.tensor([0.2]))
    prior_theta = Uniform(low=torch.tensor([0.0]), high=torch.tensor([0.1]))
    prior_L = DiscreteUniform(low=torch.tensor([100.0]), high=torch.tensor([10000.0]))
    prior = MultipleIndependent(
        dists=[prior_rho, prior_theta, prior_L],
        validate_args=False,
        device=torch_device
    )

    seed = 100
    num_posterior_samples=1000
    learning_rate = 0.0005

    inference = NPE_C(prior=prior, density_estimator="nsf", device=torch_device)
    torch.manual_seed(seed)
    np.random.seed(seed)

    density_estimator = inference.append_simulations(theta, x).train(
        max_num_epochs=500, learning_rate=learning_rate
    )
    posterior = inference.build_posterior(density_estimator)

    print("End training and begin sampling...")

    theta_dim = prior.event_shape[0]
    fixed_dim = 1                     # Dimensional index for theta
    fixed_value = 0.04755576699972153 # Fixed value for theta

    dims_to_sample = [d for d in range(theta_dim) if d != fixed_dim]

    theta_conditional = np.full((x_obs_torch.shape[0], num_posterior_samples, 2), np.nan)
    num_obs = x_obs_torch.shape[0]
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(process_mcmc, i): i for i in range(num_obs)}
        # Process results as they finish with a single, clean progress bar
        for future in tqdm(as_completed(futures), total=num_obs, desc="Running MCMC Chains"):
            i, result = future.result()
            theta_conditional[i, :, :] = result
    
    print("End sampling.")

np.save(str(data_path / "bacillus" / "conditional_posterior_samples.npy"), theta_conditional)
