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


def process_mcmc(i, x_obs_torch, density_estimator, prior, theta_dim, fixed_dim, fixed_value,
                 dims_to_sample, num_posterior_samples):
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
    # Plot the staph clonal tree
    phylo_tree = Phylo.read(str(data_path / "staph" / "saureus_clonal.nwk"), "newick")
    Phylo.draw_ascii(phylo_tree)

    # Load the observed data and simulation data
    x_obs_df = pd.read_csv(str(data_path / "staph" / "core_gene_summary_stats.csv"), index_col=0)
    x_obs_np = x_obs_df.to_numpy()

    no_signal_id = np.where(x_obs_np[:, 33] == 0)[0]
    x_obs_np = np.delete(x_obs_np, no_signal_id, axis=0)
    x_obs_torch = torch.tensor(x_obs_np, device=torch_device)
    x_obs_torch = x_obs_torch.to(torch.float32)

    theta1 = np.loadtxt(str(data_path / "staph" / "ClonalOrigin_sim" / "change_theta" / "theta1.csv"), delimiter=",")
    x1 = np.loadtxt(str(data_path / "staph" / "ClonalOrigin_sim" / "change_theta" / "x1.csv"), delimiter=",")
    theta2 = np.loadtxt(str(data_path / "staph" / "ClonalOrigin_sim" / "change_theta" / "theta2.csv"), delimiter=",")
    x2 = np.loadtxt(str(data_path / "staph" / "ClonalOrigin_sim" / "change_theta" / "x2.csv"), delimiter=",")
    theta3 = np.loadtxt(str(data_path / "staph" / "ClonalOrigin_sim" / "change_theta" / "theta3.csv"), delimiter=",")
    x3 = np.loadtxt(str(data_path / "staph" / "ClonalOrigin_sim" / "change_theta" / "x3.csv"), delimiter=",")
    theta4 = np.loadtxt(str(data_path / "staph" / "ClonalOrigin_sim" / "change_theta" / "theta4.csv"), delimiter=",")
    x4 = np.loadtxt(str(data_path / "staph" / "ClonalOrigin_sim" / "change_theta" / "x4.csv"), delimiter=",")
    theta5 = np.loadtxt(str(data_path / "staph" / "ClonalOrigin_sim" / "change_theta" / "theta5.csv"), delimiter=",")
    x5 = np.loadtxt(str(data_path / "staph" / "ClonalOrigin_sim" / "change_theta" / "x5.csv"), delimiter=",")
    x = np.vstack([x1, x2, x3, x4, x5])
    theta = np.vstack([theta1, theta2, theta3, theta4, theta5])

    print("Load training...")

    checkpoint = torch.load(
        str(data_path / "staph" / "trained_npe_density_estimator.pt"),
        map_location=torch_device,
        weights_only=False,  # may be needed in newer PyTorch versions
    )

    density_estimator = checkpoint["density_estimator"]
    prior = checkpoint["prior"]

    density_estimator = density_estimator.to(torch_device)

    inference = NPE_C(
        prior=prior,
        density_estimator="nsf",
        device=torch_device,
    )

    posterior = inference.build_posterior(density_estimator)

    num_posterior_samples = 1000

    print("Begin sampling...")

    theta_dim = prior.event_shape[0]
    fixed_dim = 1                      # Dimensional index for theta
    fixed_value = 0.001823019701987505 # Fixed value for theta

    dims_to_sample = [d for d in range(theta_dim) if d != fixed_dim]

    theta_conditional = np.full((x_obs_torch.shape[0], num_posterior_samples, 2), np.nan)
    num_obs = x_obs_torch.shape[0]
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = {}
        for i in range(num_obs):
            future = executor.submit(
                process_mcmc,
                i,
                x_obs_torch,
                density_estimator,
                prior,
                theta_dim,
                fixed_dim,
                fixed_value,
                dims_to_sample,
                num_posterior_samples
            )
            futures[future] = i
            
        # Catch the results as they finish, updating the progress bar
        for future in tqdm(as_completed(futures), total=num_obs, desc="MCMC Progress"):
            i, result = future.result()
            theta_conditional[i, :, :] = result
    
    print("End sampling.")

    np.save(str(data_path / "staph" / "conditional_posterior_samples.npy"), theta_conditional)
