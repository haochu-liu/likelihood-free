import numpy as np
from scipy import stats
import torch
import torch.nn as nn
from torch.distributions import Uniform
from tqdm.auto import tqdm
import os
import warnings
from sbi.utils.user_input_checks import MultipleIndependent
from sbi.neural_nets import posterior_nn
from sbi.inference import NPE_C
import sys
from pathlib import Path
current_path = Path(__file__).resolve()
project_root = current_path.parent.parent.parent
pysimARG_path = project_root / 'pysimARG'
data_path = project_root / 'data'

sys.path.append(str(pysimARG_path))
from discrete_uniform import DiscreteUniform
from LeaveLengthOut_NN import LeaveLengthOut_NN

torch_device = "cpu"

warnings.filterwarnings("ignore", category=UserWarning)


if __name__ == "__main__":
    drop_col = range(16, 32)
    theta_test = np.loadtxt(str(data_path / 'ClonalOrigin' / 'rho_and_theta' / 'theta_sbc.csv'), delimiter=",")
    x_test = np.loadtxt(str(data_path / 'ClonalOrigin' / 'rho_and_theta' / 'x_sbc.csv'), delimiter=",")

    theta_test = torch.tensor(theta_test, device=torch_device)
    theta_test = theta_test.to(torch.float32)
    theta_test_numpy = theta_test.cpu().numpy()

    x_test = np.delete(x_test, drop_col, axis=1)
    x_test = torch.tensor(x_test, device=torch_device)
    x_test = x_test.to(torch.float32)
    x_test_numpy = x_test.cpu().numpy()

    print(theta_test.shape, x_test.shape)

    nan_row = np.where(np.isnan(x_test_numpy))[0]
    theta_test = theta_test[~np.isnan(x_test_numpy).any(axis=1)]
    x_test = x_test[~np.isnan(x_test_numpy).any(axis=1)]

    theta_test_numpy = theta_test.cpu().numpy()
    x_test_numpy = x_test.cpu().numpy()

    theta1 = np.loadtxt(str(data_path / 'ClonalOrigin' / 'rho_and_theta' / 'theta1.csv'), delimiter=",")
    x1 = np.loadtxt(str(data_path / 'ClonalOrigin' / 'rho_and_theta' / 'x1.csv'), delimiter=",")
    theta2 = np.loadtxt(str(data_path / 'ClonalOrigin' / 'rho_and_theta' / 'theta2.csv'), delimiter=",")
    x2 = np.loadtxt(str(data_path / 'ClonalOrigin' / 'rho_and_theta' / 'x2.csv'), delimiter=",")
    theta3 = np.loadtxt(str(data_path / 'ClonalOrigin' / 'rho_and_theta' / 'theta3.csv'), delimiter=",")
    x3 = np.loadtxt(str(data_path / 'ClonalOrigin' / 'rho_and_theta' / 'x3.csv'), delimiter=",")
    theta4 = np.loadtxt(str(data_path / 'ClonalOrigin' / 'rho_and_theta' / 'theta4.csv'), delimiter=",")
    x4 = np.loadtxt(str(data_path / 'ClonalOrigin' / 'rho_and_theta' / 'x4.csv'), delimiter=",")

    x = np.vstack([x1, x2, x3, x4])
    x = np.delete(x, drop_col, axis=1)
    theta = np.vstack([theta1, theta2, theta3, theta4])

    drop_indices = np.unique(np.concatenate([np.where(np.isnan(x))[0], np.where(np.isinf(x))[0]]))
    theta = np.delete(theta, drop_indices, axis=0)
    x = np.delete(x, drop_indices, axis=0)

    theta = torch.tensor(theta, device=torch_device)
    theta = theta.to(torch.float32)
    theta_numpy = theta.cpu().numpy()

    x = torch.tensor(x, device=torch_device)
    x = x.to(torch.float32)
    x_numpy = x.cpu().numpy()

    budgets = [49948]
    seeds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    num_posterior_samples=500
    learning_rate = 0.0005

    prior_rho = Uniform(low=torch.tensor([0.0]), high=torch.tensor([0.1]))
    prior_theta = Uniform(low=torch.tensor([0.0]), high=torch.tensor([0.1]))
    prior_L = DiscreteUniform(low=torch.tensor([100.0]), high=torch.tensor([10000.0]))
    prior = MultipleIndependent(
        dists=[prior_rho, prior_theta, prior_L],
        validate_args=False,
        device=torch_device
    )

    npe_sbc_D = np.full((3, len(budgets)), np.nan)
    npe_sbc_p_values = np.full((3, len(budgets)), np.nan)
    npe_multidim_coverage_results = np.full((2, len(budgets)), np.nan)
    npe_marginal_coverage_results = np.full((3, len(budgets), 2), np.nan)
    npe_mean_error_results = np.full((theta_test.shape[0], len(budgets), 4), np.nan)

    torch.set_num_threads(1)
    for index_seed in range(len(seeds)):
        seed = seeds[index_seed]
        print("-" * 50)
        print(f"\nRunning NPE with seed {seed}")
        
        for i in range(len(budgets)):
            torch.manual_seed(seed)
            np.random.seed(seed)
            print("-" * 50)
            n_sim = budgets[i]
            x_train = x[:n_sim]
            theta_train = theta[:n_sim]

            print(f"\nTraining NPE with {n_sim} simulations")
            
            embedding_net = LeaveLengthOut_NN(
                input_dim=30,
                num_hiddens=48,
                num_hidden_layers=2,
                num_outputs=4)
            neural_posterior = posterior_nn(
                model="nsf",
                embedding_net=embedding_net,
                z_score_theta="independent",
                z_score_x="independent"
            )
            inference_benchmark = NPE_C(prior=prior, density_estimator=neural_posterior, device=torch_device)
            density_estimator_benchmark = inference_benchmark.append_simulations(theta_train, x_train).train(
                max_num_epochs=500, learning_rate=learning_rate
            )
            posterior_benchmark = inference_benchmark.build_posterior(density_estimator_benchmark)

            print(f"\nSampling posterior for NPE, n={n_sim}")

            theta_est_post = np.full((theta_test.shape[0], num_posterior_samples, 3), np.nan)
            for j in tqdm(range(theta_test.shape[0]), desc="Sampling posterior"):
                theta_post = posterior_benchmark.sample((num_posterior_samples,), x=x_test[j, :],
                                                        show_progress_bars=False, reject_outside_prior=False)
                theta_est_post[j, :, :] = theta_post.detach().numpy()

            print(f"\nSave results for n={n_sim}")

            file_name = f"npe_post_stage1_seed{seed}_sim{n_sim}.npy"
            np.save(str(data_path / "benchmark_extension" / file_name), theta_est_post)

        print("-" * 50)
