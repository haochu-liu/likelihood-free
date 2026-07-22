import numpy as np
from scipy import stats
import torch
from tqdm.auto import tqdm
import itertools
import torch
import torch.nn as nn
from torch.distributions import Uniform
import sbi
from sbi.utils.user_input_checks import MultipleIndependent
from sbi.neural_nets import posterior_nn
from sbi.neural_nets.embedding_nets import FCEmbedding
from sbi.inference import NPE_C
from sbi.diagnostics import run_sbc, check_sbc
import warnings
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


def SBC_KStest(ranks, num_posterior_samples, parameter_labels):
    num_dimensions = ranks.shape[1] 

    ks_results = []
    p_values = []
    for dim in range(num_dimensions):
        normalized_ranks = ranks[:, dim] / num_posterior_samples
        ks_stat, p_value = stats.kstest(normalized_ranks, 'uniform')
        ks_results.append(ks_stat)
        p_values.append(p_value)
    
    return ks_results, p_values


def mahalanobis_error(theta_est_post, theta_test_numpy):
    maha_errors = np.full((theta_est_post.shape[0]), np.nan)
    for i in range(theta_est_post.shape[0]):
        samples = theta_est_post[i]
        truth = theta_test_numpy[i]
        post_mean = np.mean(samples, axis=0)
        cov_matrix = np.cov(samples, rowvar=False)

        try:
            inv_cov = np.linalg.inv(cov_matrix)
        except np.linalg.LinAlgError:
            print(f"Warning: Singular covariance matrix at index {i}, returning NaN.")
            continue
        
        diff = post_mean - truth
        maha_dist_sq = np.dot(np.dot(diff, inv_cov), diff)
        maha_errors[i] = np.sqrt(maha_dist_sq)
    return maha_errors


if __name__ == "__main__":
    print("Loading data...")
    drop_col = range(16, 32)
    theta_test = np.loadtxt(str(data_path / 'ClonalOrigin/rho_and_theta/theta_sbc.csv'), delimiter=",")
    x_test = np.loadtxt(str(data_path / 'ClonalOrigin/rho_and_theta/x_sbc.csv'), delimiter=",")
    x_test = np.delete(x_test, drop_col, axis=1)
    print(theta_test.shape, x_test.shape)

    nan_row_test = np.where(np.isnan(x_test) | np.isinf(x_test))[0]
    print(nan_row_test)

    theta_test = np.delete(theta_test, nan_row_test, axis=0)
    theta_test = torch.tensor(theta_test, device=torch_device)
    theta_test = theta_test.to(torch.float32)
    theta_test_numpy = theta_test.cpu().numpy()

    x_test = np.delete(x_test, nan_row_test, axis=0)
    x_test = torch.tensor(x_test, device=torch_device)
    x_test = x_test.to(torch.float32)
    x_test_numpy = x_test.cpu().numpy()

    print(theta_test.shape, x_test.shape)

    theta1 = np.loadtxt(str(data_path / 'ClonalOrigin/rho_and_theta/theta1.csv'), delimiter=",")
    x1 = np.loadtxt(str(data_path / 'ClonalOrigin/rho_and_theta/x1.csv'), delimiter=",")
    theta2 = np.loadtxt(str(data_path / 'ClonalOrigin/rho_and_theta/theta2.csv'), delimiter=",")
    x2 = np.loadtxt(str(data_path / 'ClonalOrigin/rho_and_theta/x2.csv'), delimiter=",")

    x = np.vstack([x1, x2])
    x = np.delete(x, drop_col, axis=1)
    theta = np.vstack([theta1, theta2])
    print(theta.shape, x.shape)

    nan_row = np.where(np.isnan(x) | np.isinf(x))[0]
    print(nan_row)

    theta = np.delete(theta, nan_row, axis=0)
    theta = torch.tensor(theta[:10000, :], device=torch_device)
    theta = theta.to(torch.float32)
    theta_numpy = theta.cpu().numpy()

    x = np.delete(x, nan_row, axis=0)
    x = torch.tensor(x[:10000, :], device=torch_device)
    x = x.to(torch.float32)
    x_numpy = x.cpu().numpy()

    print(theta.shape, x.shape)

    print("Model tuning setup...")
    seeds = [1, 2, 3, 4, 5]
    num_posterior_samples = 1000
    stage5_array = np.array([[2, 1, 128],
                            [2, 4, 256],
                            [2, 4, 48]])
    print(stage5_array)

    stage5_indices = [0, 1, 2]
    flow_transforms = [3, 5, 8]
    hidden_features = [30, 50, 80, 120]
    spline_bins = [5, 10, 20]
    stage6_array = np.array(list(itertools.product(stage5_indices, flow_transforms, hidden_features, spline_bins)))
    print(stage6_array)

    stage6_p_values = np.full((stage6_array.shape[0], 3), np.nan)
    stage6_D_stats = np.full((stage6_array.shape[0], 3), np.nan)
    stage6_maha_errors = np.full((stage6_array.shape[0]), np.nan)
    stage6_nll = np.full((stage6_array.shape[0]), np.nan)

    print("Starting model tuning...")
    prior_rho = Uniform(low=torch.tensor([0.0]), high=torch.tensor([0.1]))
    prior_theta = Uniform(low=torch.tensor([0.0]), high=torch.tensor([0.1]))
    prior_L = DiscreteUniform(low=torch.tensor([100.0]), high=torch.tensor([10000.0]))
    prior = MultipleIndependent(
        dists=[prior_rho, prior_theta, prior_L],
        validate_args=False,
        device=torch_device
    )

    for k in range(stage6_array.shape[0]):
        # Stage 5 configuration
        stage5_config = stage5_array[stage6_array[k, 0]]
        num_outputs = stage5_config[0]
        num_hidden_layers = stage5_config[1]
        num_hiddens = stage5_config[2]

        # Stage 6 configuration
        num_transforms = stage6_array[k, 1]
        num_features = stage6_array[k, 2]
        num_bins = stage6_array[k, 3]

        embedding_net = LeaveLengthOut_NN(
            input_dim=30,
            num_hiddens=num_hiddens,
            num_hidden_layers=num_hidden_layers,
            num_outputs=num_outputs)
        neural_posterior = posterior_nn(
            model="nsf",
            embedding_net=embedding_net,
            hidden_features=num_features,
            num_transforms=num_transforms,
            num_bins=num_bins
        )
        print(f"Running iteration {k}")
        print(f"Stage 5 output dimension {num_outputs}, hidden layers {num_hidden_layers}, hidden units {num_hiddens}.")
        print(f"Stage 6 transforms {num_transforms}, features {num_features}, bins {num_bins}.")
        print("-" * 50)
        
        seed = seeds[0]
        torch.manual_seed(seed)
        np.random.seed(seed)

        inference_baseline = NPE_C(prior=prior, density_estimator=neural_posterior, device=torch_device)
        density_estimator_baseline = inference_baseline.append_simulations(theta, x).train(
            max_num_epochs=500
        )
        posterior_baseline = inference_baseline.build_posterior(density_estimator_baseline)

        theta_est_post = np.full((theta_test.shape[0], num_posterior_samples, 3), np.nan)
        for j in tqdm(range(theta_test.shape[0]), desc="Sampling posterior"):
            theta_post = posterior_baseline.sample((num_posterior_samples,), x=x_test[j, :],
                                                show_progress_bars=False, reject_outside_prior=False)
            theta_est_post[j, :, :] = theta_post.detach().numpy()

        parameter_labels = [r"for $\rho_s$", r"for $\theta_s$", r"for L"]
        theta_test_expanded = theta_test.unsqueeze(1)
        theta_est_post_tensor = torch.tensor(theta_est_post, device=torch_device)
        theta_est_post_tensor = theta_est_post_tensor.to(torch.float32)
        is_less_than_truth = theta_est_post_tensor < theta_test_expanded
        ranks = torch.sum(is_less_than_truth, dim=1)

        ks_results, p_values = SBC_KStest(ranks, num_posterior_samples, parameter_labels)
        stage6_p_values[k, :] = p_values
        stage6_D_stats[k, :] = ks_results

        stage6_maha_errors[k] = np.mean(mahalanobis_error(theta_est_post, theta_test_numpy))

        lp = density_estimator_baseline.log_prob(theta_test, x_test)
        stage6_nll[k] = -lp.detach().cpu().mean().item()

    print("Saving results...")
    np.save(str(data_path / 'NPE_tuning/stage6_p_values.npy'), stage6_p_values)
    np.save(str(data_path / 'NPE_tuning/stage6_D_stats.npy'), stage6_D_stats)
    np.save(str(data_path / 'NPE_tuning/stage6_maha_errors.npy'), stage6_maha_errors)
    np.save(str(data_path / 'NPE_tuning/stage6_nll.npy'), stage6_nll)
