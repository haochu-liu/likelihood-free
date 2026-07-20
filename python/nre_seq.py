import numpy as np
import time
from scipy import stats
import torch
from tqdm.auto import tqdm
from joblib import Parallel, delayed
import os
import warnings
from torch.distributions import Uniform
from sbi.utils.user_input_checks import MultipleIndependent
from sbi.inference import NRE_B
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


def SBC_KStest(ranks, num_posterior_samples, parameter_labels):
    print("Kolmogorov-Smirnov Test Results (SBC):")
    print("-" * 40)

    num_dimensions = ranks.shape[1] 

    ks_results = []
    p_values = []
    for dim in range(num_dimensions):
        normalized_ranks = ranks[:, dim] / num_posterior_samples
        ks_stat, p_value = stats.kstest(normalized_ranks, 'uniform')
        ks_results.append(ks_stat)
        p_values.append(p_value)

        print(f"Parameter:        {parameter_labels[dim]}")
        print(f"KS Statistic (D): {ks_stat:.4f}")
        print(f"p-value:          {p_value:.4e}")

        if p_value < 0.05:
            print("Status: MISCALIBRATED (reject null)")
        else:
            print("Status: CALIBRATED (fail to reject null)")
        print("-" * 40)
    
    return ks_results, p_values


def multidim_coverage95(theta_est_post, theta_test_numpy):
    print("\nJoint 3D Coverage Metric:")
    print("-" * 40)

    num_test_samples = theta_test_numpy.shape[0]
    joint_covered_count = 0

    for i in range(num_test_samples):
        samples = theta_est_post[i]  # Shape: (1000, 3)
        truth = theta_test_numpy[i]  # Shape: (3,)

        mean_vec = np.mean(samples, axis=0)
        cov_matrix = np.cov(samples, rowvar=False)
        
        try:
            inv_cov = np.linalg.inv(cov_matrix)
        except np.linalg.LinAlgError:
            print(f"Warning: Singular covariance matrix at index {i}, skipping.")
            continue

        diff_samples = samples - mean_vec
        dist_samples = np.sum(np.dot(diff_samples, inv_cov) * diff_samples, axis=1)

        diff_truth = truth - mean_vec
        dist_truth = np.dot(np.dot(diff_truth, inv_cov), diff_truth)

        threshold_95 = np.percentile(dist_samples, 95)

        if dist_truth <= threshold_95:
            joint_covered_count += 1

    joint_coverage_95 = joint_covered_count / num_test_samples

    print(f"Target Joint Coverage: 95.0%")
    print(f"Actual Joint Coverage: {joint_coverage_95 * 100:.1f}%")
    print("-" * 40)

    return joint_coverage_95


def marginal_coverage95(theta_est_post, theta_test_numpy):
    print("\nMarginal 1D Coverage Metrics (95% Interval):")
    print("-" * 40)

    lower_bounds = np.percentile(theta_est_post, 2.5, axis=1)
    upper_bounds = np.percentile(theta_est_post, 97.5, axis=1)
    is_covered = (theta_test_numpy >= lower_bounds) & (theta_test_numpy <= upper_bounds)
    marginal_coverages = np.mean(is_covered, axis=0)

    for d, coverage in enumerate(marginal_coverages):
        print(f"Dimension {d} Coverage: {coverage * 100:.1f}%")
    print("-" * 40)
    
    return marginal_coverages


def L2_error(theta_est_post, theta_test_numpy):
    post_mean_error = np.full((theta_est_post.shape[0]), np.nan)
    for i in range(theta_est_post.shape[0]):
        post_mean = np.mean(theta_est_post[i, :, :], axis=0)
        error = post_mean - theta_test_numpy[i, :]
        post_mean_error[i] = np.linalg.norm(error)
    return post_mean_error


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

    budgets = [200, 500, 1000, 2000, 5000, 10000, 49948]
    seed = 1
    num_posterior_samples=1000
    learning_rate = 0.0005

    prior_rho = Uniform(low=torch.tensor([0.0]), high=torch.tensor([0.1]))
    prior_theta = Uniform(low=torch.tensor([0.0]), high=torch.tensor([0.1]))
    prior_L = DiscreteUniform(low=torch.tensor([100.0]), high=torch.tensor([10000.0]))
    prior = MultipleIndependent(
        dists=[prior_rho, prior_theta, prior_L],
        validate_args=False,
        device=torch_device
    )

    nre_sbc_D = np.full((3, len(budgets)), np.nan)
    nre_sbc_p_values = np.full((3, len(budgets)), np.nan)
    nre_multidim_coverage_results = np.full((2, len(budgets)), np.nan)
    nre_marginal_coverage_results = np.full((3, len(budgets), 2), np.nan)
    nre_mean_error_results = np.full((theta_test.shape[0], len(budgets), 4), np.nan)

    torch.set_num_threads(1)
    for i in range(len(budgets)):
        torch.manual_seed(seed)
        np.random.seed(seed)
        print("-" * 50)
        n_sim = budgets[i]
        x_train = x[:n_sim]
        theta_train = theta[:n_sim]

        print(f"\nTraining NRE with {n_sim} simulations")
        
        inference_benchmark = NRE_B(prior=prior, classifier="resnet", device=torch_device)
        density_estimator_benchmark = inference_benchmark.append_simulations(theta_train, x_train).train(
            max_num_epochs=500, learning_rate=learning_rate
        )
        posterior_benchmark = inference_benchmark.build_posterior(density_estimator_benchmark)

        print(f"Sampling posterior for NRE, n={n_sim}")

        theta_est_post = posterior_benchmark.sample_batched((num_posterior_samples,), x=x_test,
                                                            num_chains=10,
                                                            show_progress_bars=True)
        theta_est_post = theta_est_post.permute(1, 0, 2)
        theta_est_post = theta_est_post.numpy()
        print(theta_est_post.shape)
        
        print(f"Calculating metrics for NRE, n={n_sim}")

        sbc_results = []
        parameter_labels = [r"for $\rho_s$", r"for $\theta_s$", r"for L"]
        theta_test_expanded = theta_test.unsqueeze(1)
        theta_est_post_tensor = torch.tensor(theta_est_post, device=torch_device)
        theta_est_post_tensor = theta_est_post_tensor.to(torch.float32)
        is_less_than_truth = theta_est_post_tensor < theta_test_expanded
        ranks = torch.sum(is_less_than_truth, dim=1)

        ks_results, p_values = SBC_KStest(ranks, num_posterior_samples, parameter_labels)
        nre_sbc_D[:, i] = ks_results
        nre_sbc_p_values[:, i] = p_values

        nre_multidim_coverage_results[0, i] = multidim_coverage95(theta_est_post, theta_test_numpy)
        nre_multidim_coverage_results[1, i] = multidim_coverage95(theta_est_post[:, :, :2], theta_test_numpy[:, :2])
        nre_marginal_coverage_results[:, i, 0] = marginal_coverage95(theta_est_post, theta_test_numpy)
        nre_marginal_coverage_results[:, i, 1] = marginal_coverage95(theta_est_post[:, :, :2], theta_test_numpy[:, :2])

        nre_mean_error_results[:, i, 0] = L2_error(theta_est_post, theta_test_numpy)
        nre_mean_error_results[:, i, 1] = L2_error(theta_est_post[:, :, :2], theta_test_numpy[:, :2])
        nre_mean_error_results[:, i, 2] = mahalanobis_error(theta_est_post, theta_test_numpy)
        nre_mean_error_results[:, i, 3] = mahalanobis_error(theta_est_post[:, :, :2], theta_test_numpy[:, :2])

        print(f"Save results for n={n_sim}")

        np.savetxt(str(data_path / "benchmark" / "nre_sbc_D1.csv"), nre_sbc_D, delimiter=",")
        np.savetxt(str(data_path / "benchmark" / "nre_sbc_p_values1.csv"), nre_sbc_p_values, delimiter=",")
        np.save(str(data_path / "benchmark" / "nre_multidim_coverage_results1.npy"), nre_multidim_coverage_results)
        np.save(str(data_path / "benchmark" / "nre_marginal_coverage_results1.npy"), nre_marginal_coverage_results)
        np.save(str(data_path / "benchmark" / "nre_mean_error_results1.npy"), nre_mean_error_results)

        print("-" * 50)
