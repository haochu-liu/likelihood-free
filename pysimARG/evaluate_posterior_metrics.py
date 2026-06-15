import numpy as np


def evaluate_posterior_metrics(true_value, posterior_matrix):
    """
    Computes evaluation metrics for posterior samples across multiple runs.
    
    Parameters:
    -----------
    true_value : float
        The true parameter value (phi). Cannot be exactly zero.
    posterior_matrix : numpy.ndarray
        A 2D array of shape (N_runs, N_samples) containing the sampled data.
        
    Returns:
    --------
    list
        A list containing coverage, bias, posterior_variance, relative_bias, and relative_rmse.
    """
    if true_value == 0:
        raise ValueError("true_value cannot be 0, as it is used in the denominator for relative metrics.")

    nrows, nsamples = posterior_matrix.shape
    posterior_means = np.mean(posterior_matrix, axis=1)
    posterior_variances = np.var(posterior_matrix, axis=1, ddof=1)
    lower_bounds = np.percentile(posterior_matrix, 2.5, axis=1)
    upper_bounds = np.percentile(posterior_matrix, 97.5, axis=1)

    is_covered = (true_value >= lower_bounds) & (true_value <= upper_bounds)
    coverage = np.mean(is_covered)

    bias = np.mean(np.abs(posterior_means - true_value))

    avg_posterior_variance = np.mean(posterior_variances)

    relative_bias = []
    relative_rmse = []
    for i in range(nrows):
        relative_bias.append(np.mean((posterior_matrix[i, :] - true_value) / true_value))
        relative_rmse.append(np.sqrt(np.mean(((posterior_matrix[i, :] - true_value) / true_value) ** 2)))

    return [coverage, bias, avg_posterior_variance, np.mean(relative_bias), np.mean(relative_rmse)]
