import torch


class DiscretePosteriorWrapper:
    def __init__(self, continuous_posterior, discrete_indices):
        """
        Wraps an sbi posterior to automatically round specific dimensions.
        discrete_indices: A list of column indices that should be integers.
        """
        self.posterior = continuous_posterior
        self.discrete_indices = discrete_indices

    def sample(self, sample_shape=torch.Size(), **kwargs):
        samples = self.posterior.sample(sample_shape, **kwargs)
        for idx in self.discrete_indices:
            samples[..., idx] = torch.round(samples[..., idx])
            
        return samples

    def log_prob(self, theta, **kwargs):
        return self.posterior.log_prob(theta, **kwargs)
        
    def set_default_x(self, x):
        self.posterior.set_default_x(x)
