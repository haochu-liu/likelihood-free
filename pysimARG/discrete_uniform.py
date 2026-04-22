import torch
from torch.distributions import Distribution, constraints


class DiscreteUniform(Distribution):
    arg_constraints = {'low': constraints.real, 'high': constraints.real}

    def __init__(self, low, high, validate_args=None):
        self.low = torch.as_tensor(low, dtype=torch.float32)
        self.high = torch.as_tensor(high, dtype=torch.float32)
        self._num_categories = (self.high - self.low) + 1
        
        batch_shape = self.low.size()
        super().__init__(batch_shape=batch_shape, validate_args=validate_args)

    def sample(self, sample_shape=torch.Size()):
        low_int = int(self.low.min().item())
        high_int = int(self.high.max().item()) + 1
        expanded_shape = sample_shape + self.batch_shape
        samples = torch.randint(low_int, high_int, expanded_shape)
        return samples.float()

    def log_prob(self, value):
        log_p = -torch.log(self._num_categories)
        in_bounds = (value >= self.low) & (value <= self.high)
        return torch.where(in_bounds, log_p, torch.tensor(float('-inf')))
    
    @property
    def support(self):
        return constraints.interval(self.low, self.high)

    @property
    def mean(self):
        return (self.high + self.low) / 2.0

    @property
    def variance(self):
        n = self._num_categories
        return (n**2 - 1) / 12.0
    
    @property
    def stddev(self):
        return torch.sqrt(self.variance)
