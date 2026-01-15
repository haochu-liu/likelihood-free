import torch
import numpy as np
from config import torch_device
from sbi.utils.torchutils import BoxUniform


mean_radius = 0.1
sd_radius = 0.01
baseoffset = 0.25

prior_torch = BoxUniform(
    low=-1 * torch.ones(2), high=1 * torch.ones(2), device=torch_device
)

x_o = torch.tensor([0, 0], device=torch_device)
x_o = x_o.flatten()


def simulator_torch(theta):
    theta = theta.reshape(-1)
    a = torch.tensor(torch.pi * (torch.rand(1) - 0.5), device=torch_device)
    r = mean_radius + torch.randn(1, device=torch_device) * sd_radius
    p = torch.tensor(
        [r * torch.cos(a) + baseoffset, r * torch.sin(a)], device=torch_device
    )
    z0 = -torch.abs(theta[0] + theta[1]) / torch.sqrt(torch.tensor(2.0))
    z1 = (-theta[0] + theta[1]) / torch.sqrt(torch.tensor(2.0))
    return p + torch.tensor([z0, z1], device=torch_device)


def simulator_torch_batched(theta):
    n_sim, data_dim = theta.shape
    x_samples = torch.zeros((n_sim, data_dim))

    for i in range(n_sim):
        x_samples[i,] = simulator_torch(theta[i,])

    return x_samples


def simulator_numpy(theta):
    theta = theta.reshape(-1)
    a = np.array(np.pi * (np.random.random(1) - 0.5))
    r = mean_radius + np.random.normal(loc=0, scale=1, size=1) * sd_radius
    p = np.array([r * np.cos(a) + baseoffset, r * np.sin(a)])
    z0 = -np.abs(theta[0] + theta[1]) / np.sqrt(2.0)
    z1 = (-theta[0] + theta[1]) / np.sqrt(2.0)
    return p + np.array([z0, z1])


def analytic_posterior_numpy(x_o, n_samples=1):
    theta = np.zeros((n_samples, 2))
    for i in range(n_samples):
        p = simulator_numpy(np.zeros(2))
        q = np.zeros(2)
        q[0] = (p[0] - x_o[0]) * np.sqrt(2)
        q[1] = (x_o[1] - p[1]) * np.sqrt(2)

        if np.random.rand() < 0.5:
            q[0] = -q[0]

        theta[i, 0] = (q[0] - q[1]) / 2
        theta[i, 1] = (q[0] + q[1]) / 2
    return theta
