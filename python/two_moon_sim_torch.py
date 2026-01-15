import torch
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
