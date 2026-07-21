import torch
import torch.nn as nn


class LeaveLengthOut_NN(nn.Module):
    """
    Configurable embedding network for sbi.

    The network maps the first input_dim features through an MLP and then appends
    the last element of x to the embedding output. Therefore, the final output
    dimension is num_outputs + 1.

    Parameters
    ----------
    input_dim : int
        Number of input features used by the MLP. If x has shape [batch, 30]
        and you also append x[:, -1:], then input_dim should usually be 30
        only if the last feature is also intentionally used inside the MLP.
    num_hiddens : int
        Width of each hidden layer when hidden_dims is not provided.
    num_hidden_layers : int
        Number of hidden layers when hidden_dims is not provided.
    num_outputs : int
        Output dimension of the MLP before appending the last input element.
    activation : type[nn.Module]
        PyTorch activation class, e.g. nn.ReLU, nn.ELU, nn.GELU, nn.SiLU.
    hidden_dims : list[int] | None
        Optional explicit hidden-layer sizes, e.g. [64, 128, 64]. If provided,
        this overrides num_hiddens and num_hidden_layers.
    """

    def __init__(
        self,
        input_dim=30,
        num_hiddens=48,
        num_hidden_layers=2,
        num_outputs=4,
        activation=nn.ReLU,
        hidden_dims=None
    ):
        super().__init__()

        if hidden_dims is None:
            hidden_dims = [num_hiddens] * num_hidden_layers

        layers = []
        previous_dim = input_dim

        for hidden_dim in hidden_dims:
            layers.append(nn.Linear(previous_dim, hidden_dim))
            layers.append(activation())
            previous_dim = hidden_dim
        layers.append(nn.Linear(previous_dim, num_outputs))

        self.mlp = nn.Sequential(*layers)

    def forward(self, x):
        mlp_out = self.mlp(x)
        last_element = x[:, -1:]
        out = torch.cat([mlp_out, last_element], dim=-1)
        return out
