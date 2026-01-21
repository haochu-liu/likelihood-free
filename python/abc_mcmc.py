import torch
import numpy as np
from sbi.inference import MCABC
from sbi.utils.user_input_checks import (
    check_sbi_inputs,
    process_prior,
    process_simulator,
)
import two_moon_sim
from config import torch_device

torch.set_default_device(torch_device)

def two_moons_abc_mcmc(
    simulation_budget, seed, prior, x_obs, simulator, eps=1, dir_prefix=""
):
    prior, num_parameters, prior_returns_numpy = process_prior(prior)
    simulator = process_simulator(simulator, prior, prior_returns_numpy)
    check_sbi_inputs(simulator, prior)

    inference = MCABC(simulator, prior)
    samples = inference(x_obs, num_simulations=simulation_budget, eps=eps)

    np.savetxt(
        f"output/two_moons/{dir_prefix}abc_mcmc_post_sims{simulation_budget}_seed{seed}.csv",
        samples.cpu().numpy(),
        delimiter=",",
    )


if __name__ == "__main__":
    import config

    for seed in config.seeds:
        for simulation_budget in [5000]:
            two_moons_abc_mcmc(
                simulation_budget=simulation_budget,
                seed=seed,
                prior=two_moon_sim.prior_torch,
                x_obs=two_moon_sim.x_o,
                simulator=two_moon_sim.simulator_torch,
                eps=config.abc_mcmc_eps,
                dir_prefix="",
            )
