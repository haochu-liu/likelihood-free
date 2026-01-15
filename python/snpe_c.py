import sbi.inference
import torch
import numpy as np
from sbi.inference import NPE_C, simulate_for_sbi
from sbi.utils.user_input_checks import (
    check_sbi_inputs,
    process_prior,
    process_simulator,
)
import two_moon_sim_torch
from config import torch_device

torch.set_default_device(torch_device)


def two_moons_snpe_c(
    num_simulations_per_round,
    num_rounds,
    seed,
    prior,
    x_obs,
    simulator,
    num_posterior_samples=1000,
    dir_prefix=""
):
    prior, num_parameters, prior_returns_numpy = process_prior(prior)
    simulator = process_simulator(simulator, prior, prior_returns_numpy)
    check_sbi_inputs(simulator, prior)

    inference = NPE_C(prior=prior, density_estimator="nsf", device=torch_device)

    learning_rate = 0.0005  # default value

    torch.manual_seed(seed)
    np.random.seed(seed)

    # posteriors = []
    proposal = prior

    for i in range(num_rounds):
        theta, x = simulate_for_sbi(
            simulator=simulator,
            proposal=proposal,
            num_simulations=num_simulations_per_round
        )

        density_estimator = inference.append_simulations(theta, x, proposal=proposal).train(
            max_num_epochs=100, learning_rate=learning_rate
        )
        posterior = inference.build_posterior(density_estimator).set_default_x(x_obs)
        # posteriors.append(posterior)
        proposal = posterior

    # for i in range(num_rounds):
    #     theta_trained = posteriors[i].sample((num_posterior_samples,), x=x_obs)
    #     theta_trained = theta_trained.reshape((num_posterior_samples, 2))

    #     np.savetxt(
    #         f"output/two_moons/{dir_prefix}snpec_post_round{i+1}_seed{seed}.csv",
    #         theta_trained.detach().cpu().numpy(),
    #         delimiter=",",
    #     )
    theta_trained = proposal.sample((num_posterior_samples,), x=x_obs)
    theta_trained = theta_trained.reshape((num_posterior_samples, 2))

    simulation_budget = int(num_simulations_per_round*num_rounds)
    np.savetxt(
        f"output/two_moons/{dir_prefix}snpec_post_sims{simulation_budget}_seed{seed}.csv",
            theta_trained.detach().cpu().numpy(),
            delimiter=",",
        )


if __name__ == "__main__":
    import config

    for seed in config.seeds:
        two_moons_snpe_c(
            num_simulations_per_round=int(5000/config.num_rounds),
            num_rounds=config.num_rounds,
            seed=seed,
            prior=two_moon_sim_torch.prior_torch,
            x_obs=two_moon_sim_torch.x_o,
            simulator=two_moon_sim_torch.simulator_torch,
            num_posterior_samples=config.num_posterior_samples,
            dir_prefix="",
        )
