import torch
import numpy as np
from sbi.inference import NLE_A, simulate_for_sbi
from sbi.neural_nets import likelihood_nn
from sbi.utils.user_input_checks import (
    check_sbi_inputs,
    process_prior,
    process_simulator,
)
import two_moon_sim
from config import torch_device

torch.set_default_device(torch_device)


def two_moons_snle(
    num_simulations_per_round,
    num_rounds,
    seed,
    prior,
    x_obs,
    simulator,
    num_posterior_samples=1000,
    dir_prefix="",
):
    prior, num_parameters, prior_returns_numpy = process_prior(prior)
    simulator = process_simulator(simulator, prior, prior_returns_numpy)
    check_sbi_inputs(simulator, prior)

    density_estimator_fun = likelihood_nn(
        model="maf",
        hidden_features=50,
        z_score_x="independent",
        z_score_theta="independent",
    )

    inference = NLE_A(
        prior=prior,
        density_estimator=density_estimator_fun,
        device=torch_device,
    )

    learning_rate = 0.0005

    torch.manual_seed(seed)
    np.random.seed(seed)

    # posteriors = []
    proposal = prior

    for i in range(num_rounds):
        theta, x = simulate_for_sbi(
            simulator=simulator,
            proposal=proposal,
            num_simulations=num_simulations_per_round,
        )

        density_estimator = inference.append_simulations(theta, x).train(
            max_num_epochs=100, learning_rate=learning_rate
        )
        posterior = inference.build_posterior(density_estimator).set_default_x(x_obs)
        # posteriors.append(posterior)
        proposal = posterior

    # Posterior inference
    # for i in range(num_rounds):
    #     theta_trained = posteriors[i].sample((num_posterior_samples,), x=x_obs)
    #     theta_trained = theta_trained.reshape((num_posterior_samples, 2))

    #     simulation_budget = int(num_simulations_per_round*num_rounds)
    #     np.savetxt(
    #         f"output/{dir_prefix}snle_post_sims{i+1}_seed{seed}.csv",
    #         theta_trained.detach().cpu().numpy(),
    #         delimiter=",",
    #     )
    theta_trained = proposal.sample((num_posterior_samples,), x=x_obs)
    theta_trained = theta_trained.reshape((num_posterior_samples, 2))

    simulation_budget = int(num_simulations_per_round*num_rounds)
    np.savetxt(
        f"output/two_moons/{dir_prefix}snle_post_sims{simulation_budget}_seed{seed}.csv",
            theta_trained.detach().cpu().numpy(),
            delimiter=",",
        )


if __name__ == "__main__":
    import config

    for seed in config.seeds:
        two_moons_snle(
            num_simulations_per_round=int(5000/config.num_rounds),
            num_rounds=config.num_rounds,
            seed=seed,
            prior=two_moon_sim_torch.prior_torch,
            x_obs=two_moon_sim_torch.x_o,
            simulator=two_moon_sim_torch.simulator_torch,
            num_posterior_samples=config.num_posterior_samples,
            dir_prefix="",
        )
