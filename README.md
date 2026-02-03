# Likelihood-Free Methods

A project-based implementation and comparison of likelihood-free inference methods for Bayesian inference, including Approximate Bayesian Computation (ABC), Bayesian Synthetic Likelihood (BSL), and Simulation-Based Inference (SBI).

## Overview

This repository implements and compares various likelihood-free inference methods for complex statistical models where the likelihood function is intractable or computationally expensive to evaluate. The methods are applied to multiple benchmark problems and real-world applications in population genetics.

## Methods Implemented

### 1. Approximate Bayesian Computation (ABC)

ABC methods approximate the posterior distribution by accepting parameters that generate simulated data "close" to the observed data, measured by a distance metric and tolerance threshold.

**R Implementations** ([R/ABC](R/ABC)):
- **ABC-MCMC** ([ABC_MCMC.R](R/ABC/ABC_MCMC.R)): Metropolis-Hastings MCMC with ABC acceptance criterion
- **ABC-AM** ([ABC_AM.R](R/ABC/ABC_AM.R)): Adaptive Metropolis algorithm with covariance adaptation during burn-in
- **ABC-SMC** ([ABC_SMC.R](R/ABC/ABC_SMC.R)): Sequential Monte Carlo with adaptive tolerance scheduling
- **Kernel Functions** ([gaussian_kernel.R](R/ABC/gaussian_kernel.R), [gaussian_kernel2.R](R/ABC/gaussian_kernel2.R)): Gaussian kernels for smooth ABC

**Python Implementations**:
- **ABC-MCMC** Using the `sbi` library's MCABC implementation

### 2. Bayesian Synthetic Likelihood (BSL)

BSL constructs a synthetic likelihood by approximating the intractable likelihood with a multivariate normal distribution fitted to summary statistics from simulated data.

**R Implementations** ([R/BSL](R/BSL)):
- **SL-MCMC** ([SL_MCMC.R](R/BSL/SL_MCMC.R), [SL_MCMC2.R](R/BSL/SL_MCMC2.R)): MCMC sampling from the synthetic likelihood posterior
- **SL-AM** ([SL_AM.R](R/BSL/SL_AM.R)): Adaptive Metropolis for BSL
- **SL-SMC** ([SL_SMC.R](R/BSL/SL_SMC.R)): Sequential Monte Carlo for BSL with tempering
- **SL-IBIS** ([SL_IBIS.R](R/BSL/SL_IBIS.R)): Iterated Batch Importance Sampling
- **Waste-Free Variants** ([SL_WF_SMC.R](R/BSL/SL_WF_SMC.R), [SL_WF_IBIS.R](R/BSL/SL_WF_IBIS.R)): Waste-free SMC and IBIS algorithms that reuse simulations
- **Resampling Utilities** ([SL_SMC_resample.R](R/BSL/SL_SMC_resample.R)): Systematic resampling
- **Parallel Implementations** ([SL_WF_SMC_par.R](R/BSL/SL_WF_SMC_par.R), [SL_WF_IBIS_par.R](R/BSL/SL_WF_IBIS_par.R)): Parallelized versions for computational efficiency
- **Effective sample size** ([ESS_weight.R](R/BSL/ESS_weight.R), [ESS_weight2.R](R/BSL/ESS_weight2.R), [CESS_weight.R](R/BSL/CESS_weight.R)): Effective and conditional ESS calculations

### 3. Simulation-Based Inference (SBI)

Neural network-based methods that learn the posterior distribution or likelihood function using neural density estimation.

**Python Implementations** ([python](python)):

#### Neural Posterior Estimation (NPE)
- **NPE-C** Neural Posterior Estimation using Normalizing Flows (NSF)
- **SNPE-C** Sequential NPE with multiple rounds of refinement

#### Neural Likelihood Estimation (NLE)
- **NLE-A** Neural Likelihood Estimation using Masked Autoregressive Flows (MAF)
- **SNLE** Sequential NLE with proposal refinement

All SBI methods use the [sbi](https://github.com/sbi-dev/sbi) library with PyTorch.

## Applications and Examples

### 1. Two Moons Problem ([python/two_moons](python/two_moons), [output/two_moons](output/two_moons))
A challenging 2D bimodal posterior benchmark problem for comparing different inference methods:
- Simulator: [two_moon_sim.py](python/two_moons/two_moon_sim.py)
- Notebook: [two_moon.ipynb](python/two_moons/two_moon.ipynb)
- Comparison of NPE-C, SNPE-C, NLE, SNLE, and ABC-MCMC with 5000 simulations across 10 random seeds

### 2. Poisson Toy Example ([R/Poisson_toy](R/Poisson_toy), [matlab](matlab))
Reproduce simple Poisson model with BSL-MCMC from Price et al. (2018):
- Implementation: [poisson_toy.R](R/Poisson_toy/poisson_toy.R)
- Parallel versions: [poisson_toy_par.R](R/Poisson_toy/poisson_toy_par.R), [poisson_toy_par2.R](R/Poisson_toy/poisson_toy_par2.R)
- Parameter estimation: [poisson_toy_phat.R](R/Poisson_toy/poisson_toy_phat.R)

### 3. Ricker Model ([R/Ricker](R/Ricker), [matlab](matlab))
Reproduce a discrete-time population dynamics model with BSL-MCMC from Price et al. (2018):
- Simulator: [simulate_ricker.R](R/Ricker/simulate_ricker.R), [simulate_ricker.m](matlab/simulate_ricker.m)
- Summary statistics: [ricker_summstats.R](R/Ricker/ricker_summstats.R), [ricker_summstats.m](matlab/ricker_summstats.m)
- Synthetic likelihood parameter estimation: [ricker_phat.R](R/Ricker/ricker_phat.R), [phat_ricker.m](matlab/phat_ricker.m)
- Reproducing Price et al. (2018): [ricker_price2018Bayesian.R](R/Ricker/ricker_price2018Bayesian.R)

### 4. Birth-Death Process ([python](python), [R/simARG](R/simARG))
Continuous-time Markov chain for population dynamics:
- Python simulator: [birth_death_sim.py](python/birth_death_sim.py)
- SBI notebook: [birth_death_sbi.ipynb](python/birth_death_sbi.ipynb)
- R implementation: [run_birth_death_sim.R](R/simARG/run_birth_death_sim.R)
- Hitting time estimation using NPE-C

### 5. ClonalOrigin Model ([python](python), [R/simARG](R/simARG))
Bacterial recombination inference from genomic data:
- Simulator: [ClonalOrigin_sim.R](R/simARG/ClonalOrigin_sim.R)
- SBI notebook: [ClonalOrigin_sbi.ipynb](python/ClonalOrigin_sbi.ipynb)
- Inference for recombination rate, delta, and mutation rate parameters

### 6. Additional Examples ([R/example](R/example), [R/qmd](R/qmd))
- **Lotka-Volterra model**: Predator-prey dynamics ([Lotka-Volterra.R](R/example/Lotka-Volterra.R), [Lotka-Volterra.qmd](R/qmd/Lotka-Volterra.qmd))
- **Normal, Cauchy, Skew-normal distributions**: Methodological tests ([R/qmd](R/qmd))
- **Negative binomial model**: [negative_binomial.qmd](R/qmd/negative_binomial.qmd)

## Evaluation Metrics

- **Maximum Mean Discrepancy (MMD)** ([MMD.R](R/MMD.R)): Kernel-based distance between distributions using RBF kernel
- **Wasserstein Distance** ([wasserstein1.R](R/wasserstein1.R)): 1-dimensional Wasserstein distance for comparing posteriors
- **Effective Sample Size (ESS)**: For assessing particle quality in SMC methods
- **Acceptance Rate**: For MCMC convergence diagnostics

## Repository Structure

```
likelihood-free/
├── R/ # R implementations
│ ├── ABC/ # ABC methods
│ ├── BSL/ # BSL methods
│ ├── Ricker/ # Ricker model
│ ├── Poisson_toy/ # Poisson examples
│ ├── simARG/ # Population genetics examples
│ ├── example/ # Additional examples
│ ├── qmd/ # Quarto markdown documents
│ ├── MMD.R # MMD metric
│ └── wasserstein1.R # Wasserstein distance
├── python/ # Python implementations
│ ├── two_moons/ # Two moons benchmark
│ ├── birth_death_sbi.ipynb # Birth-death SBI
│ ├── ClonalOrigin_sbi.ipynb# ClonalOrigin SBI
│ └── birth_death_sim.py # Birth-death simulator
├── matlab/ # Reproduce Price et al. (2018) examples
│ ├── bayes_sl_ricker_wood.m # BSL for Ricker
│ └── run_*.m # Run scripts
├── data/ # Results and datasets
├── output/ # Output files and results
└── sbi_env/ # Python virtual environment
```
