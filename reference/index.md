# Package index

## Package Overview

- [`IRTsimrel-package`](https://joonho112.github.io/IRTsimrel/reference/IRTsimrel-package.md)
  [`IRTsimrel`](https://joonho112.github.io/IRTsimrel/reference/IRTsimrel-package.md)
  : IRTsimrel: Reliability-Targeted Simulation for Item Response Data

## Calibration Algorithms

Core functions for reliability-targeted calibration

- [`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
  : Empirical Quadrature Calibration (Algorithm 1: EQC/SQC)
- [`spc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/spc_calibrate.md)
  [`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/spc_calibrate.md)
  : Stochastic Approximation Calibration (Algorithm 2: SPC/SAC)

## Simulation Functions

Generate latent abilities and item parameters

- [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md)
  : Simulate Latent Ability Distribution for IRT Studies (G-family)
- [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md)
  : Simulate Item Parameters for IRT Studies
- [`simulate_response_data()`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md)
  : Simulate Item Response Data from EQC Results

## Reliability Utilities

Compute reliability metrics

- [`compute_rho_bar()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)
  [`compute_rho_tilde()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)
  : Compute Marginal Reliability from Simulated Test Information
- [`compute_apc_init()`](https://joonho112.github.io/IRTsimrel/reference/compute_apc_init.md)
  : Analytic Pre-Calibration (APC) Initialization
- [`compute_reliability_tam()`](https://joonho112.github.io/IRTsimrel/reference/compute_reliability_tam.md)
  : Compute WLE and EAP Reliability Using TAM

## Comparison & Visualization

Compare distributions and calibration results

- [`compare_shapes()`](https://joonho112.github.io/IRTsimrel/reference/compare_shapes.md)
  : Compare Multiple Distribution Shapes
- [`compare_eqc_spc()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_spc.md)
  : Compare EQC and SPC Calibration Results

## S3 Methods

Print, plot, and summary methods

- [`print(`*`<latent_G>`*`)`](https://joonho112.github.io/IRTsimrel/reference/print.latent_G.md)
  : Print Method for latent_G Objects
- [`plot(`*`<item_params>`*`)`](https://joonho112.github.io/IRTsimrel/reference/plot.item_params.md)
  : Plot method for item_params objects
- [`plot(`*`<latent_G>`*`)`](https://joonho112.github.io/IRTsimrel/reference/plot.latent_G.md)
  : Plot Method for latent_G Objects
- [`plot(`*`<spc_result>`*`)`](https://joonho112.github.io/IRTsimrel/reference/plot.spc_result.md)
  : Plot SPC Convergence Trajectory
- [`summary(`*`<latent_G>`*`)`](https://joonho112.github.io/IRTsimrel/reference/summary.latent_G.md)
  : Summary Method for latent_G Objects
- [`as.data.frame(`*`<item_params>`*`)`](https://joonho112.github.io/IRTsimrel/reference/as.data.frame.item_params.md)
  : Extract item parameters as data frame
- [`as.numeric(`*`<latent_G>`*`)`](https://joonho112.github.io/IRTsimrel/reference/as.numeric.latent_G.md)
  : Extract theta values from latent_G object
