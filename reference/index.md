# Package index

## Package Overview

- [`IRTsimrel-package`](https://joonho112.github.io/IRTsimrel/reference/IRTsimrel-package.md)
  [`IRTsimrel`](https://joonho112.github.io/IRTsimrel/reference/IRTsimrel-package.md)
  : IRTsimrel: Reliability-Targeted Simulation for Item Response Data

## Calibration Algorithms

Core functions for reliability-targeted calibration

- [`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
  [`print(`*`<eqc_result>`*`)`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
  [`summary(`*`<eqc_result>`*`)`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
  : Empirical Quadrature Calibration (Algorithm 1: EQC/SQC)
- [`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  [`spc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  [`print(`*`<sac_result>`*`)`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  [`summary(`*`<sac_result>`*`)`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  : Stochastic Approximation Calibration (Algorithm 2: SAC)

## Simulation Functions

Generate latent abilities and item parameters

- [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md)
  : Simulate Latent Ability Distribution for IRT Studies (G-family)
- [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md)
  [`print(`*`<item_params>`*`)`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md)
  [`summary(`*`<item_params>`*`)`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md)
  : Simulate Item Parameters for IRT Studies
- [`simulate_response_data()`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md)
  : Simulate Item Response Data from Calibration Results

## Reliability Utilities

Compute reliability metrics and screening tools

- [`compute_rho_bar()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)
  [`compute_rho_tilde()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)
  : Compute Marginal Reliability from Simulated Test Information
- [`compute_rho_both()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_both.md)
  : Compute Both Reliability Metrics in a Single Pass
- [`compute_apc_init()`](https://joonho112.github.io/IRTsimrel/reference/compute_apc_init.md)
  : Analytic Pre-Calibration (APC) Initialization
- [`compute_reliability_tam()`](https://joonho112.github.io/IRTsimrel/reference/compute_reliability_tam.md)
  : Compute WLE and EAP Reliability Using TAM
- [`check_feasibility()`](https://joonho112.github.io/IRTsimrel/reference/check_feasibility.md)
  [`print(`*`<feasibility_check>`*`)`](https://joonho112.github.io/IRTsimrel/reference/check_feasibility.md)
  : Check Feasibility of Target Reliability
- [`rho_curve()`](https://joonho112.github.io/IRTsimrel/reference/rho_curve.md)
  [`print(`*`<rho_curve>`*`)`](https://joonho112.github.io/IRTsimrel/reference/rho_curve.md)
  : Compute Reliability as a Function of Scaling Factor

## Comparison & Visualization

Compare distributions and calibration results

- [`compare_shapes()`](https://joonho112.github.io/IRTsimrel/reference/compare_shapes.md)
  : Compare Multiple Distribution Shapes
- [`compare_eqc_sac()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md)
  [`compare_eqc_spc()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md)
  : Compare EQC and SAC Calibration Results

## S3 Methods

Print, plot, summary, coef, and predict methods

- [`check_feasibility()`](https://joonho112.github.io/IRTsimrel/reference/check_feasibility.md)
  [`print(`*`<feasibility_check>`*`)`](https://joonho112.github.io/IRTsimrel/reference/check_feasibility.md)
  : Check Feasibility of Target Reliability
- [`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
  [`print(`*`<eqc_result>`*`)`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
  [`summary(`*`<eqc_result>`*`)`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
  : Empirical Quadrature Calibration (Algorithm 1: EQC/SQC)
- [`print(`*`<latent_G>`*`)`](https://joonho112.github.io/IRTsimrel/reference/print.latent_G.md)
  : Print Method for latent_G Objects
- [`print(`*`<summary.eqc_result>`*`)`](https://joonho112.github.io/IRTsimrel/reference/print.summary.eqc_result.md)
  : Print Method for summary.eqc_result Objects
- [`print(`*`<summary.item_params>`*`)`](https://joonho112.github.io/IRTsimrel/reference/print.summary.item_params.md)
  : Print Method for summary.item_params Objects
- [`print(`*`<summary.latent_G>`*`)`](https://joonho112.github.io/IRTsimrel/reference/print.summary.latent_G.md)
  : Print Method for summary.latent_G Objects
- [`print(`*`<summary.sac_result>`*`)`](https://joonho112.github.io/IRTsimrel/reference/print.summary.sac_result.md)
  : Print Method for summary.sac_result Objects
- [`rho_curve()`](https://joonho112.github.io/IRTsimrel/reference/rho_curve.md)
  [`print(`*`<rho_curve>`*`)`](https://joonho112.github.io/IRTsimrel/reference/rho_curve.md)
  : Compute Reliability as a Function of Scaling Factor
- [`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  [`spc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  [`print(`*`<sac_result>`*`)`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  [`summary(`*`<sac_result>`*`)`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  : Stochastic Approximation Calibration (Algorithm 2: SAC)
- [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md)
  [`print(`*`<item_params>`*`)`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md)
  [`summary(`*`<item_params>`*`)`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md)
  : Simulate Item Parameters for IRT Studies
- [`plot(`*`<item_params>`*`)`](https://joonho112.github.io/IRTsimrel/reference/plot.item_params.md)
  : Plot method for item_params objects
- [`plot(`*`<latent_G>`*`)`](https://joonho112.github.io/IRTsimrel/reference/plot.latent_G.md)
  : Plot Method for latent_G Objects
- [`plot(`*`<sac_result>`*`)`](https://joonho112.github.io/IRTsimrel/reference/plot.sac_result.md)
  : Plot SAC Convergence Trajectory
- [`summary(`*`<latent_G>`*`)`](https://joonho112.github.io/IRTsimrel/reference/summary.latent_G.md)
  : Summary Method for latent_G Objects
- [`coef(`*`<eqc_result>`*`)`](https://joonho112.github.io/IRTsimrel/reference/coef.eqc_result.md)
  : Extract Calibrated Item Parameters from EQC Results
- [`coef(`*`<sac_result>`*`)`](https://joonho112.github.io/IRTsimrel/reference/coef.sac_result.md)
  : Extract Calibrated Item Parameters from SAC Results
- [`predict(`*`<eqc_result>`*`)`](https://joonho112.github.io/IRTsimrel/reference/predict.eqc_result.md)
  : Predict Reliability at New Scaling Factor Values (EQC)
- [`predict(`*`<sac_result>`*`)`](https://joonho112.github.io/IRTsimrel/reference/predict.sac_result.md)
  : Predict Reliability at New Scaling Factor Values (SAC)
- [`as.data.frame(`*`<item_params>`*`)`](https://joonho112.github.io/IRTsimrel/reference/as.data.frame.item_params.md)
  : Extract item parameters as data frame
- [`as.numeric(`*`<latent_G>`*`)`](https://joonho112.github.io/IRTsimrel/reference/as.numeric.latent_G.md)
  : Extract theta values from latent_G object
