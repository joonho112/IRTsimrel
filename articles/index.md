# Articles

### Getting Started

- [IRTsimrel: Reliability-Targeted IRT
  Simulation](https://joonho112.github.io/IRTsimrel/articles/introduction.md):

  Gateway vignette introducing the IRTsimrel package for
  reliability-targeted simulation of item response data. Covers the core
  challenge, how the package works, the two calibration algorithms (EQC
  and SAC), reliability metrics, exported functions, and a reading
  roadmap for both applied and methodological researchers.

- [Quick Start: Your First Calibration in 5
  Minutes](https://joonho112.github.io/IRTsimrel/articles/quick-start.md):

  Learn to calibrate item parameters for target reliability, generate
  response data, and extract results in just 5 minutes.

### Applied Researchers Track

- [Applied Guide: Reliability-Targeted IRT
  Simulation](https://joonho112.github.io/IRTsimrel/articles/applied-guide.md):

  A comprehensive practical workflow guide for reliability-targeted IRT
  simulation using IRTsimrel, covering the full 6-step pipeline from
  feasibility screening through response data generation.

- [Working with Latent
  Distributions](https://joonho112.github.io/IRTsimrel/articles/latent-distributions.md):

  A comprehensive guide to generating latent ability distributions using
  sim_latentG(), covering all 12 built-in shapes, custom mixtures,
  covariate effects, and the connection to reliability-targeted
  simulation.

- [Generating Realistic Item
  Parameters](https://joonho112.github.io/IRTsimrel/articles/item-parameters.md):

  A guide to generating realistic item parameters with
  sim_item_params(), covering parametric, IRW, hierarchical, and custom
  sources, the copula method for correlated parameters, and extracting
  calibrated items with coef().

- [Designing IRT Simulation Studies with Reliability
  Control](https://joonho112.github.io/IRTsimrel/articles/simulation-design.md):

  A comprehensive guide to designing factorial IRT simulation studies
  with reliability as a controlled design factor using IRTsimrel.

- [Case Studies: Publication-Ready Simulation
  Templates](https://joonho112.github.io/IRTsimrel/articles/case-studies.md):

  Three self-contained case studies demonstrating reliability-targeted
  IRT simulation for model comparison, robustness analysis, and sample
  size planning.

### Methodological Researchers Track

- [Mathematical Foundations: Reliability Theory and the Inverse Design
  Problem](https://joonho112.github.io/IRTsimrel/articles/theory-reliability.md):

  A rigorous treatment of marginal reliability theory under the 2PL IRT
  model, the global discrimination scaling framework, Jensen’s
  inequality and its consequences, monotonicity and existence-uniqueness
  results, and the theoretical foundations of EQC and SAC calibration
  algorithms.

- [Algorithm 1: Empirical Quadrature Calibration
  (EQC)](https://joonho112.github.io/IRTsimrel/articles/algorithm-eqc.md):

  Detailed treatment of the EQC algorithm for reliability-targeted IRT
  simulation, including formal algorithm statement, Brent’s root-finding
  method, convergence theory, Monte Carlo error analysis, numerical
  stability considerations, and diagnostic examples.

- [Algorithm 2: Stochastic Approximation Calibration
  (SAC)](https://joonho112.github.io/IRTsimrel/articles/algorithm-sac.md):

  Detailed treatment of the SAC algorithm for reliability-targeted IRT
  simulation, including the Robbins-Monro framework, step-size analysis,
  Polyak-Ruppert averaging, initialization strategies, convergence
  theory, and diagnostic examples.

- [Validating Calibration
  Results](https://joonho112.github.io/IRTsimrel/articles/validation.md):

  A comprehensive guide to validating reliability-targeted IRT
  calibration results at three levels: internal prediction checks,
  cross-algorithm comparison (EQC vs SAC), and external validation with
  TAM. Covers Jensen’s inequality diagnostics, the 960-condition
  validation study, misspecification analysis, and troubleshooting.

### Reference

- [Complete API Reference with
  Examples](https://joonho112.github.io/IRTsimrel/articles/api-reference.md):

  Comprehensive reference for all exported functions in IRTsimrel,
  organized by category with full signatures, parameter descriptions,
  return values, and working examples.
