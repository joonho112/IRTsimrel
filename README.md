
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IRTsimrel <img src="man/figures/logo.png" align="right" height="139" alt="IRTsimrel hex sticker showing a reliability curve with a target crosshair" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/joonho112/IRTsimrel/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/joonho112/IRTsimrel/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![CRAN
status](https://www.r-pkg.org/badges/version/IRTsimrel)](https://CRAN.R-project.org/package=IRTsimrel)
<!-- badges: end -->

**IRTsimrel** is an R package for **reliability-targeted simulation** of
Item Response Theory (IRT) data. Rather than treating reliability as an
afterthought, IRTsimrel lets you specify a **target marginal
reliability** as an explicit input—and calibrates the data-generating
process to achieve it exactly.

> *Marginal reliability in IRT serves the same role as ICC in multilevel
> modeling.*

## Why This Matters

In Monte Carlo IRT studies, researchers routinely vary sample size, test
length, and item parameters—but **marginal reliability is almost never
directly controlled**. This leads to three problems:

| Problem | Consequence |
|----|----|
| **Ecological validity** | Real assessments often have $\rho = 0.5$–$0.7$, yet simulations may implicitly generate $\rho > 0.90$ |
| **Confounded comparisons** | Conclusions about estimator superiority may only hold within a narrow reliability regime |
| **Limited replicability** | Without reporting implied reliability, exact replication is impossible |

Just as multilevel modelers always specify the ICC, IRT simulators
should always specify the target reliability. IRTsimrel makes this easy.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("joonho112/IRTsimrel")
```

## Quick Start

Calibrate a 25-item Rasch simulation to achieve $\rho = 0.80$:

``` r
library(IRTsimrel)

result <- eqc_calibrate(
  target_rho  = 0.80,
  n_items     = 25,
  model       = "rasch",
  item_source = "parametric",
  M           = 10000,
  seed        = 42
)

result
#> 
#> =======================================================
#>   Empirical Quadrature Calibration (EQC) Results
#> =======================================================
#> 
#> Calibration Summary:
#>   Model                        : RASCH
#>   Target reliability (rho*)    : 0.8000
#>   Achieved reliability         : 0.8000
#>   Absolute error               : 1.42e-07
#>   Scaling factor (c*)          : 0.9059
#> 
#> Design Parameters:
#>   Number of items (I)          : 25
#>   Quadrature points (M)        : 10000
#>   Reliability metric           : Average-information (tilde)
#>   Latent variance              : 1.0123
#> 
#> Convergence:
#>   Root status                  : uniroot_success
#>   Search bracket               : [0.300, 3.000]
#>   Bracket reliabilities        : [0.3540, 0.9532]
#> 
#> Parameter Summaries:
#>   theta:        mean = -0.011, sd = 1.006
#>   beta:         mean = -0.000, sd = 0.909, range = [-1.60, 1.73]
#>   lambda_base:  mean = 1.000, sd = 0.000
#>   lambda_scaled: mean = 0.906, sd = 0.000
```

Generate a response matrix and you’re ready for analysis:

``` r
sim_data <- simulate_response_data(result, n_persons = 1000, seed = 123)
dim(sim_data$response_matrix)
#> [1] 1000   25
```

See `vignette("quick-start")` for a 5-minute walkthrough, or
`vignette("applied-guide")` for the full 6-step workflow.

## How It Works

IRTsimrel separates **structure** from **scale**:

- **Structure** — Realistic item difficulties, discriminations, and
  latent distributions are drawn from empirically-grounded generators.
- **Scale** — A single global scaling factor $c^*$ is calibrated so that
  $\lambda_i^* = c^* \cdot \lambda_{i,0}$ yields the target reliability.

This means you can vary reliability independently of all other design
features—just like varying ICC in a multilevel simulation.

### Two Calibration Algorithms

| Algorithm | Function | Method | Recommended For |
|----|----|----|----|
| **EQC** | `eqc_calibrate()` | Deterministic root-finding (Brent’s method) | Default: fast, exact, no tuning |
| **SAC** | `sac_calibrate()` | Stochastic approximation (Robbins–Monro) | Cross-validation, complex DGPs |

When both algorithms converge to the same $c^*$, calibration is highly
trustworthy.

### 12 Pre-Standardized Latent Distributions

``` r
compare_shapes(
  n = 3000,
  shapes = c("normal", "bimodal", "skew_pos", "heavy_tail"),
  seed = 42
)
```

<img src="man/figures/README-latent-shapes-1.png" width="100%" />

All 12 built-in shapes (normal, bimodal, trimodal, skewed, heavy-tailed,
uniform, floor/ceiling effects) are **pre-standardized** to mean 0 and
variance 1. This ensures that distributional shape is cleanly separated
from scale—enabling rigorous factorial manipulation.

### Empirically Realistic Item Parameters

``` r
items <- sim_item_params(
  n_items = 30, model = "2pl",
  source = "irw", method = "copula",
  discrimination_params = list(rho = -0.3)
)
```

Generate item parameters from **four sources** (parametric, IRW,
hierarchical, custom) with empirically-observed
difficulty–discrimination correlation ($\rho \approx -0.3$) preserved
via copula methods.

### Feasibility Screening

Before calibrating, verify that your target reliability is achievable
and visualize the reliability–scaling curve:

``` r
check_feasibility(n_items = 25, model = "rasch", seed = 42)
rho_curve(n_items = 25, model = "rasch", seed = 42, plot = TRUE)
```

### Three-Level Validation

| Level | Tool | What It Checks |
|----|----|----|
| **Internal** | `predict()`, `compute_rho_both()` | Consistency of calibrated reliability; Jensen’s gap between $\tilde{\rho}$ and $\bar{w}$ |
| **Cross-algorithm** | `compare_eqc_sac()` | Agreement between EQC and SAC on $c^*$ |
| **External** | `compute_reliability_tam()` | TAM-fitted WLE/EAP reliability on generated response data |

## Documentation

IRTsimrel includes [12
vignettes](https://joonho112.github.io/IRTsimrel/articles/) organized in
two tracks:

**Getting Started** — `vignette("introduction")` \|
`vignette("quick-start")`

**Applied Researchers Track** — For running simulation studies:

- `vignette("applied-guide")` — 6-step workflow from feasibility check
  to data export
- `vignette("latent-distributions")` — All 12 shapes with examples
- `vignette("item-parameters")` — Parametric, IRW, hierarchical, and
  custom item generation
- `vignette("simulation-design")` — Designing factorial Monte Carlo
  studies with reliability as a factor
- `vignette("case-studies")` — Three complete publication-ready
  simulation templates

**Methodological Researchers Track** — For understanding the theory:

- `vignette("theory-reliability")` — Formal reliability definitions,
  Jensen’s inequality, and the inverse design problem
- `vignette("algorithm-eqc")` — EQC convergence theory and diagnostics
- `vignette("algorithm-sac")` — SAC step-size analysis and
  Polyak–Ruppert averaging
- `vignette("validation")` — Three-level validation framework
- `vignette("api-reference")` — Complete function reference with
  examples

## Citation

    Lee, J.-H. (2025). Reliability-targeted simulation for item response data.
    arXiv preprint arXiv:2512.16012. https://arxiv.org/abs/2512.16012

## Related Work

- [TAM](https://CRAN.R-project.org/package=TAM) — Test Analysis Modules
  for IRT
- [mirt](https://CRAN.R-project.org/package=mirt) — Multidimensional
  Item Response Theory
- [irw](https://github.com/ben-domingue/irw) — Item Response Warehouse

## License

MIT
