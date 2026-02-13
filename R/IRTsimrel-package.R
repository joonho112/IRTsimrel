#' IRTsimrel: Reliability-Targeted Simulation for Item Response Data
#'
#' @description
#' A framework for reliability-targeted simulation studies in Item Response
#' Theory (IRT). Researchers specify a target reliability level, and the
#' package calibrates item discriminations to achieve that target.
#'
#' @details
#' Main functions:
#' - \code{\link{eqc_calibrate}}: Empirical Quadrature Calibration (recommended)
#' - \code{\link{sac_calibrate}}: Stochastic Approximation Calibration
#' - \code{\link{sim_latentG}}: Generate latent ability distributions
#' - \code{\link{sim_item_params}}: Generate item parameters
#' - \code{\link{simulate_response_data}}: Generate response data
#'
#' @author JoonHo Lee
#' @docType package
#' @name IRTsimrel-package
#' @aliases IRTsimrel
#'
#' @importFrom stats var sd plogis rbinom rnorm rgamma rt runif uniroot optimize
#' @importFrom stats quantile median cor qnorm dnorm pnorm predict
#' @importFrom stats coef lm
#' @importFrom utils modifyList data head
#' @importFrom graphics hist par abline curve legend lines plot
#' @importFrom grDevices rgb
"_PACKAGE"

utils::globalVariables(c(
  "theta", "density", "c", "rho", "iteration", "phase",
  "beta", "lambda", "shape", "after_stat",
  "form_id", "item_id", "lambda_unscaled",
  "metric"
))
