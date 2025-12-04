# =============================================================================
# validation_utils.R
# =============================================================================
# Validation and Comparison Utilities for EQC/SPC Calibration
#
# Contents:
#   - compute_reliability_tam(): Compute WLE/EAP reliability using TAM
#   - simulate_response_data(): Generate item response data from EQC results
#   - compare_eqc_spc(): Compare EQC and SPC calibration results
#
# Author: JoonHo Lee
# Date: December 2025
# =============================================================================


#' Compute WLE and EAP Reliability Using TAM
#'
#' @description
#' Fits a Rasch or 2PL model using \pkg{TAM} and computes WLE and EAP
#' reliability using the official \code{WLErel()} and \code{EAPrel()} functions.
#'
#' @param resp Matrix or data.frame of item responses (0/1).
#' @param model Character. \code{"rasch"} or \code{"2pl"}.
#' @param verbose Logical. If TRUE, print fitting messages.
#' @param ... Additional arguments passed to TAM fitting functions.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{rel_wle}}{WLE reliability.}
#'   \item{\code{rel_eap}}{EAP reliability.}
#'   \item{\code{mod}}{Fitted TAM model object.}
#'   \item{\code{wle}}{Output from \code{TAM::tam.wle()}.}
#' }
#'
#' @details
#' ## WLE vs EAP Reliability
#'
#' TAM defines these reliability coefficients differently:
#'
#' - **WLE reliability**: \eqn{1 - \bar{s}^2 / V_{WLE}}, based on design effect
#' - **EAP reliability**: \eqn{V_{EAP} / (V_{EAP} + \bar{\sigma}^2)}, based on posterior variance
#'
#' Mathematically, \eqn{\rho_{EAP} \geq \rho_{WLE}} always holds under TAM's definitions.
#' EAP reliability more closely corresponds to MSEM-based population reliability.
#' For conservative inference, treat WLE as a lower bound.
#'
#' @examples
#' \dontrun{
#' # Simulate response data from EQC results
#' sim_data <- simulate_response_data(eqc_result, n_persons = 500)
#'
#' # Compute TAM reliability
#' tam_rel <- compute_reliability_tam(sim_data$response_matrix, model = "rasch")
#' cat(sprintf("WLE reliability: %.4f\n", tam_rel$rel_wle))
#' cat(sprintf("EAP reliability: %.4f\n", tam_rel$rel_eap))
#' }
#'
#' @seealso
#' \code{\link{simulate_response_data}} for generating test data,
#' \code{\link{eqc_calibrate}} for calibration.
#'
#' @export
compute_reliability_tam <- function(resp,
                                    model = c("rasch", "2pl"),
                                    verbose = FALSE,
                                    ...) {
  if (!requireNamespace("TAM", quietly = TRUE)) {
    stop("Package 'TAM' is required. Please install it first.")
  }

  model <- match.arg(model)
  resp  <- as.matrix(resp)

  # Fit IRT model
  if (model == "rasch") {
    mod <- TAM::tam.mml(resp = resp, verbose = verbose, ...)
  } else {
    mod <- TAM::tam.mml.2pl(resp = resp, irtmodel = "2PL", verbose = verbose, ...)
  }

  # --- WLE reliability ---
  wle <- TAM::tam.wle(mod)
  if (!all(c("theta", "error") %in% colnames(wle))) {
    stop("`TAM::tam.wle()` output must contain 'theta' and 'error'.")
  }
  rel_wle <- TAM::WLErel(theta = wle$theta, error = wle$error)

  # --- EAP reliability ---
  person <- mod$person
  if (!all(c("EAP", "SD.EAP") %in% colnames(person))) {
    stop("TAM model object does not contain 'EAP' and 'SD.EAP' in $person.")
  }
  rel_eap <- TAM::EAPrel(theta = person$EAP, error = person$SD.EAP)

  list(
    rel_wle = rel_wle,
    rel_eap = rel_eap,
    mod     = mod,
    wle     = wle
  )
}


#' Simulate Item Response Data from EQC Results
#'
#' @description
#' Generates item response data using the calibrated parameters from EQC.
#'
#' @param eqc_result An \code{eqc_result} object from \code{eqc_calibrate()}.
#' @param n_persons Integer. Number of persons to simulate.
#' @param latent_shape Character. Shape argument for \code{sim_latentG()}.
#' @param latent_params List. Additional arguments for \code{sim_latentG()}.
#' @param seed Optional integer for reproducibility.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{response_matrix}}{N x I matrix of binary responses}
#'   \item{\code{theta}}{True abilities (N x 1)}
#'   \item{\code{beta}}{Item difficulties (I x 1)}
#'   \item{\code{lambda}}{Item discriminations (I x 1)}
#' }
#'
#' @examples
#' \dontrun{
#' # First, run EQC calibration
#' eqc_result <- eqc_calibrate(
#'   target_rho = 0.80,
#'   n_items = 25,
#'   model = "rasch",
#'   seed = 42
#' )
#'
#' # Generate response data
#' sim_data <- simulate_response_data(
#'   eqc_result = eqc_result,
#'   n_persons = 1000,
#'   latent_shape = "normal",
#'   seed = 123
#' )
#'
#' # Use with TAM for validation
#' tam_rel <- compute_reliability_tam(sim_data$response_matrix, model = "rasch")
#' }
#'
#' @seealso
#' \code{\link{eqc_calibrate}}, \code{\link{compute_reliability_tam}}
#'
#' @export
simulate_response_data <- function(eqc_result,
                                   n_persons,
                                   latent_shape = "normal",
                                   latent_params = list(),
                                   seed = NULL) {

  if (!inherits(eqc_result, "eqc_result")) {
    stop("eqc_result must be an object of class 'eqc_result'")
  }

  if (!is.null(seed)) set.seed(seed)

  # Generate theta
  latent_args <- modifyList(
    list(n = n_persons, shape = latent_shape),
    latent_params
  )
  theta <- do.call(sim_latentG, latent_args)$theta

  beta   <- eqc_result$beta_vec
  lambda <- eqc_result$lambda_scaled
  I <- length(beta)

  # Generate responses (vectorized)
  lin_pred <- outer(theta, beta, "-")
  lin_pred <- sweep(lin_pred, 2, lambda, "*")
  prob_mat <- plogis(lin_pred)

  response_matrix <- matrix(
    rbinom(n_persons * I, size = 1, prob = as.vector(prob_mat)),
    nrow = n_persons, ncol = I
  )

  colnames(response_matrix) <- paste0("item", 1:I)

  list(
    response_matrix = response_matrix,
    theta = theta,
    beta = beta,
    lambda = lambda
  )
}


# =============================================================================
# Comparison Function for EQC and SPC Results
# =============================================================================

#' Compare EQC and SPC Calibration Results
#'
#' @description
#' Compares the calibration results from EQC and SPC algorithms.
#'
#' @param eqc_result An object of class \code{"eqc_result"}.
#' @param spc_result An object of class \code{"spc_result"}.
#' @param verbose Logical. If TRUE, print comparison summary.
#'
#' @return A list with comparison statistics (invisibly):
#' \describe{
#'   \item{\code{c_eqc}}{Calibrated c* from EQC.}
#'   \item{\code{c_spc}}{Calibrated c* from SPC.}
#'   \item{\code{diff_abs}}{Absolute difference between c* values.}
#'   \item{\code{diff_pct}}{Percent difference relative to EQC.}
#'   \item{\code{agreement}}{Logical. TRUE if difference is < 5%.}
#'   \item{\code{target_rho}}{Target reliability.}
#' }
#'
#' @examples
#' \dontrun{
#' # Run both algorithms
#' eqc_result <- eqc_calibrate(target_rho = 0.80, n_items = 25, seed = 42)
#' spc_result <- spc_calibrate(target_rho = 0.80, n_items = 25,
#'                             c_init = eqc_result, seed = 42)
#'
#' # Compare results
#' compare_eqc_spc(eqc_result, spc_result)
#' }
#'
#' @seealso
#' \code{\link{eqc_calibrate}}, \code{\link{spc_calibrate}}
#'
#' @export
compare_eqc_spc <- function(eqc_result, spc_result, verbose = TRUE) {

  if (!inherits(eqc_result, "eqc_result")) {
    stop("eqc_result must be an object of class 'eqc_result'.")
  }
  if (!inherits(spc_result, "spc_result")) {
    stop("spc_result must be an object of class 'spc_result'.")
  }

  # Extract key values
  c_eqc <- eqc_result$c_star
  c_spc <- spc_result$c_star

  diff_abs <- abs(c_eqc - c_spc)
  diff_pct <- 100 * diff_abs / c_eqc

  agreement <- diff_pct < 5  # Within 5%

  comparison <- list(
    c_eqc = c_eqc,
    c_spc = c_spc,
    diff_abs = diff_abs,
    diff_pct = diff_pct,
    agreement = agreement,
    target_rho = spc_result$target_rho
  )

  if (verbose) {
    cat("\n")
    cat("=======================================================\n")
    cat("  EQC vs SPC Comparison\n")
    cat("=======================================================\n\n")
    cat(sprintf("  Target reliability  : %.4f\n", comparison$target_rho))
    cat(sprintf("  EQC c*              : %.6f\n", c_eqc))
    cat(sprintf("  SPC c*              : %.6f\n", c_spc))
    cat(sprintf("  Absolute difference : %.6f\n", diff_abs))
    cat(sprintf("  Percent difference  : %.2f%%\n", diff_pct))
    cat(sprintf("  Agreement (< 5%%)    : %s\n", ifelse(agreement, "YES", "NO")))
    cat("\n")
  }

  invisible(comparison)
}
