# =============================================================================
# validation_utils.R
# =============================================================================
# Validation and Comparison Utilities for EQC/SAC Calibration
#
# Contents:
#   - compute_reliability_tam(): Compute WLE/EAP reliability using TAM
#   - simulate_response_data(): Generate item response data from calibration results
#   - compare_eqc_sac(): Compare EQC and SAC calibration results (primary)
#   - compare_eqc_spc(): Deprecated alias for compare_eqc_sac
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
#' # Simulate response data from calibration results
#' sim_data <- simulate_response_data(result = eqc_result, n_persons = 500)
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


#' Simulate Item Response Data from Calibration Results
#'
#' @description
#' Generates item response data using the calibrated parameters from
#' \code{eqc_calibrate()} or \code{sac_calibrate()} (formerly
#' \code{spc_calibrate()}).
#'
#' @param result A calibration result object of class \code{"eqc_result"},
#'   \code{"sac_result"}, or \code{"spc_result"} (for backward compatibility),
#'   as returned by \code{eqc_calibrate()} or \code{sac_calibrate()}.
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
#' # Example 1: Using EQC calibration result
#' eqc_result <- eqc_calibrate(
#'   target_rho = 0.80,
#'   n_items = 25,
#'   model = "rasch",
#'   seed = 42
#' )
#'
#' sim_data <- simulate_response_data(
#'   result = eqc_result,
#'   n_persons = 1000,
#'   latent_shape = "normal",
#'   seed = 123
#' )
#'
#' # Example 2: Using SAC calibration result
#' sac_result <- sac_calibrate(
#'   target_rho = 0.80,
#'   n_items = 25,
#'   model = "rasch",
#'   n_iter = 200,
#'   seed = 42
#' )
#'
#' sim_data2 <- simulate_response_data(
#'   result = sac_result,
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
#' \code{\link{eqc_calibrate}}, \code{\link{sac_calibrate}},
#' \code{\link{compute_reliability_tam}}
#'
#' @export
simulate_response_data <- function(result,
                                   n_persons,
                                   latent_shape = "normal",
                                   latent_params = list(),
                                   seed = NULL) {

  if (!inherits(result, c("eqc_result", "sac_result", "spc_result"))) {
    stop("'result' must be an 'eqc_result' or 'sac_result' object from eqc_calibrate() or sac_calibrate().")
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = globalenv(), inherits = FALSE)
      } else {
        assign(".Random.seed", old_seed, envir = globalenv())
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  # Generate theta
  latent_args <- modifyList(
    list(n = n_persons, shape = latent_shape),
    latent_params
  )
  theta <- do.call(sim_latentG, latent_args)$theta

  beta   <- result$beta_vec
  lambda <- result$lambda_scaled
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
# Comparison Function for EQC and SAC Results
# =============================================================================

#' Compare EQC and SAC Calibration Results
#'
#' @description
#' Compares the calibration results from EQC and SAC algorithms.
#'
#' @param eqc_result An object of class \code{"eqc_result"}.
#' @param sac_result An object of class \code{"sac_result"} (or \code{"spc_result"}
#'   for backward compatibility).
#' @param verbose Logical. If TRUE, print comparison summary.
#'
#' @return A list with comparison statistics (invisibly):
#' \describe{
#'   \item{\code{c_eqc}}{Calibrated c* from EQC.}
#'   \item{\code{c_sac}}{Calibrated c* from SAC.}
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
#' sac_result <- sac_calibrate(target_rho = 0.80, n_items = 25,
#'                             c_init = eqc_result, seed = 42)
#'
#' # Compare results
#' compare_eqc_sac(eqc_result, sac_result)
#' }
#'
#' @seealso
#' \code{\link{eqc_calibrate}}, \code{\link{sac_calibrate}}
#'
#' @export
compare_eqc_sac <- function(eqc_result, sac_result, verbose = TRUE) {

  if (!inherits(eqc_result, "eqc_result")) {
    stop("eqc_result must be an object of class 'eqc_result'.")
  }
  if (!inherits(sac_result, c("sac_result", "spc_result"))) {
    stop("sac_result must be an object of class 'sac_result' (or 'spc_result' for backward compatibility).")
  }

  # Configuration validation
  if (!isTRUE(all.equal(eqc_result$target_rho, sac_result$target_rho))) {
    stop(sprintf(
      "target_rho differs between EQC (%.4f) and SAC (%.4f). Results are not comparable.",
      eqc_result$target_rho, sac_result$target_rho
    ))
  }
  if (!identical(eqc_result$model, sac_result$model)) {
    warning(sprintf(
      "Model differs between EQC ('%s') and SAC ('%s').",
      eqc_result$model, sac_result$model
    ))
  }
  if (!identical(eqc_result$n_items, sac_result$n_items)) {
    warning(sprintf(
      "n_items differs between EQC (%d) and SAC (%d).",
      eqc_result$n_items, sac_result$n_items
    ))
  }
  if (!identical(eqc_result$metric, sac_result$metric)) {
    warning(sprintf(
      "Reliability metric differs between EQC ('%s') and SAC ('%s').",
      eqc_result$metric, sac_result$metric
    ))
  }

  # Extract key values
  c_eqc <- eqc_result$c_star
  c_sac <- sac_result$c_star

  diff_abs <- abs(c_eqc - c_sac)
  diff_pct <- 100 * diff_abs / c_eqc

  agreement <- diff_pct < 5  # Within 5%

  comparison <- list(
    c_eqc = c_eqc,
    c_sac = c_sac,
    diff_abs = diff_abs,
    diff_pct = diff_pct,
    agreement = agreement,
    target_rho = sac_result$target_rho
  )

  if (verbose) {
    message("")
    message("=======================================================")
    message("  EQC vs SAC Comparison")
    message("=======================================================")
    message("")
    message(sprintf("  Target reliability  : %.4f", comparison$target_rho))
    message(sprintf("  EQC c*              : %.6f", c_eqc))
    message(sprintf("  SAC c*              : %.6f", c_sac))
    message(sprintf("  Absolute difference : %.6f", diff_abs))
    message(sprintf("  Percent difference  : %.2f%%", diff_pct))
    message(sprintf("  Agreement (< 5%%)    : %s", ifelse(agreement, "YES", "NO")))
    message("")
  }

  invisible(comparison)
}


#' @rdname compare_eqc_sac
#' @param ... Arguments passed to \code{compare_eqc_sac()}.
#' @export
compare_eqc_spc <- function(...) {
  .Deprecated("compare_eqc_sac")
  compare_eqc_sac(...)
}
