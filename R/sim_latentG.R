#' Simulate Latent Ability Distribution for IRT Studies (G-family)
#'
#' @description
#' `sim_latentG()` generates latent abilities (person parameters \eqn{\theta}) for
#' Item Response Theory (IRT) simulation studies. It implements the population
#' model \eqn{\theta_p \sim G} where \eqn{G} is a flexible distribution family.
#'
#' The function is designed with two key principles:
#' \enumerate{
#'   \item \strong{Pre-standardization:} Each distribution shape is mathematically
#'     constructed to have mean 0 and variance 1, ensuring that changing the
#'     shape does not inadvertently change the scale.
#'   \item \strong{Separation of Structure and Scale:} The `sigma` parameter
#'     directly controls the standard deviation of the latent trait, independent
#'     of the distributional shape.
#' }
#'
#' The generated abilities follow:
#' \deqn{\theta_p = \mu + X_p^\top \beta + \sigma \cdot z_p}
#' where \eqn{z_p \sim G_0} with \eqn{E[z]=0} and \eqn{Var[z]=1}.
#'
#' @param n Integer. Number of persons (latent abilities) to generate.
#'
#' @param shape Character. The distributional shape of the standardized
#'   component. One of:
#'   \describe{
#'     \item{\code{"normal"}}{Standard normal \eqn{N(0,1)}}
#'     \item{\code{"bimodal"}}{Symmetric two-component Gaussian mixture with
#'       analytically standardized parameters}
#'     \item{\code{"trimodal"}}{Symmetric three-component Gaussian mixture}
#'     \item{\code{"multimodal"}}{Four-component Gaussian mixture}
#'     \item{\code{"skew_pos"}}{Right-skewed distribution via standardized Gamma}
#'     \item{\code{"skew_neg"}}{Left-skewed distribution (negated Gamma)}
#'     \item{\code{"heavy_tail"}}{Heavy-tailed distribution via standardized Student-t}
#'     \item{\code{"light_tail"}}{Light-tailed (platykurtic) mixture distribution}
#'     \item{\code{"uniform"}}{Uniform distribution on \eqn{[-\sqrt{3}, \sqrt{3}]}}
#'     \item{\code{"floor"}}{Distribution with floor effect (left-truncated feel)}
#'     \item{\code{"ceiling"}}{Distribution with ceiling effect (right-truncated feel)}
#'     \item{\code{"custom"}}{User-specified mixture distribution}
#'   }
#'
#' @param sigma Numeric. Scale (standard deviation) of the residual latent trait.
#'   Since the standardized component has variance 1, `sigma` directly equals
#'   the marginal SD of the residual term. Default is 1.
#'
#' @param mu Numeric. Grand mean of the latent ability distribution. In Rasch
#'   models this is often fixed to 0 for identification. Default is 0.
#'
#' @param xcov Matrix or data.frame. Optional covariate matrix with `n` rows.
#'   If supplied, person-specific covariate effects are added as \eqn{\eta = X\beta}.
#'
#' @param beta Numeric vector. Regression coefficients for `xcov`. Must have
#'   length equal to `ncol(xcov)`. Ignored if `xcov` is NULL.
#'
#' @param shape_params List. Additional parameters controlling the shape.
#'   See Details for shape-specific parameters.
#'
#' @param mixture_spec List. For `shape = "custom"`, specifies the mixture:
#'   \describe{
#'     \item{\code{weights}}{Numeric vector of mixing proportions (must sum to 1)}
#'     \item{\code{means}}{Numeric vector of component means}
#'     \item{\code{sds}}{Numeric vector of component standard deviations}
#'   }
#'   The custom mixture is automatically standardized to have mean 0 and variance 1.
#'
#' @param standardize_custom Logical. If TRUE (default), custom mixtures are
#'   post-standardized to ensure mean 0 and variance 1. If FALSE, the raw
#'   mixture is used (user must ensure proper standardization).
#'
#' @param seed Integer. Random seed for reproducibility. If NULL (default),
#'   the current RNG state is used.
#'
#' @param return_z Logical. If TRUE, include the standardized draws `z` in
#'   the output. Default is TRUE.
#'
#' @details
#' ## Pre-standardization Mathematics
#'
#' Each built-in shape is constructed to have exactly mean 0 and variance 1:
#'
#' \strong{Bimodal:} Two-component mixture with modes at \eqn{\pm\delta}:
#' \deqn{z = s \cdot \delta + \epsilon, \quad s \sim \text{Rademacher}, \quad \epsilon \sim N(0, 1-\delta^2)}
#' where the component variance \eqn{1-\delta^2} ensures \eqn{Var[z] = \delta^2 + (1-\delta^2) = 1}.
#'
#' \strong{Trimodal:} Three-component mixture with weights \eqn{(w_L, w_0, w_R)}
#' and means \eqn{(-m, 0, m)}. Component variance is \eqn{\sigma_c^2 = 1 - (1-w_0)m^2}
#' to ensure unit total variance.
#'
#' \strong{Skewed:} Standardized Gamma distribution:
#' \deqn{z = \frac{\Gamma(k, 1) - k}{\sqrt{k}}}
#' which has \eqn{E[z]=0} and \eqn{Var[z]=1} for any \eqn{k > 0}.
#'
#' \strong{Heavy-tailed:} Standardized Student-t:
#' \deqn{z = \frac{t_\nu}{\sqrt{\nu/(\nu-2)}}}
#' which has \eqn{Var[z]=1} for \eqn{\nu > 2}.
#'
#' ## Shape-Specific Parameters
#'
#' \describe{
#'   \item{\code{delta}}{For "bimodal": mode separation, must satisfy
#'     \eqn{0 < \delta < 1}. Default: 0.8}
#'   \item{\code{w0}}{For "trimodal": weight of central component,
#'     must satisfy \eqn{0 < w_0 < 1}. Default: 1/3}
#'   \item{\code{m}}{For "trimodal": magnitude of side component means.
#'     Must satisfy \eqn{(1-w_0)m^2 < 1}. Default: 1.2}
#'   \item{\code{k}}{For "skew_pos"/"skew_neg": Gamma shape parameter,
#'     controls skewness magnitude. Default: 4}
#'   \item{\code{df}}{For "heavy_tail": degrees of freedom, must be > 2.
#'     Default: 5}
#' }
#'
#' ## Connection to IRT Framework
#'
#' In the Rasch/2PL model, the latent distribution \eqn{G} affects:
#' \itemize{
#'   \item Marginal reliability: \eqn{\bar{w} = \sigma_\theta^2 / (\sigma_\theta^2 + \text{MSEM})}
#'   \item Expected test information: \eqn{\bar{\mathcal{J}} = E_G[\mathcal{J}(\theta)]}
#'   \item Identifiability (see Appendix F of the manuscript)
#' }
#'
#' This function serves as the generator for \eqn{G} in reliability-targeted
#' simulation studies, allowing researchers to examine how distributional
#' shape affects model performance while holding scale constant.
#'
#' @return An object of class `"latent_G"` (a list) containing:
#' \describe{
#'   \item{\code{theta}}{Numeric vector of length `n`, the simulated latent abilities}
#'   \item{\code{z}}{Standardized draws (if `return_z = TRUE`)}
#'   \item{\code{eta_cov}}{Covariate linear predictor (0 if no covariates)}
#'   \item{\code{mu}}{Grand mean used}
#'   \item{\code{sigma}}{Scale parameter used}
#'   \item{\code{shape}}{Shape label}
#'   \item{\code{shape_params}}{Shape parameters used}
#'   \item{\code{n}}{Sample size}
#'   \item{\code{sample_moments}}{List with sample mean, sd, skewness, kurtosis}
#' }
#'
#' @examples
#' # Basic usage: standard normal abilities
#' sim1 <- sim_latentG(n = 1000, shape = "normal")
#' mean(sim1$theta)  # approximately 0
#' sd(sim1$theta)    # approximately 1
#'
#' # Bimodal distribution for heterogeneous population
#' sim2 <- sim_latentG(n = 1000, shape = "bimodal",
#'                     shape_params = list(delta = 0.9))
#'
#' # Skewed distribution with larger scale
#' sim3 <- sim_latentG(n = 1000, shape = "skew_pos", sigma = 1.5)
#'
#' # With covariate effects (e.g., group differences)
#' group <- rbinom(1000, 1, 0.5)
#' sim4 <- sim_latentG(n = 1000, shape = "normal",
#'                     xcov = data.frame(group = group),
#'                     beta = 0.5)
#'
#' # Custom mixture distribution
#' sim5 <- sim_latentG(n = 1000, shape = "custom",
#'                     mixture_spec = list(
#'                       weights = c(0.3, 0.5, 0.2),
#'                       means = c(-1.5, 0, 2),
#'                       sds = c(0.5, 0.7, 0.5)
#'                     ))
#'
#' @seealso
#' \code{\link{summary.latent_G}} for summary statistics,
#' \code{\link{plot.latent_G}} for visualization,
#' \code{\link{compare_shapes}} for comparing multiple shapes.
#'
#' @references
#' Baker, F. B., & Kim, S.-H. (2004). \emph{Item Response Theory: Parameter
#'   Estimation Techniques} (2nd ed.). Marcel Dekker.
#'
#' Paganin, S., et al. (2022). Computational strategies and estimation
#'   performance with Bayesian semiparametric item response theory models.
#'   \emph{Journal of Educational and Behavioral Statistics, 48}(2), 147-188.
#'
#' @export
sim_latentG <- function(n,
                        shape = c("normal", "bimodal", "trimodal", "multimodal",
                                  "skew_pos", "skew_neg", "heavy_tail",
                                  "light_tail", "uniform", "floor", "ceiling",
                                  "custom"),
                        sigma = 1,
                        mu = 0,
                        xcov = NULL,
                        beta = NULL,
                        shape_params = list(),
                        mixture_spec = NULL,
                        standardize_custom = TRUE,
                        seed = NULL,
                        return_z = TRUE) {

  # ===========================================================================
  # Input Validation
  # ===========================================================================

  # Validate n

  if (!is.numeric(n) || length(n) != 1L || n <= 0 || n != as.integer(n)) {
    stop("`n` must be a positive integer.")
  }
  n <- as.integer(n)

  # Match and validate shape
  shape <- match.arg(shape)

  # Validate sigma
  if (!is.numeric(sigma) || length(sigma) != 1L || sigma <= 0) {
    stop("`sigma` must be a positive scalar.")
  }

  # Validate mu
  if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu)) {
    stop("`mu` must be a finite numeric scalar.")
  }

  # Validate shape_params
  if (!is.list(shape_params)) {
    stop("`shape_params` must be a list.")
  }

  # Validate custom mixture specification
  if (shape == "custom") {
    if (is.null(mixture_spec) || !is.list(mixture_spec)) {
      stop("For shape = 'custom', `mixture_spec` must be provided as a list.")
    }
    required_fields <- c("weights", "means", "sds")
    if (!all(required_fields %in% names(mixture_spec))) {
      stop("`mixture_spec` must contain 'weights', 'means', and 'sds'.")
    }
    if (abs(sum(mixture_spec$weights) - 1) > 1e-10) {
      stop("mixture_spec$weights must sum to 1.")
    }
    K <- length(mixture_spec$weights)
    if (length(mixture_spec$means) != K || length(mixture_spec$sds) != K) {
      stop("All components of mixture_spec must have the same length.")
    }
    if (any(mixture_spec$sds <= 0)) {
      stop("All mixture_spec$sds must be positive.")
    }
  }

  # Set seed if provided
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L) {
      stop("`seed` must be a single numeric value.")
    }
    set.seed(as.integer(seed))
  }

  # ===========================================================================
  # Generate Standardized Component z ~ G0 with E[z]=0, Var[z]=1
  # ===========================================================================

  z <- switch(shape,

              # --- Normal: N(0,1) ---
              "normal" = {
                rnorm(n, mean = 0, sd = 1)
              },

              # --- Bimodal: Symmetric two-component mixture ---
              # z = s * delta + epsilon, where s in {-1, +1}, epsilon ~ N(0, 1-delta^2)
              # Var[z] = delta^2 + (1-delta^2) = 1
              "bimodal" = {
                delta <- shape_params$delta %||% 0.8
                if (!is.numeric(delta) || length(delta) != 1L ||
                    delta <= 0 || delta >= 1) {
                  stop("For 'bimodal', `delta` must satisfy 0 < delta < 1.")
                }
                sd_comp <- sqrt(1 - delta^2)
                sign_vec <- sample(c(-1, 1), size = n, replace = TRUE)
                rnorm(n, mean = sign_vec * delta, sd = sd_comp)
              },

              # --- Trimodal: Three-component mixture ---
              # Components at -m, 0, +m with weights (w_L, w_0, w_R) where w_L = w_R
              # Var[z] = (1-w0)*m^2 + sigma_c^2 = 1, so sigma_c^2 = 1 - (1-w0)*m^2
              "trimodal" = {
                w0 <- shape_params$w0 %||% (1/3)
                m <- shape_params$m %||% 1.2

                if (!is.numeric(w0) || length(w0) != 1L || w0 <= 0 || w0 >= 1) {
                  stop("For 'trimodal', `w0` must satisfy 0 < w0 < 1.")
                }
                if (!is.numeric(m) || length(m) != 1L || m <= 0) {
                  stop("For 'trimodal', `m` must be positive.")
                }

                # Check variance constraint
                var_from_means <- (1 - w0) * m^2
                if (var_from_means >= 1) {
                  stop("For 'trimodal', must have (1 - w0) * m^2 < 1. ",
                       "Current value: ", round(var_from_means, 3))
                }

                sd_comp <- sqrt(1 - var_from_means)
                w_side <- (1 - w0) / 2
                probs <- c(w_side, w0, w_side)
                comp <- sample(c("left", "mid", "right"), size = n,
                               replace = TRUE, prob = probs)
                means <- ifelse(comp == "mid", 0,
                                ifelse(comp == "left", -m, m))
                rnorm(n, mean = means, sd = sd_comp)
              },

              # --- Multimodal: Four-component mixture ---
              # Symmetric: components at -m2, -m1, +m1, +m2
              "multimodal" = {
                m1 <- shape_params$m1 %||% 0.5
                m2 <- shape_params$m2 %||% 1.3
                w_inner <- shape_params$w_inner %||% 0.30

                w_outer <- (1 - 2*w_inner) / 2
                if (w_outer < 0) {
                  stop("For 'multimodal', w_inner must be < 0.5")
                }

                # Var = 2*w_inner*m1^2 + 2*w_outer*m2^2 + sigma_c^2
                # We need: 2*w_inner*m1^2 + 2*w_outer*m2^2 < 1
                var_from_means <- 2*w_inner*m1^2 + 2*w_outer*m2^2
                if (var_from_means >= 1) {
                  stop("For 'multimodal', mean variance contribution must be < 1. ",
                       "Try reducing m1, m2, or adjusting w_inner. Current: ",
                       round(var_from_means, 3))
                }
                sd_comp <- sqrt(1 - var_from_means)

                probs <- c(w_outer, w_inner, w_inner, w_outer)
                comp <- sample(1:4, size = n, replace = TRUE, prob = probs)
                means <- c(-m2, -m1, m1, m2)[comp]
                rnorm(n, mean = means, sd = sd_comp)
              },

              # --- Skew positive: Standardized Gamma ---
              # z = (Gamma(k,1) - k) / sqrt(k) has E[z]=0, Var[z]=1
              "skew_pos" = {
                k <- shape_params$k %||% 4
                if (!is.numeric(k) || length(k) != 1L || k <= 0) {
                  stop("For 'skew_pos', `k` must be positive.")
                }
                g <- rgamma(n, shape = k, rate = 1)
                (g - k) / sqrt(k)
              },

              # --- Skew negative: Negated standardized Gamma ---
              "skew_neg" = {
                k <- shape_params$k %||% 4
                if (!is.numeric(k) || length(k) != 1L || k <= 0) {
                  stop("For 'skew_neg', `k` must be positive.")
                }
                g <- rgamma(n, shape = k, rate = 1)
                -(g - k) / sqrt(k)
              },

              # --- Heavy-tailed: Standardized Student-t ---
              # z = t_df / sqrt(df/(df-2)) has Var[z]=1 for df > 2
              "heavy_tail" = {
                df <- shape_params$df %||% 5
                if (!is.numeric(df) || length(df) != 1L || df <= 2) {
                  stop("For 'heavy_tail', `df` must be > 2.")
                }
                t_raw <- rt(n, df = df)
                t_raw / sqrt(df / (df - 2))
              },

              # --- Light-tailed: Mixture approximating platykurtic ---
              # Four-point mixture approximating uniform-ish but with smoother tails
              "light_tail" = {
                # Symmetric 4-component mixture
                probs <- rep(0.25, 4)
                means <- c(-1.3, -0.43, 0.43, 1.3)
                # Var = sum(w * mu^2) + sigma_c^2
                # = 0.25*(1.69 + 0.1849 + 0.1849 + 1.69) + sigma_c^2
                # = 0.9375 + sigma_c^2 = 1 => sigma_c^2 = 0.0625
                sd_comp <- sqrt(0.0625)  # = 0.25
                comp <- sample(1:4, size = n, replace = TRUE, prob = probs)
                rnorm(n, mean = means[comp], sd = sd_comp)
              },

              # --- Uniform: U(-sqrt(3), sqrt(3)) ---
              # Has E[z]=0, Var[z]=1
              "uniform" = {
                runif(n, min = -sqrt(3), max = sqrt(3))
              },

              # --- Floor effect: Left-skewed mixture ---
              "floor" = {
                # Heavy component near lower bound
                w_floor <- shape_params$w_floor %||% 0.3
                m_floor <- shape_params$m_floor %||% (-1.5)

                w_main <- 1 - w_floor
                # E[z] = w_floor * m_floor + w_main * m_main = 0
                m_main <- -w_floor * m_floor / w_main

                # Var = w_floor*(m_floor^2 + s_f^2) + w_main*(m_main^2 + s_m^2) = 1
                # Simplify: use same sd_comp
                var_from_means <- w_floor * m_floor^2 + w_main * m_main^2
                if (var_from_means >= 1) {
                  stop("Floor effect parameters result in invalid variance.")
                }
                sd_comp <- sqrt((1 - var_from_means))

                is_floor <- runif(n) < w_floor
                means <- ifelse(is_floor, m_floor, m_main)
                rnorm(n, mean = means, sd = sd_comp)
              },

              # --- Ceiling effect: Right-skewed mixture ---
              "ceiling" = {
                w_ceil <- shape_params$w_ceil %||% 0.3
                m_ceil <- shape_params$m_ceil %||% 1.5

                w_main <- 1 - w_ceil
                m_main <- -w_ceil * m_ceil / w_main

                var_from_means <- w_ceil * m_ceil^2 + w_main * m_main^2
                if (var_from_means >= 1) {
                  stop("Ceiling effect parameters result in invalid variance.")
                }
                sd_comp <- sqrt(1 - var_from_means)

                is_ceil <- runif(n) < w_ceil
                means <- ifelse(is_ceil, m_ceil, m_main)
                rnorm(n, mean = means, sd = sd_comp)
              },

              # --- Custom mixture ---
              "custom" = {
                weights <- mixture_spec$weights
                means <- mixture_spec$means
                sds <- mixture_spec$sds
                K <- length(weights)

                # Sample component memberships
                comp <- sample(1:K, size = n, replace = TRUE, prob = weights)
                z_raw <- rnorm(n, mean = means[comp], sd = sds[comp])

                # Optionally standardize
                if (standardize_custom) {
                  z_mean <- sum(weights * means)
                  z_var <- sum(weights * (sds^2 + means^2)) - z_mean^2
                  (z_raw - z_mean) / sqrt(z_var)
                } else {
                  z_raw
                }
              },

              stop("Unknown shape: ", shape)
  )

  # ===========================================================================
  # Covariate Effects
  # ===========================================================================

  if (is.null(xcov)) {
    eta_cov <- rep(0, n)
  } else {
    xmat <- as.matrix(xcov)
    if (nrow(xmat) != n) {
      stop("`xcov` must have ", n, " rows.")
    }
    if (is.null(beta)) {
      stop("`beta` must be provided when `xcov` is not NULL.")
    }
    if (!is.numeric(beta) || length(beta) != ncol(xmat)) {
      stop("`beta` must be numeric with length equal to ncol(xcov).")
    }
    eta_cov <- as.numeric(xmat %*% beta)
  }

  # ===========================================================================
  # Assemble Final Latent Abilities
  # ===========================================================================

  theta <- mu + eta_cov + sigma * z

  # ===========================================================================
  # Compute Sample Moments
  # ===========================================================================

  theta_mean <- mean(theta)
  theta_sd <- sd(theta)
  z_std <- (theta - theta_mean) / theta_sd
  skewness <- mean(z_std^3)
  kurtosis <- mean(z_std^4) - 3  # Excess kurtosis

  # ===========================================================================
  # Construct Output Object
  # ===========================================================================

  out <- list(
    theta = theta,
    mu = mu,
    sigma = sigma,
    eta_cov = eta_cov,
    shape = shape,
    shape_params = shape_params,
    n = n,
    sample_moments = list(
      mean = theta_mean,
      sd = theta_sd,
      skewness = skewness,
      kurtosis = kurtosis
    )
  )

  if (return_z) {
    out$z <- z
  }

  if (shape == "custom") {
    out$mixture_spec <- mixture_spec
  }

  class(out) <- c("latent_G", "list")
  return(out)
}


# =============================================================================
# Helper: Null-coalescing operator
# =============================================================================
`%||%` <- function(x, y) if (is.null(x)) y else x


# =============================================================================
# S3 Methods
# =============================================================================

#' Print Method for latent_G Objects
#'
#' @param x A `latent_G` object from `sim_latentG()`.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.latent_G <- function(x, ...) {
  cat("Latent Ability Distribution (G-family)\n")
  cat("=======================================\n")
  cat(sprintf("  Shape     : %s\n", x$shape))
  cat(sprintf("  n         : %d\n", x$n))
  cat(sprintf("  Target mu : %.3f\n", x$mu))
  cat(sprintf("  Target sigma: %.3f\n", x$sigma))
  cat("\nSample Moments:\n")
  cat(sprintf("  Mean      : %.4f\n", x$sample_moments$mean))
  cat(sprintf("  SD        : %.4f\n", x$sample_moments$sd))
  cat(sprintf("  Skewness  : %.4f\n", x$sample_moments$skewness))
  cat(sprintf("  Kurtosis  : %.4f (excess)\n", x$sample_moments$kurtosis))
  invisible(x)
}


#' Summary Method for latent_G Objects
#'
#' @param object A `latent_G` object from `sim_latentG()`.
#' @param ... Additional arguments (ignored).
#'
#' @return A list with summary statistics.
#'
#' @export
summary.latent_G <- function(object, ...) {
  theta <- object$theta

  quants <- quantile(theta, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))

  out <- list(
    shape = object$shape,
    n = object$n,
    target_mu = object$mu,
    target_sigma = object$sigma,
    sample_mean = mean(theta),
    sample_sd = sd(theta),
    sample_median = median(theta),
    skewness = object$sample_moments$skewness,
    kurtosis = object$sample_moments$kurtosis,
    min = min(theta),
    max = max(theta),
    quantiles = quants,
    has_covariates = !all(object$eta_cov == 0)
  )

  class(out) <- "summary.latent_G"
  return(out)
}


#' @export
print.summary.latent_G <- function(x, digits = 4, ...) {
  cat("Summary: Latent Ability Distribution\n")
  cat("====================================\n")
  cat(sprintf("  Shape      : %s\n", x$shape))
  cat(sprintf("  n          : %d\n", x$n))
  cat(sprintf("  Target     : mu = %.2f, sigma = %.2f\n",
              x$target_mu, x$target_sigma))
  cat(sprintf("  Covariates : %s\n", ifelse(x$has_covariates, "Yes", "No")))
  cat("\nSample Statistics:\n")
  cat(sprintf("  Mean       : %.*f\n", digits, x$sample_mean))
  cat(sprintf("  SD         : %.*f\n", digits, x$sample_sd))
  cat(sprintf("  Median     : %.*f\n", digits, x$sample_median))
  cat(sprintf("  Skewness   : %.*f\n", digits, x$skewness))
  cat(sprintf("  Kurtosis   : %.*f (excess)\n", digits, x$kurtosis))
  cat(sprintf("  Range      : [%.*f, %.*f]\n", digits, x$min, digits, x$max))
  cat("\nQuantiles:\n")
  print(round(x$quantiles, digits))
  invisible(x)
}


#' Plot Method for latent_G Objects
#'
#' @param x A `latent_G` object from `sim_latentG()`.
#' @param type Character. Type of plot: "histogram", "density", or "both".
#' @param show_normal Logical. Overlay a normal reference density?
#' @param bins Integer. Number of histogram bins.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A ggplot2 object (if available) or base R plot.
#'
#' @export
plot.latent_G <- function(x,
                          type = c("both", "histogram", "density"),
                          show_normal = TRUE,
                          bins = 50,
                          ...) {
  type <- match.arg(type)
  theta <- x$theta
  mu_sample <- x$sample_moments$mean
  sd_sample <- x$sample_moments$sd

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    df <- data.frame(theta = theta)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = theta))

    if (type %in% c("histogram", "both")) {
      p <- p + ggplot2::geom_histogram(
        ggplot2::aes(y = ggplot2::after_stat(density)),
        bins = bins,
        fill = "steelblue",
        color = "white",
        alpha = 0.7
      )
    }

    if (type %in% c("density", "both")) {
      p <- p + ggplot2::geom_density(
        color = "darkblue",
        linewidth = 1
      )
    }

    if (show_normal) {
      p <- p + ggplot2::stat_function(
        fun = dnorm,
        args = list(mean = mu_sample, sd = sd_sample),
        color = "darkred",
        linewidth = 0.8,
        linetype = "dashed"
      )
    }

    p <- p +
      ggplot2::labs(
        title = paste0("Latent Ability Distribution (", x$shape, ")"),
        subtitle = sprintf("n = %d, μ = %.3f, σ = %.3f",
                           x$n, mu_sample, sd_sample),
        x = expression(theta ~ "(Latent Ability)"),
        y = "Density"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

    return(p)

  } else {
    # Base R fallback
    hist(theta, breaks = bins, freq = FALSE,
         main = paste0("Latent Distribution (", x$shape, ")"),
         xlab = expression(theta),
         col = "steelblue", border = "white")
    if (show_normal) {
      curve(dnorm(x, mean = mu_sample, sd = sd_sample),
            add = TRUE, col = "darkred", lwd = 2, lty = 2)
    }
    invisible(NULL)
  }
}


#' Compare Multiple Distribution Shapes
#'
#' @description
#' Generates and compares multiple latent ability distributions side-by-side.
#'
#' @param n Integer. Sample size for each distribution.
#' @param shapes Character vector. Shapes to compare.
#' @param sigma Numeric. Common scale parameter.
#' @param seed Integer. Random seed.
#'
#' @return A ggplot2 object with faceted density plots.
#'
#' @examples
#' \dontrun{
#' compare_shapes(n = 3000)
#' }
#'
#' @export
compare_shapes <- function(n = 2000,
                           shapes = c("normal", "bimodal", "trimodal",
                                      "skew_pos", "skew_neg", "heavy_tail",
                                      "uniform"),
                           sigma = 1,
                           seed = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.")
  }

  if (!is.null(seed)) set.seed(seed)

  # Generate all distributions
  all_data <- lapply(shapes, function(s) {
    sim <- sim_latentG(n = n, shape = s, sigma = sigma)
    data.frame(
      theta = sim$theta,
      shape = s,
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, all_data)

  # Format shape names
  df$shape <- factor(df$shape, levels = shapes)

  # Create faceted plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = theta)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 45,
      fill = "steelblue",
      color = "white",
      alpha = 0.7
    ) +
    ggplot2::stat_function(
      fun = dnorm,
      args = list(mean = 0, sd = sigma),
      color = "darkred",
      linewidth = 0.7,
      linetype = "dashed"
    ) +
    ggplot2::facet_wrap(~ shape, ncol = 3, scales = "free_y") +
    ggplot2::labs(
      title = "Comparison of Latent Ability Distribution Shapes",
      subtitle = sprintf("n = %d per shape; dashed = N(0, %.1f) reference", n, sigma),
      x = expression(theta ~ "(Latent Ability)"),
      y = "Density"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      strip.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::coord_cartesian(xlim = c(-4*sigma, 4*sigma))

  return(p)
}


#' Extract theta values from latent_G object
#'
#' @description
#' Convenience function to extract the theta vector for use with other functions.
#'
#' @param x A `latent_G` object.
#'
#' @return Numeric vector of latent abilities.
#'
#' @export
as.numeric.latent_G <- function(x) {
  x$theta
}
