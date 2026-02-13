#' Simulate Item Parameters for IRT Studies
#'
#' @description
#' `sim_item_params()` generates item parameters (difficulty \eqn{\beta} and
#' discrimination \eqn{\lambda}) for Item Response Theory (IRT) simulation studies.
#' It wraps the IRW `irw_simu_diff()` function for realistic difficulty distributions
#' and provides multiple methods for generating correlated discriminations.
#'
#' The function is designed with four key principles:
#' \enumerate{
#'   \item \strong{Realistic difficulties:} Integration with Item Response Warehouse (IRW)
#'     for empirically-grounded difficulty distributions.
#'   \item \strong{Correlated parameters:} Support for the empirically observed negative
#'     correlation between difficulty and discrimination (Sweeney et al., 2022).
#'   \item \strong{Marginal preservation:} Copula method preserves exact marginal
#'     distributions while achieving target correlation.
#'   \item \strong{Reliability targeting:} Scale factor for subsequent calibration.
#' }
#'
#' @param n_items Integer. Number of items to generate per form.
#' @param model Character. The data-generating model: "rasch" or "2pl".
#' @param source Character. Source for generating difficulties:
#'   \describe{
#'     \item{\code{"irw"}}{Use IRW difficulty pool (realistic, empirical)}
#'     \item{\code{"parametric"}}{Generate from parametric distribution}
#'     \item{\code{"hierarchical"}}{Joint MVN for both parameters (Glas & van der Linden)}
#'     \item{\code{"custom"}}{User-supplied parameters or function}
#'   }
#' @param method Character. Method for generating discriminations (when model = "2pl"):
#'   \describe{
#'     \item{\code{"copula"}}{Gaussian copula - preserves marginals exactly (RECOMMENDED)}
#'     \item{\code{"conditional"}}{Conditional normal regression on difficulty}
#'     \item{\code{"independent"}}{Independent generation (no correlation)}
#'   }
#' @param n_forms Integer. Number of test forms to generate. Default is 1.
#'   When > 1, returns a data frame with form_id column.
#' @param difficulty_params List. Parameters for difficulty generation:
#'   \describe{
#'     \item{For \code{source = "irw"}:}{\code{pool} - difficulty pool data frame}
#'     \item{For \code{source = "parametric"}:}{\code{mu}, \code{sigma}, \code{distribution}}
#'   }
#' @param discrimination_params List. Parameters for discrimination generation:
#'   \describe{
#'     \item{\code{mu_log}}{Mean of log-discrimination (default: 0)}
#'     \item{\code{sigma_log}}{SD of log-discrimination (default: 0.3)}
#'     \item{\code{rho}}{Target correlation between \eqn{\beta} and \eqn{\log(\lambda)} (default: -0.3)}
#'   }
#' @param hierarchical_params List. For source = "hierarchical":
#'   \describe{
#'     \item{\code{mu}}{2-vector: means of \eqn{(\log\lambda, \beta)}}
#'     \item{\code{tau}}{2-vector: SDs}
#'     \item{\code{rho}}{Correlation}
#'   }
#' @param custom_params List. For source = "custom":
#'   \describe{
#'     \item{\code{beta}}{Vector or function returning difficulties}
#'     \item{\code{lambda}}{Vector or function returning discriminations}
#'   }
#' @param scale Numeric. Global discrimination scaling factor for reliability targeting.
#'   Final discriminations are \eqn{\lambda_i^* = c \cdot \lambda_i}. Default is 1.
#' @param center_difficulties Logical. If TRUE, center difficulties to sum to zero
#'   for identification. Default is TRUE.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @details
#' ## Why Copula Method is Recommended
#'
#' When difficulties come from the IRW pool (which has realistic, often non-normal

#' marginal distributions), the conditional normal method can distort the achieved
#' correlation because it assumes linearity. The Gaussian copula method:
#'
#' \enumerate{
#'   \item Transforms difficulties to uniform scale via empirical CDF
#'   \item Generates correlated uniforms through Gaussian copula
#'   \item Transforms back to desired marginals (log-normal for discrimination)
#' }
#'
#' This guarantees:
#' \itemize{
#'   \item Exact preservation of difficulty marginal (whatever IRW provides)
#'   \item Exact log-normal marginal for discriminations
#'   \item Spearman correlation \eqn{\approx \rho} (rank-based, robust to non-normality
#' }
#'
#' ## Connection to Reliability-Targeted Framework
#'
#' The \code{scale} parameter implements "separation of structure and scale":
#' \itemize{
#'   \item \strong{Structure}: Realistic item characteristics from IRW + correlation
#'   \item \strong{Scale}: Global factor \eqn{c} calibrated for target reliability
#' }
#'
#' @return An object of class \code{"item_params"} containing:
#' \describe{
#'   \item{\code{data}}{Data frame with columns: form_id, item_id, beta, lambda, lambda_unscaled}
#'   \item{\code{model}}{Model type used}
#'   \item{\code{source}}{Source used for generation}
#'   \item{\code{method}}{Method used for discrimination generation}
#'   \item{\code{n_items}}{Number of items per form}
#'   \item{\code{n_forms}}{Number of forms generated}
#'   \item{\code{scale}}{Scale factor applied}
#'   \item{\code{centered}}{Whether difficulties were centered}
#'   \item{\code{params}}{Parameters used for generation}
#'   \item{\code{achieved}}{Achieved statistics (correlations, moments)}
#' }
#'
#' @examples
#' # Example 1: Rasch with parametric difficulties
#' items1 <- sim_item_params(n_items = 25, model = "rasch",
#'                           source = "parametric", seed = 42)
#'
#' # Example 2: 2PL with copula method (recommended)
#' items2 <- sim_item_params(
#'   n_items = 30, model = "2pl", source = "parametric",
#'   method = "copula",
#'   discrimination_params = list(rho = -0.3),
#'   seed = 42
#' )
#'
#' # Example 3: Hierarchical 2PL
#' items3 <- sim_item_params(
#'   n_items = 25, model = "2pl", source = "hierarchical",
#'   hierarchical_params = list(mu = c(0, 0), tau = c(0.25, 1), rho = -0.3),
#'   seed = 42
#' )
#'
#' \dontrun{
#' # Example 4: Using IRW difficulty pool (requires irw package)
#' items4 <- sim_item_params(n_items = 25, model = "rasch", source = "irw")
#'
#' # Example 5: Multiple forms with IRW
#' items5 <- sim_item_params(
#'   n_items = 20, model = "2pl", n_forms = 5,
#'   source = "irw", method = "copula"
#' )
#' }
#'
#' @references
#' Glas, C. A. W., & van der Linden, W. J. (2003). Computerized adaptive testing
#'   with item cloning. \emph{Applied Psychological Measurement, 27}(4), 247-261.
#'
#' Sweeney, S. M., et al. (2022). An investigation of the nature and consequence
#'   of the relationship between IRT difficulty and discrimination.
#'   \emph{EM:IP, 41}(4), 50-67.
#'
#' Zhang, L., et al. (2025). Realistic simulation of item difficulties. \emph{PsyArXiv}.
#'
#' @export
sim_item_params <- function(n_items,
                            model = c("rasch", "2pl"),
                            source = c("irw", "parametric", "hierarchical", "custom"),
                            method = c("copula", "conditional", "independent"),
                            n_forms = 1L,
                            difficulty_params = list(),
                            discrimination_params = list(),
                            hierarchical_params = list(),
                            custom_params = list(),
                            scale = 1,
                            center_difficulties = TRUE,
                            seed = NULL) {

  # ===========================================================================
  # Input Validation
  # ===========================================================================

  if (!is.numeric(n_items) || length(n_items) != 1L || n_items <= 0) {
    stop("`n_items` must be a positive integer.")
  }
  n_items <- as.integer(n_items)

  model <- match.arg(model)
  source <- match.arg(source)
  method <- match.arg(method)

  if (!is.numeric(n_forms) || length(n_forms) != 1L || n_forms <= 0) {
    stop("`n_forms` must be a positive integer.")
  }
  n_forms <- as.integer(n_forms)

  if (!is.numeric(scale) || length(scale) != 1L || scale <= 0) {
    stop("`scale` must be a positive scalar.")
  }

  # Set seed (save/restore RNG state)
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

  # ===========================================================================
  # Set Default Parameters
  # ===========================================================================

  # Difficulty defaults
  diff_defaults <- list(
    mu = 0, sigma = 1, distribution = "normal", pool = NULL
  )
  difficulty_params <- modifyList(diff_defaults, difficulty_params)

  # Discrimination defaults
  disc_defaults <- list(
    mu_log = 0,        # E[lambda] ≈ exp(0) = 1
    sigma_log = 0.3,   # Moderate variation
    rho = -0.3         # Negative correlation (Sweeney et al., 2022)
  )
  discrimination_params <- modifyList(disc_defaults, discrimination_params)

  # Hierarchical defaults (Glas & van der Linden style)
  hier_defaults <- list(
    mu = c(0, 0),      # (log-lambda, beta)
    tau = c(0.3, 1),   # SDs
    rho = -0.3         # Correlation
  )
  hierarchical_params <- modifyList(hier_defaults, hierarchical_params)

  # ===========================================================================
  # Generate Item Parameters by Source
  # ===========================================================================

  # Container for all forms
  all_data <- vector("list", n_forms)

  for (f in seq_len(n_forms)) {

    result <- switch(source,

                     # -----------------------------------------------------------------------
                     # Source: IRW
                     # -----------------------------------------------------------------------
                     "irw" = {
                       beta <- .generate_irw_difficulties(n_items, difficulty_params)

                       if (model == "2pl") {
                         lambda <- .generate_discriminations(
                           beta = beta,
                           method = method,
                           params = discrimination_params
                         )
                       } else {
                         lambda <- rep(1, n_items)
                       }

                       list(beta = beta, lambda = lambda)
                     },

                     # -----------------------------------------------------------------------
                     # Source: Parametric
                     # -----------------------------------------------------------------------
                     "parametric" = {
                       beta <- .generate_parametric_difficulties(n_items, difficulty_params)

                       if (model == "2pl") {
                         lambda <- .generate_discriminations(
                           beta = beta,
                           method = method,
                           params = discrimination_params
                         )
                       } else {
                         lambda <- rep(1, n_items)
                       }

                       list(beta = beta, lambda = lambda)
                     },

                     # -----------------------------------------------------------------------
                     # Source: Hierarchical (Joint MVN)
                     # -----------------------------------------------------------------------
                     "hierarchical" = {
                       pars <- .generate_hierarchical_2pl(n_items, hierarchical_params)
                       list(beta = pars$beta, lambda = pars$lambda)
                     },

                     # -----------------------------------------------------------------------
                     # Source: Custom
                     # -----------------------------------------------------------------------
                     "custom" = {
                       beta <- .process_custom_param(custom_params$beta, n_items, "beta")

                       if (model == "2pl") {
                         if (is.null(custom_params$lambda)) {
                           stop("For model = '2pl' with source = 'custom', `custom_params$lambda` is required.")
                         }
                         lambda <- .process_custom_param(custom_params$lambda, n_items, "lambda")
                       } else {
                         lambda <- rep(1, n_items)
                       }

                       list(beta = beta, lambda = lambda)
                     }
    )

    beta <- result$beta
    lambda <- result$lambda

    # Center difficulties if requested
    if (center_difficulties) {
      beta <- beta - mean(beta)
    }

    # Store unscaled lambda
    lambda_unscaled <- lambda

    # Apply scaling
    lambda_scaled <- scale * lambda

    # Create data frame for this form
    all_data[[f]] <- data.frame(
      form_id = f,
      item_id = seq_len(n_items),
      beta = beta,
      lambda = lambda_scaled,
      lambda_unscaled = lambda_unscaled
    )
  }

  # Combine all forms
  data <- do.call(rbind, all_data)
  rownames(data) <- NULL

  # ===========================================================================
  # Compute Achieved Statistics
  # ===========================================================================

  achieved <- list()

  # Per-form statistics
  achieved$by_form <- lapply(split(data, data$form_id), function(df) {
    list(
      beta_mean = mean(df$beta),
      beta_sd = sd(df$beta),
      beta_range = range(df$beta),
      lambda_mean = mean(df$lambda),
      lambda_sd = sd(df$lambda),
      lambda_range = range(df$lambda),
      cor_pearson = if (model == "2pl") cor(df$beta, log(df$lambda_unscaled)) else NA,
      cor_spearman = if (model == "2pl") cor(df$beta, log(df$lambda_unscaled), method = "spearman") else NA
    )
  })

  # Overall statistics
  achieved$overall <- list(
    n_total = nrow(data),
    beta_mean = mean(data$beta),
    beta_sd = sd(data$beta),
    lambda_mean = mean(data$lambda),
    lambda_sd = sd(data$lambda),
    cor_pearson_pooled = if (model == "2pl") cor(data$beta, log(data$lambda_unscaled)) else NA,
    cor_spearman_pooled = if (model == "2pl") cor(data$beta, log(data$lambda_unscaled), method = "spearman") else NA
  )

  # ===========================================================================
  # Construct Output Object
  # ===========================================================================

  output <- list(
    data = data,
    model = model,
    source = source,
    method = if (model == "2pl" && source != "hierarchical") method else NA,
    n_items = n_items,
    n_forms = n_forms,
    scale = scale,
    centered = center_difficulties,
    params = list(
      difficulty = difficulty_params,
      discrimination = discrimination_params,
      hierarchical = hierarchical_params
    ),
    achieved = achieved
  )

  class(output) <- c("item_params", "list")

  return(output)
}


# =============================================================================
# Internal Helper Functions
# =============================================================================

#' Generate difficulties from IRW pool
#' @noRd
.generate_irw_difficulties <- function(n_items, params) {

  irw_available <- nzchar(system.file(package = "irw"))
  if (!irw_available) {
    stop("Package 'irw' is required for source = 'irw'. ",
         "Install with: devtools::install_github('itemresponsewarehouse/Rpkg')")
  }

  pool <- params$pool
  if (is.null(pool)) {
    utils::data("diff_long", package = "irw", envir = environment())
    pool <- get("diff_long", envir = environment())
  }

  # Use IRW function with num_replications = 1
  irw_simu_diff <- utils::getFromNamespace("irw_simu_diff", "irw")
  result <- irw_simu_diff(
    num_items = n_items,
    num_replications = 1,
    difficulty_pool = pool
  )

  # Extract difficulty column
  if (is.data.frame(result) && "difficulty" %in% names(result)) {
    beta <- result$difficulty
  } else if (is.matrix(result)) {
    beta <- as.numeric(result[1, ])
  } else {
    beta <- as.numeric(result)
  }

  return(beta)
}


#' Generate difficulties from parametric distribution
#' @noRd
.generate_parametric_difficulties <- function(n_items, params) {

  dist <- params$distribution
  mu <- params$mu
  sigma <- params$sigma

  beta <- switch(dist,
                 "normal" = rnorm(n_items, mean = mu, sd = sigma),
                 "uniform" = {
                   half_range <- sqrt(3) * sigma
                   runif(n_items, min = mu - half_range, max = mu + half_range)
                 },
                 rnorm(n_items, mean = mu, sd = sigma)  # default
  )

  return(beta)
}


#' Generate discriminations using specified method
#' @noRd
.generate_discriminations <- function(beta, method, params) {

  n <- length(beta)
  mu_log <- params$mu_log
  sigma_log <- params$sigma_log
  rho <- params$rho

  # Validate rho
  if (!is.numeric(rho) || abs(rho) > 1) {
    stop("`rho` must be between -1 and 1.")
  }

  lambda <- switch(method,

                   # =========================================================================
                   # Method: Copula (RECOMMENDED)
                   # Preserves marginals exactly, achieves target Spearman correlation
                   # =========================================================================
                   "copula" = {
                     .generate_copula(beta, mu_log, sigma_log, rho)
                   },

                   # =========================================================================
                   # Method: Conditional Normal
                   # log(lambda) | beta ~ N(conditional_mean, conditional_var)
                   # =========================================================================
                   "conditional" = {
                     .generate_conditional(beta, mu_log, sigma_log, rho)
                   },

                   # =========================================================================
                   # Method: Independent
                   # log(lambda) ~ N(mu_log, sigma_log^2) independent of beta
                   # =========================================================================
                   "independent" = {
                     exp(rnorm(n, mean = mu_log, sd = sigma_log))
                   }
  )

  return(lambda)
}


#' Generate discriminations using Gaussian copula
#'
#' This method preserves exact marginals while achieving target correlation.
#'
#' Algorithm:
#' 1. Transform beta to uniform via empirical CDF: u = rank(beta)/(n+1)
#' 2. Transform to normal: z_beta = qnorm(u)
#' 3. Generate correlated normal: z_lambda = rho * z_beta + sqrt(1-rho^2) * z_indep
#' 4. Transform to uniform: v = pnorm(z_lambda)
#' 5. Transform to log-normal: lambda = exp(mu + sigma * qnorm(v))
#'
#' @noRd
.generate_copula <- function(beta, mu_log, sigma_log, rho) {

  n <- length(beta)

  # Step 1: Transform beta to uniform via empirical CDF (rank-based)
  # Using (rank - 0.5)/n for better behavior at boundaries
  u_beta <- (rank(beta, ties.method = "random") - 0.5) / n

  # Step 2: Transform to normal scale
  z_beta <- qnorm(u_beta)

  # Step 3: Generate correlated normal (Gaussian copula)
  z_indep <- rnorm(n)
  z_lambda <- rho * z_beta + sqrt(1 - rho^2) * z_indep

  # Step 4: Transform to uniform
  v_lambda <- pnorm(z_lambda)

  # Step 5: Transform to log-normal via quantile function
  # If log(lambda) ~ N(mu, sigma^2), then lambda = exp(mu + sigma * qnorm(v))
  log_lambda <- mu_log + sigma_log * qnorm(v_lambda)
  lambda <- exp(log_lambda)

  return(lambda)
}


#' Generate discriminations using conditional normal
#' @noRd
.generate_conditional <- function(beta, mu_log, sigma_log, rho) {

  n <- length(beta)

  if (abs(rho) < 1e-10) {
    # Independent case
    return(exp(rnorm(n, mean = mu_log, sd = sigma_log)))
  }

  # Standardize beta
  beta_z <- (beta - mean(beta)) / sd(beta)

  # Conditional mean: mu_log + rho * sigma_log * z_beta
  cond_mean <- mu_log + rho * sigma_log * beta_z

  # Conditional SD: sigma_log * sqrt(1 - rho^2)
  cond_sd <- sigma_log * sqrt(1 - rho^2)

  log_lambda <- rnorm(n, mean = cond_mean, sd = cond_sd)
  lambda <- exp(log_lambda)

  return(lambda)
}


#' Generate parameters from hierarchical 2PL (Glas & van der Linden)
#' @noRd
.generate_hierarchical_2pl <- function(n_items, params) {

  mu <- params$mu
  tau <- params$tau
  rho <- params$rho

  # Construct covariance matrix
  Omega <- matrix(c(1, rho, rho, 1), nrow = 2)
  Sigma <- diag(tau) %*% Omega %*% diag(tau)

  # Draw from bivariate normal
  if (requireNamespace("MASS", quietly = TRUE)) {
    xi <- MASS::mvrnorm(n_items, mu = mu, Sigma = Sigma)
  } else {
    # Fallback using Cholesky
    L <- t(chol(Sigma))
    Z <- matrix(rnorm(n_items * 2), nrow = n_items, ncol = 2)
    xi <- sweep(Z %*% t(L), 2, mu, "+")
  }

  log_lambda <- xi[, 1]
  beta <- xi[, 2]
  lambda <- exp(log_lambda)

  return(list(beta = beta, lambda = lambda))
}


#' Process custom parameter input
#' @noRd
.process_custom_param <- function(x, n, name) {
  if (is.null(x)) {
    stop("custom_params$", name, " is required for source = 'custom'.")
  }

  if (is.function(x)) {
    result <- x(n)
  } else if (is.numeric(x)) {
    if (length(x) != n) {
      stop("custom_params$", name, " must have length n_items = ", n)
    }
    result <- x
  } else {
    stop("custom_params$", name, " must be numeric vector or function.")
  }

  return(as.numeric(result))
}


# =============================================================================
# S3 Methods
# =============================================================================

#' @rdname sim_item_params
#' @param x An object of class \code{"item_params"}.
#' @param ... Additional arguments passed to or from other methods.
#' @return The input object, invisibly.
#' @export
print.item_params <- function(x, ...) {
  cat("Item Parameters Object\n")
  cat("======================\n")
  cat(sprintf("  Model          : %s\n", toupper(x$model)))
  cat(sprintf("  Source         : %s\n", x$source))
  if (!is.na(x$method)) {
    cat(sprintf("  Method         : %s\n", x$method))
  }
  cat(sprintf("  Items per form : %d\n", x$n_items))
  cat(sprintf("  Number of forms: %d\n", x$n_forms))
  cat(sprintf("  Scale factor   : %.3f\n", x$scale))
  cat(sprintf("  Centered       : %s\n", ifelse(x$centered, "Yes", "No")))

  cat("\nDifficulty (beta):\n")
  cat(sprintf("  Mean: %.4f, SD: %.4f, Range: [%.3f, %.3f]\n",
              x$achieved$overall$beta_mean,
              x$achieved$overall$beta_sd,
              min(x$data$beta), max(x$data$beta)))

  if (x$model == "2pl") {
    cat("\nDiscrimination (lambda, scaled):\n")
    cat(sprintf("  Mean: %.4f, SD: %.4f, Range: [%.3f, %.3f]\n",
                x$achieved$overall$lambda_mean,
                x$achieved$overall$lambda_sd,
                min(x$data$lambda), max(x$data$lambda)))
    cat(sprintf("\nCorrelation (beta, log-lambda):\n"))
    cat(sprintf("  Target (rho): %.3f\n", x$params$discrimination$rho))
    cat(sprintf("  Achieved Pearson : %.3f\n", x$achieved$overall$cor_pearson_pooled))
    cat(sprintf("  Achieved Spearman: %.3f\n", x$achieved$overall$cor_spearman_pooled))
  }

  invisible(x)
}


#' @rdname sim_item_params
#' @param object An object of class \code{"item_params"}.
#' @return An object of class \code{"summary.item_params"} containing key
#'   parameter summaries.
#' @export
summary.item_params <- function(object, ...) {

  beta_vals <- object$data$beta
  beta_summary <- list(
    mean = mean(beta_vals),
    sd   = sd(beta_vals),
    min  = min(beta_vals),
    max  = max(beta_vals),
    quantiles = quantile(beta_vals, probs = c(0.25, 0.50, 0.75))
  )

  lambda_summary <- NULL
  if (object$model == "2pl") {
    lambda_vals <- object$data$lambda
    lambda_unscaled_vals <- object$data$lambda_unscaled
    lambda_summary <- list(
      mean_unscaled = mean(lambda_unscaled_vals),
      sd_unscaled   = sd(lambda_unscaled_vals),
      mean_scaled   = mean(lambda_vals),
      sd_scaled     = sd(lambda_vals),
      min_scaled    = min(lambda_vals),
      max_scaled    = max(lambda_vals)
    )
  }

  achieved_cors <- NULL
  if (object$model == "2pl") {
    achieved_cors <- list(
      target_rho       = object$params$discrimination$rho,
      pearson_pooled   = object$achieved$overall$cor_pearson_pooled,
      spearman_pooled  = object$achieved$overall$cor_spearman_pooled
    )
  }

  out <- list(
    model          = object$model,
    source         = object$source,
    method         = object$method,
    n_items        = object$n_items,
    n_forms        = object$n_forms,
    scale          = object$scale,
    centered       = object$centered,
    beta_summary   = beta_summary,
    lambda_summary = lambda_summary,
    achieved_cors  = achieved_cors
  )
  class(out) <- "summary.item_params"
  out
}


#' Print Method for summary.item_params Objects
#'
#' @param x A \code{summary.item_params} object from \code{summary.item_params()}.
#' @param digits Integer. Number of decimal places for printing.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @export
print.summary.item_params <- function(x, digits = 4, ...) {
  cat("Summary: Item Parameters\n")
  cat("========================\n")
  cat(sprintf("  Model          : %s\n", toupper(x$model)))
  cat(sprintf("  Source         : %s\n", x$source))
  if (!is.na(x$method)) {
    cat(sprintf("  Method         : %s\n", x$method))
  }
  cat(sprintf("  Items per form : %d\n", x$n_items))
  cat(sprintf("  Number of forms: %d\n", x$n_forms))
  cat(sprintf("  Scale factor   : %.3f\n", x$scale))
  cat(sprintf("  Centered       : %s\n", ifelse(x$centered, "Yes", "No")))

  cat("\nDifficulty (beta):\n")
  cat(sprintf("  Mean     : %.*f\n", digits, x$beta_summary$mean))
  cat(sprintf("  SD       : %.*f\n", digits, x$beta_summary$sd))
  cat(sprintf("  Min      : %.*f\n", digits, x$beta_summary$min))
  cat(sprintf("  Max      : %.*f\n", digits, x$beta_summary$max))
  cat(sprintf("  Quantiles: Q25=%.*f, Q50=%.*f, Q75=%.*f\n",
              digits, x$beta_summary$quantiles[[1]],
              digits, x$beta_summary$quantiles[[2]],
              digits, x$beta_summary$quantiles[[3]]))

  if (x$model == "2pl" && !is.null(x$lambda_summary)) {
    cat("\nDiscrimination (lambda):\n")
    cat(sprintf("  Before scaling: Mean=%.*f, SD=%.*f\n",
                digits, x$lambda_summary$mean_unscaled,
                digits, x$lambda_summary$sd_unscaled))
    cat(sprintf("  After scaling (c=%.3f): Mean=%.*f, SD=%.*f\n",
                x$scale, digits, x$lambda_summary$mean_scaled,
                digits, x$lambda_summary$sd_scaled))
    cat(sprintf("  Range [%.*f, %.*f]\n",
                digits, x$lambda_summary$min_scaled,
                digits, x$lambda_summary$max_scaled))
  }

  if (x$model == "2pl" && !is.null(x$achieved_cors)) {
    cat("\nCorrelation (beta, log-lambda):\n")
    cat(sprintf("  Target (rho)     : %.*f\n", digits, x$achieved_cors$target_rho))
    cat(sprintf("  Achieved Pearson : %.*f\n", digits, x$achieved_cors$pearson_pooled))
    cat(sprintf("  Achieved Spearman: %.*f\n", digits, x$achieved_cors$spearman_pooled))
  }

  invisible(x)
}


#' Plot method for item_params objects
#'
#' @param x An `item_params` object from `sim_item_params()`.
#' @param type Character. Type of plot: `"scatter"`, `"density"`, or `"both"`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A \code{ggplot} object if \pkg{ggplot2} is available, or
#'   \code{NULL} invisibly when using base R graphics fallback.
#'
#' @export
plot.item_params <- function(x, type = c("scatter", "density", "both"), ...) {

  type <- match.arg(type)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Base R fallback
    if (x$model == "2pl") {
      oldpar <- par(mfrow = c(1, 2))
      on.exit(par(oldpar))
      hist(x$data$beta, main = "Difficulty Distribution", xlab = "beta", col = "steelblue")
      plot(x$data$beta, x$data$lambda, main = "Difficulty vs Discrimination",
           xlab = "beta", ylab = "lambda", pch = 19, col = rgb(0, 0, 0.5, 0.5))
      abline(lm(lambda ~ beta, data = x$data), col = "red", lty = 2)
    } else {
      hist(x$data$beta, main = "Difficulty Distribution", xlab = "beta", col = "steelblue")
    }
    return(invisible(NULL))
  }

  df <- x$data

  # Scatter plot
  p_scatter <- NULL
  if (x$model == "2pl") {
    rho_achieved <- x$achieved$overall$cor_spearman_pooled
    p_scatter <- ggplot2::ggplot(df, ggplot2::aes(x = beta, y = lambda)) +
      ggplot2::geom_point(alpha = 0.6, color = "darkblue", size = 2) +
      ggplot2::geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
      ggplot2::labs(
        title = sprintf("Difficulty vs Discrimination (Spearman r = %.3f)", rho_achieved),
        subtitle = sprintf("Method: %s | Target rho: %.2f", x$method, x$params$discrimination$rho),
        x = expression(beta ~ "(Difficulty)"),
        y = expression(lambda ~ "(Discrimination)")
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  }

  # Density plots
  p_beta <- ggplot2::ggplot(df, ggplot2::aes(x = beta)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 30,
                   fill = "steelblue", alpha = 0.5, color = "white") +
    ggplot2::geom_density(color = "darkblue", linewidth = 1) +
    ggplot2::labs(title = "Difficulty Distribution", x = expression(beta), y = "Density") +
    ggplot2::theme_minimal()

  p_lambda <- NULL
  if (x$model == "2pl") {
    p_lambda <- ggplot2::ggplot(df, ggplot2::aes(x = lambda)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 30,
                     fill = "coral", alpha = 0.5, color = "white") +
      ggplot2::geom_density(color = "darkred", linewidth = 1) +
      ggplot2::labs(title = "Discrimination Distribution", x = expression(lambda), y = "Density") +
      ggplot2::theme_minimal()
  }

  # Return based on type
  if (type == "scatter") {
    if (x$model == "2pl") return(p_scatter) else return(p_beta)
  } else if (type == "density") {
    if (x$model == "2pl" && requireNamespace("patchwork", quietly = TRUE)) {
      return(patchwork::wrap_plots(p_beta, p_lambda, ncol = 2))
    } else {
      return(p_beta)
    }
  } else {  # both
    if (x$model == "2pl" && requireNamespace("patchwork", quietly = TRUE)) {
      combined_top <- patchwork::wrap_plots(p_beta, p_lambda, ncol = 2)
      return(patchwork::wrap_plots(combined_top, p_scatter, ncol = 1))
    } else {
      return(p_beta)
    }
  }
}


#' Extract item parameters as data frame
#' @param x An `item_params` object.
#' @param row.names NULL or a character vector giving the row names.
#' @param optional Logical. If TRUE, setting row names is optional.
#' @param ... Additional arguments (ignored).
#' @return A data frame of item parameters.
#' @export
as.data.frame.item_params <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$data
}
