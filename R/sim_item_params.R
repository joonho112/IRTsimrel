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
#' # Example 1: Rasch with IRW difficulties
#' items1 <- sim_item_params(n_items = 25, model = "rasch", source = "irw")
#'
#' # Example 2: 2PL with copula method (recommended)
#' items2 <- sim_item_params(
#'   n_items = 30, model = "2pl", source = "irw",
#'   method = "copula",
#'   discrimination_params = list(rho = -0.3)
#' )
#'
#' # Example 3: Multiple forms
#' items3 <- sim_item_params(
#'   n_items = 20, model = "2pl", n_forms = 5,
#'   source = "irw", method = "copula"
#' )
#'
#' # Example 4: Hierarchical 2PL
#' items4 <- sim_item_params(
#'   n_items = 25, model = "2pl", source = "hierarchical",
#'   hierarchical_params = list(mu = c(0, 0), tau = c(0.25, 1), rho = -0.3)
#' )
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

  # Set seed
  if (!is.null(seed)) {
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

  if (!requireNamespace("irw", quietly = TRUE)) {
    stop("Package 'irw' is required for source = 'irw'. ",
         "Install with: devtools::install_github('itemresponsewarehouse/Rpkg')")
  }

  pool <- params$pool
  if (is.null(pool)) {
    data("diff_long", package = "irw", envir = environment())
    pool <- get("diff_long", envir = environment())
  }

  # Use IRW function with num_replications = 1
  result <- irw::irw_simu_diff(
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


#' @export
summary.item_params <- function(object, ...) {

  cat("="
      , rep("=", 60), "\n", sep = "")
  cat("Summary of Item Parameters\n")
  cat(rep("=", 61), "\n\n", sep = "")

  cat("Generation Settings:\n")
  cat(sprintf("  Model          : %s\n", toupper(object$model)))
  cat(sprintf("  Source         : %s\n", object$source))
  if (!is.na(object$method)) {
    cat(sprintf("  Method         : %s\n", object$method))
  }
  cat(sprintf("  Items/form     : %d\n", object$n_items))
  cat(sprintf("  Forms          : %d\n", object$n_forms))
  cat(sprintf("  Total items    : %d\n", nrow(object$data)))
  cat(sprintf("  Scale factor   : %.3f\n", object$scale))
  cat(sprintf("  Centered       : %s\n\n", ifelse(object$centered, "Yes", "No")))

  # Difficulty summary
  cat("Difficulty (beta):\n")
  cat(sprintf("  Mean     : %8.4f\n", object$achieved$overall$beta_mean))
  cat(sprintf("  SD       : %8.4f\n", object$achieved$overall$beta_sd))
  cat(sprintf("  Min      : %8.4f\n", min(object$data$beta)))
  cat(sprintf("  Max      : %8.4f\n", max(object$data$beta)))
  cat(sprintf("  Quantiles: Q25=%.3f, Q50=%.3f, Q75=%.3f\n\n",
              quantile(object$data$beta, 0.25),
              quantile(object$data$beta, 0.50),
              quantile(object$data$beta, 0.75)))

  if (object$model == "2pl") {
    cat("Discrimination (lambda):\n")
    cat("  Before scaling:\n")
    cat(sprintf("    Mean   : %8.4f\n", mean(object$data$lambda_unscaled)))
    cat(sprintf("    SD     : %8.4f\n", sd(object$data$lambda_unscaled)))
    cat(sprintf("  After scaling (c = %.3f):\n", object$scale))
    cat(sprintf("    Mean   : %8.4f\n", object$achieved$overall$lambda_mean))
    cat(sprintf("    SD     : %8.4f\n", object$achieved$overall$lambda_sd))
    cat(sprintf("    Min    : %8.4f\n", min(object$data$lambda)))
    cat(sprintf("    Max    : %8.4f\n\n", max(object$data$lambda)))

    cat("Correlation (beta, log-lambda):\n")
    cat(sprintf("  Target (rho)    : %8.4f\n", object$params$discrimination$rho))
    cat(sprintf("  Achieved Pearson: %8.4f\n", object$achieved$overall$cor_pearson_pooled))
    cat(sprintf("  Achieved Spearman:%8.4f\n", object$achieved$overall$cor_spearman_pooled))
  }

  invisible(object)
}


#' Plot method for item_params objects
#' @export
plot.item_params <- function(x, type = c("scatter", "density", "both"), ...) {

  type <- match.arg(type)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Base R fallback
    if (x$model == "2pl") {
      par(mfrow = c(1, 2))
      hist(x$data$beta, main = "Difficulty Distribution", xlab = "beta", col = "steelblue")
      plot(x$data$beta, x$data$lambda, main = "Difficulty vs Discrimination",
           xlab = "beta", ylab = "lambda", pch = 19, col = rgb(0, 0, 0.5, 0.5))
      abline(lm(lambda ~ beta, data = x$data), col = "red", lty = 2)
      par(mfrow = c(1, 1))
    } else {
      hist(x$data$beta, main = "Difficulty Distribution", xlab = "beta", col = "steelblue")
    }
    return(invisible(NULL))
  }

  library(ggplot2)

  df <- x$data

  # Scatter plot
  p_scatter <- NULL
  if (x$model == "2pl") {
    rho_achieved <- x$achieved$overall$cor_spearman_pooled
    p_scatter <- ggplot(df, aes(x = beta, y = lambda)) +
      geom_point(alpha = 0.6, color = "darkblue", size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
      labs(
        title = sprintf("Difficulty vs Discrimination (Spearman r = %.3f)", rho_achieved),
        subtitle = sprintf("Method: %s | Target rho: %.2f", x$method, x$params$discrimination$rho),
        x = expression(beta ~ "(Difficulty)"),
        y = expression(lambda ~ "(Discrimination)")
      ) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
  }

  # Density plots
  p_beta <- ggplot(df, aes(x = beta)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30,
                   fill = "steelblue", alpha = 0.5, color = "white") +
    geom_density(color = "darkblue", linewidth = 1) +
    labs(title = "Difficulty Distribution", x = expression(beta), y = "Density") +
    theme_minimal()

  p_lambda <- NULL
  if (x$model == "2pl") {
    p_lambda <- ggplot(df, aes(x = lambda)) +
      geom_histogram(aes(y = after_stat(density)), bins = 30,
                     fill = "coral", alpha = 0.5, color = "white") +
      geom_density(color = "darkred", linewidth = 1) +
      labs(title = "Discrimination Distribution", x = expression(lambda), y = "Density") +
      theme_minimal()
  }

  # Return based on type
  if (type == "scatter") {
    if (x$model == "2pl") return(p_scatter) else return(p_beta)
  } else if (type == "density") {
    if (x$model == "2pl" && requireNamespace("patchwork", quietly = TRUE)) {
      return(p_beta + p_lambda)
    } else {
      return(p_beta)
    }
  } else {  # both
    if (x$model == "2pl" && requireNamespace("patchwork", quietly = TRUE)) {
      library(patchwork)
      return((p_beta + p_lambda) / p_scatter)
    } else {
      return(p_beta)
    }
  }
}


#' Extract item parameters as data frame
#' @export
as.data.frame.item_params <- function(x, ...) {
  x$data
}
