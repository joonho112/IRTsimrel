#' Simulate Item Parameters for IRT Studies
#'
#' @description
#' `sim_item_params()` generates item parameters (difficulty \eqn{\beta} and
#' discrimination \eqn{\lambda}, and optional lower asymptote \eqn{g}) for
#' Item Response Theory (IRT) simulation studies.
#' It provides parametric, hierarchical, custom, and optional Item Response
#' Warehouse (IRW) sources, plus multiple methods for generating correlated
#' discriminations.
#'
#' The function is designed with five key principles:
#' \enumerate{
#'   \item \strong{Portable defaults:} Parametric item generation works without
#'     external item-pool data.
#'   \item \strong{Empirical integration:} Optional IRW integration is available
#'     when the Item Response Warehouse package is installed.
#'   \item \strong{Correlated parameters:} Support for the empirically observed negative
#'     correlation between difficulty and discrimination (Sweeney et al., 2022).
#'   \item \strong{Marginal control:} The copula method leaves sampled
#'     difficulties unchanged and approximately targets the requested
#'     discrimination marginal and dependence in finite forms.
#'   \item \strong{Reliability targeting:} Scale factor for subsequent calibration.
#' }
#'
#' @param n_items Integer. Number of items to generate per form.
#' @param model Character. The data-generating model: "rasch", "2pl", or
#'   "3pl". The 3PL uses a fixed logistic scaling constant of \eqn{D=1}.
#' @param source Character. Source for generating difficulties:
#'   \describe{
#'     \item{\code{"parametric"}}{Generate from a parametric difficulty distribution (default).}
#'     \item{\code{"irw"}}{Use an IRW difficulty pool when the external IRW package is installed.}
#'     \item{\code{"hierarchical"}}{Joint MVN for both parameters (Glas & van der Linden)}
#'     \item{\code{"custom"}}{User-supplied parameters or function}
#'   }
#' @param method Character. Method for generating discriminations (when model
#'   is "2pl" or "3pl"):
#'   \describe{
#'     \item{\code{"copula"}}{Rank-based Gaussian-copula construction
#'       (recommended; finite-form marginal and correlation are approximate).}
#'     \item{\code{"conditional"}}{Conditional normal regression on difficulty}
#'     \item{\code{"independent"}}{Independent generation (no correlation)}
#'   }
#' @param n_forms Integer. Number of test forms to generate. Default is 1.
#'   When > 1, returns a data frame with form_id column.
#' @param difficulty_params List. Parameters for difficulty generation:
#'   \describe{
#'     \item{For \code{source = "irw"}:}{\code{pool} - optional difficulty pool data frame.}
#'     \item{For \code{source = "parametric"}:}{\code{mu} - finite mean
#'       (default 0), \code{sigma} - positive finite SD (default 1), and
#'       \code{distribution} - \code{"normal"} or \code{"uniform"} (default
#'       \code{"normal"}).}
#'   }
#' @param discrimination_params List. Parameters for discrimination generation:
#'   \describe{
#'     \item{\code{mu_log}}{Finite mean of log-discrimination (default: 0).}
#'     \item{\code{sigma_log}}{Positive finite SD of log-discrimination (default: 0.3).}
#'     \item{\code{rho}}{Finite latent-Gaussian dependence parameter for
#'       \eqn{\beta} and \eqn{\log(\lambda)} in the closed interval from -1
#'       to 1 (default: -0.3). Realized finite-form Pearson and Spearman
#'       correlations are stochastic and need not equal this value.}
#'   }
#' @param hierarchical_params List. For source = "hierarchical":
#'   \describe{
#'     \item{\code{mu}}{Finite 2-vector: means of \eqn{(\log\lambda, \beta)}.}
#'     \item{\code{tau}}{Positive finite 2-vector of SDs.}
#'     \item{\code{rho}}{Finite correlation in (-1, 1).}
#'   }
#' @param custom_params List. For source = "custom":
#'   \describe{
#'     \item{\code{beta}}{Finite numeric vector of length \code{n_items}, or
#'       function returning one.}
#'     \item{\code{lambda}}{Positive finite numeric vector of length
#'       \code{n_items}, or function returning one; required for
#'       \code{model = "2pl"} or \code{"3pl"}.}
#'   }
#' @param scale Numeric. Global discrimination scaling factor for reliability targeting.
#'   Final discriminations are \eqn{\lambda_i^* = c \cdot \lambda_i}. Default is 1.
#' @param center_difficulties Logical. If TRUE, center difficulties to sum to zero
#'   for identification. Default is TRUE.
#' @param seed Integer. Random seed for reproducibility.
#' @param guessing_params List. Lower-asymptote generator for
#'   \code{model = "3pl"}. Set \code{distribution} to \code{"fixed"}
#'   (with scalar or item-length \code{value}, default 0.20), \code{"beta"}
#'   (with positive \code{shape1} and \code{shape2}, defaults 5 and 17), or
#'   \code{"uniform"} (with \code{min} and \code{max}, defaults 0.10 and
#'   0.30). Values must be finite and satisfy \eqn{0 \le g_i < 1}. For a
#'   custom source, supply \code{custom_params$guessing} instead.
#'
#' @details
#' ## Why the Copula Method is Recommended
#'
#' When difficulties come from realistic or otherwise non-normal marginal
#' distributions, the conditional normal method can distort the achieved
#' correlation because it assumes linearity. The Gaussian copula method:
#'
#' \enumerate{
#'   \item Transforms difficulties to uniform scale via empirical CDF
#'   \item Generates correlated uniforms through Gaussian copula
#'   \item Transforms back to desired marginals (log-normal for discrimination)
#' }
#'
#' Its finite-sample properties are:
#' \itemize{
#'   \item Exact preservation of the sampled difficulty marginal
#'   \item Approximate log-normal discrimination marginal, converging toward
#'     the requested marginal as the form size grows
#'   \item The supplied \eqn{\rho} is a latent Gaussian dependence parameter;
#'     realized Pearson and Spearman correlations are stochastic and need not
#'     equal \eqn{\rho}
#' }
#'
#' ## Connection to Reliability-Targeted Framework
#'
#' The \code{scale} parameter implements "separation of structure and scale":
#' \itemize{
#'   \item \strong{Structure}: Item characteristics and optional difficulty-discrimination correlation
#'   \item \strong{Scale}: Global factor \eqn{c} calibrated for target reliability
#' }
#'
#' @return An object of class \code{"item_params"} containing:
#' \describe{
#'   \item{\code{data}}{Data frame with columns: form_id, item_id, beta,
#'     lambda, lambda_unscaled, plus guessing for the 3PL.}
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
#' if (requireNamespace("irw", quietly = TRUE)) {
#'   items4 <- sim_item_params(n_items = 25, model = "rasch", source = "irw")
#'
#'   # Example 5: Multiple forms with IRW
#'   items5 <- sim_item_params(
#'     n_items = 20, model = "2pl", n_forms = 5,
#'     source = "irw", method = "copula"
#'   )
#' }
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
                            model = c("rasch", "2pl", "3pl"),
                            source = c("parametric", "irw", "hierarchical", "custom"),
                            method = c("copula", "conditional", "independent"),
                            n_forms = 1L,
                            difficulty_params = list(),
                            discrimination_params = list(),
                            hierarchical_params = list(),
                            custom_params = list(),
                            scale = 1,
                            center_difficulties = TRUE,
                            seed = NULL,
                            guessing_params = list()) {

  # ===========================================================================
  # Input Validation
  # ===========================================================================

  n_items <- .irtsimrel_validate_positive_integer_scalar(n_items, "n_items")

  model <- match.arg(model)
  source <- match.arg(source)
  method <- match.arg(method)

  n_forms <- .irtsimrel_validate_positive_integer_scalar(n_forms, "n_forms")

  if (!is.numeric(scale) || length(scale) != 1L ||
      is.na(scale) || !is.finite(scale) || scale <= 0) {
    stop("`scale` must be a positive finite scalar.")
  }
  if (!is.logical(center_difficulties) || length(center_difficulties) != 1L ||
      is.na(center_difficulties)) {
    stop("`center_difficulties` must be TRUE or FALSE.")
  }

  difficulty_params <- .irtsimrel_check_list_arg(
    difficulty_params, "difficulty_params"
  )
  discrimination_params <- .irtsimrel_check_list_arg(
    discrimination_params, "discrimination_params"
  )
  hierarchical_params <- .irtsimrel_check_list_arg(
    hierarchical_params, "hierarchical_params"
  )
  custom_params <- .irtsimrel_check_list_arg(custom_params, "custom_params")
  guessing_params <- .irtsimrel_check_list_arg(
    guessing_params, "guessing_params"
  )
  guessing_spec <- .normalize_guessing_params(
    guessing_params = guessing_params,
    model = model,
    source = source,
    custom_params = custom_params,
    n_items = n_items
  )

  restore_seed <- .irtsimrel_set_seed(seed)
  if (!is.null(restore_seed)) on.exit(restore_seed(), add = TRUE)

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

                       if (model %in% c("2pl", "3pl")) {
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

                       if (model %in% c("2pl", "3pl")) {
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
                       lambda <- if (model == "rasch") {
                         rep(1, n_items)
                       } else {
                         pars$lambda
                       }
                       list(beta = pars$beta, lambda = lambda)
                     },

                     # -----------------------------------------------------------------------
                     # Source: Custom
                     # -----------------------------------------------------------------------
                     "custom" = {
                       beta <- .process_custom_param(custom_params$beta, n_items, "beta")

                       if (model %in% c("2pl", "3pl")) {
                         if (is.null(custom_params$lambda)) {
                           stop(
                             "For model = '", model,
                             "' with source = 'custom', `custom_params$lambda` is required."
                           )
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
    guessing <- if (model == "3pl") {
      if (source == "custom") {
        .process_custom_guessing(custom_params$guessing, n_items)
      } else {
        .generate_guessing(n_items, guessing_spec)
      }
    } else {
      NULL
    }

    if (!is.numeric(beta) || length(beta) != n_items ||
        anyNA(beta) || any(!is.finite(beta))) {
      stop("Generated `beta` values must be a finite numeric vector of length n_items.")
    }
    if (!is.numeric(lambda) || length(lambda) != n_items ||
        anyNA(lambda) || any(!is.finite(lambda)) || any(lambda <= 0)) {
      stop("Generated `lambda` values must be a positive finite numeric vector of length n_items.")
    }
    if (model == "3pl") {
      guessing <- .validate_guessing_values(guessing, n_items)
    }

    # Center difficulties if requested
    if (center_difficulties) {
      beta <- beta - mean(beta)
    }

    # Store unscaled lambda
    lambda_unscaled <- lambda

    # Apply scaling
    lambda_scaled <- scale * lambda

    # Create data frame for this form
    form_data <- data.frame(
      form_id = f,
      item_id = seq_len(n_items),
      beta = beta,
      lambda = lambda_scaled,
      lambda_unscaled = lambda_unscaled
    )
    if (model == "3pl") {
      form_data$guessing <- guessing
    }
    all_data[[f]] <- form_data
  }

  # Combine all forms
  data <- do.call(rbind, all_data)
  rownames(data) <- NULL

  if (model == "3pl" && any(data$guessing >= 0.5)) {
    warning(
      "Some `guessing` values are at least 0.5; this is allowed but may be ",
      "impractical for typical multiple-choice items.",
      call. = FALSE
    )
  }

  # ===========================================================================
  # Compute Achieved Statistics
  # ===========================================================================

  achieved <- .irtsimrel_recompute_item_achieved(data, model)

  # ===========================================================================
  # Construct Output Object
  # ===========================================================================

  params_used <- list(
    difficulty = difficulty_params,
    discrimination = discrimination_params,
    hierarchical = hierarchical_params
  )
  if (source == "custom") {
    params_used$custom <- custom_params
  }
  if (model == "3pl") {
    params_used$guessing <- guessing_spec
  }

  output <- list(
    data = data,
    model = model,
    source = source,
    method = if (model %in% c("2pl", "3pl") &&
                 source %in% c("parametric", "irw")) {
      method
    } else {
      NA_character_
    },
    n_items = n_items,
    n_forms = n_forms,
    scale = scale,
    centered = center_difficulties,
    params = params_used,
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

  if (!is.character(dist) || length(dist) != 1L ||
      is.na(dist) || !(dist %in% c("normal", "uniform"))) {
    stop("`difficulty_params$distribution` must be 'normal' or 'uniform'.")
  }
  if (!is.numeric(mu) || length(mu) != 1L ||
      is.na(mu) || !is.finite(mu)) {
    stop("`difficulty_params$mu` must be a finite numeric scalar.")
  }
  if (!is.numeric(sigma) || length(sigma) != 1L ||
      is.na(sigma) || !is.finite(sigma) || sigma <= 0) {
    stop("`difficulty_params$sigma` must be a positive finite numeric scalar.")
  }

  beta <- switch(dist,
                 "normal" = rnorm(n_items, mean = mu, sd = sigma),
                 "uniform" = {
                   half_range <- sqrt(3) * sigma
                   runif(n_items, min = mu - half_range, max = mu + half_range)
                 }
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

  if (!is.numeric(mu_log) || length(mu_log) != 1L ||
      is.na(mu_log) || !is.finite(mu_log)) {
    stop("`discrimination_params$mu_log` must be a finite numeric scalar.")
  }
  if (!is.numeric(sigma_log) || length(sigma_log) != 1L ||
      is.na(sigma_log) || !is.finite(sigma_log) || sigma_log <= 0) {
    stop("`discrimination_params$sigma_log` must be a positive finite numeric scalar.")
  }
  if (!is.numeric(rho) || length(rho) != 1L ||
      is.na(rho) || !is.finite(rho) || abs(rho) > 1) {
    stop("`discrimination_params$rho` must be a finite numeric scalar between -1 and 1.")
  }

  lambda <- switch(method,

                   # =========================================================================
                   # Method: Copula (RECOMMENDED)
                   # Keeps sampled beta values fixed; finite-form lambda marginal
                   # and realized dependence are approximate and stochastic.
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


#' Normalize a 3PL guessing generator specification
#' @noRd
.normalize_guessing_params <- function(guessing_params,
                                       model,
                                       source,
                                       custom_params,
                                       n_items) {
  has_generator <- length(guessing_params) > 0L
  has_custom_guessing <- !is.null(custom_params$guessing)

  if (model != "3pl") {
    if (has_generator || has_custom_guessing) {
      stop(
        "Guessing inputs are only supported when `model = \"3pl\"`.",
        call. = FALSE
      )
    }
    return(NULL)
  }

  if (source == "custom") {
    if (has_generator && has_custom_guessing) {
      stop(
        "Supply custom guessing through `custom_params$guessing` only; ",
        "using nonempty `guessing_params` at the same time is ambiguous.",
        call. = FALSE
      )
    }
    if (has_generator) {
      stop(
        "For `source = \"custom\"`, `guessing_params` must be empty and ",
        "guessing must be supplied as `custom_params$guessing`.",
        call. = FALSE
      )
    }
    return(list(distribution = "custom"))
  }

  if (has_custom_guessing) {
    stop(
      "`custom_params$guessing` is only supported when `source = \"custom\"`.",
      call. = FALSE
    )
  }

  distribution <- if (is.null(guessing_params$distribution)) {
    "fixed"
  } else {
    guessing_params$distribution
  }
  if (!is.character(distribution) || length(distribution) != 1L ||
      is.na(distribution) ||
      !(distribution %in% c("fixed", "beta", "uniform"))) {
    stop(
      "`guessing_params$distribution` must be one of 'fixed', 'beta', or ",
      "'uniform'.",
      call. = FALSE
    )
  }

  allowed <- switch(
    distribution,
    fixed = c("distribution", "value"),
    beta = c("distribution", "shape1", "shape2"),
    uniform = c("distribution", "min", "max")
  )
  unknown <- setdiff(names(guessing_params), allowed)
  if (length(unknown) > 0L) {
    stop(
      "Unknown or inapplicable `guessing_params` field(s) for distribution ",
      "'", distribution, "': ",
      .irtsimrel_backtick_collapse(unknown), ".",
      call. = FALSE
    )
  }

  if (distribution == "fixed") {
    value <- if (is.null(guessing_params$value)) 0.20 else guessing_params$value
    value <- .validate_guessing_values(value, n_items, recycle_scalar = TRUE)
    return(list(distribution = distribution, value = value))
  }

  if (distribution == "beta") {
    shape1 <- if (is.null(guessing_params$shape1)) 5 else guessing_params$shape1
    shape2 <- if (is.null(guessing_params$shape2)) 17 else guessing_params$shape2
    for (entry in c("shape1", "shape2")) {
      value <- get(entry)
      if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
          !is.finite(value) || value <= 0) {
        stop(
          "`guessing_params$", entry,
          "` must be a positive finite numeric scalar.",
          call. = FALSE
        )
      }
    }
    return(list(distribution = distribution, shape1 = shape1, shape2 = shape2))
  }

  min_value <- if (is.null(guessing_params$min)) 0.10 else guessing_params$min
  max_value <- if (is.null(guessing_params$max)) 0.30 else guessing_params$max
  if (!is.numeric(min_value) || length(min_value) != 1L ||
      is.na(min_value) || !is.finite(min_value) ||
      !is.numeric(max_value) || length(max_value) != 1L ||
      is.na(max_value) || !is.finite(max_value) ||
      min_value < 0 || min_value >= max_value || max_value >= 1) {
    stop(
      "`guessing_params$min` and `$max` must satisfy 0 <= min < max < 1.",
      call. = FALSE
    )
  }
  list(distribution = distribution, min = min_value, max = max_value)
}


#' Validate lower-asymptote values
#' @noRd
.validate_guessing_values <- function(x,
                                      n_items,
                                      recycle_scalar = FALSE) {
  if (!is.numeric(x) || length(x) == 0L || anyNA(x) || any(!is.finite(x))) {
    stop("`guessing` must contain finite numeric values.", call. = FALSE)
  }
  if (recycle_scalar && length(x) == 1L) {
    x <- rep(as.numeric(x), n_items)
  }
  if (length(x) != n_items) {
    stop(
      "`guessing` must have length 1 or n_items = ", n_items, ".",
      call. = FALSE
    )
  }
  x <- as.numeric(x)
  if (any(x < 0 | x >= 1)) {
    stop("`guessing` values must satisfy 0 <= guessing < 1.", call. = FALSE)
  }
  x
}


#' Generate 3PL lower-asymptote parameters
#' @noRd
.generate_guessing <- function(n_items, params) {
  switch(
    params$distribution,
    fixed = params$value,
    beta = stats::rbeta(n_items, shape1 = params$shape1, shape2 = params$shape2),
    uniform = stats::runif(n_items, min = params$min, max = params$max)
  )
}


#' Generate discriminations using Gaussian copula
#'
#' This rank-based construction leaves the sampled difficulties unchanged and
#' approximately targets the requested discrimination marginal and dependence.
#'
#' Algorithm:
#' 1. Transform beta to uniform scores: u = (rank(beta)-0.5)/n
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

  if (n == 1L) {
    return(exp(rnorm(1L, mean = mu_log, sd = sigma_log)))
  }

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

  if (!is.numeric(mu) || length(mu) != 2L ||
      anyNA(mu) || any(!is.finite(mu))) {
    stop("`hierarchical_params$mu` must be a finite numeric vector of length 2.")
  }
  if (!is.numeric(tau) || length(tau) != 2L ||
      anyNA(tau) || any(!is.finite(tau)) || any(tau <= 0)) {
    stop("`hierarchical_params$tau` must be a positive finite numeric vector of length 2.")
  }
  if (!is.numeric(rho) || length(rho) != 1L ||
      is.na(rho) || !is.finite(rho) || abs(rho) >= 1) {
    stop("`hierarchical_params$rho` must be a finite numeric scalar in (-1, 1).")
  }

  # Construct covariance matrix
  Omega <- matrix(c(1, rho, rho, 1), nrow = 2)
  Sigma <- diag(tau) %*% Omega %*% diag(tau)

  # Draw from bivariate normal
  if (requireNamespace("MASS", quietly = TRUE)) {
    xi <- MASS::mvrnorm(n_items, mu = mu, Sigma = Sigma)
    if (is.null(dim(xi))) {
      xi <- matrix(xi, nrow = 1L, ncol = 2L)
    }
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

  result <- as.numeric(result)
  if (length(result) != n || anyNA(result) || any(!is.finite(result))) {
    stop(
      "custom_params$", name,
      " must return a finite numeric vector of length n_items = ", n
    )
  }
  if (name == "lambda" && any(result <= 0)) {
    stop("custom_params$lambda must contain positive values.")
  }

  return(result)
}


#' Process custom 3PL lower-asymptote input
#' @noRd
.process_custom_guessing <- function(x, n) {
  if (is.null(x)) {
    stop(
      "custom_params$guessing is required for model = '3pl' with ",
      "source = 'custom'.",
      call. = FALSE
    )
  }

  result <- if (is.function(x)) x(n) else x
  if (!is.numeric(result)) {
    stop(
      "custom_params$guessing must be a numeric scalar, numeric vector, ",
      "or function.",
      call. = FALSE
    )
  }
  tryCatch(
    .validate_guessing_values(result, n, recycle_scalar = TRUE),
    error = function(e) {
      stop("custom_params$guessing: ", conditionMessage(e), call. = FALSE)
    }
  )
}


# =============================================================================
# S3 Methods
# =============================================================================

.item_target_rho <- function(x) {
  if (!(x$model %in% c("2pl", "3pl"))) {
    return(NA_real_)
  }
  if (identical(x$source, "hierarchical")) {
    return(x$params$hierarchical$rho)
  }
  if (x$source %in% c("parametric", "irw")) {
    return(x$params$discrimination$rho)
  }
  NA_real_
}

#' @rdname sim_item_params
#' @param x An object of class \code{"item_params"}.
#' @param digits Integer. Number of decimal places for printing.
#' @param ... Additional arguments passed to or from other methods.
#' @return The input object, invisibly.
#' @export
print.item_params <- function(x, digits = 4, ...) {
  cat("Item Parameters Object\n")
  cat("======================\n")
  cat(sprintf("  Model          : %s\n", toupper(x$model)))
  cat(sprintf("  Source         : %s\n", x$source))
  if (!is.na(x$method)) {
    cat(sprintf("  Method         : %s\n", x$method))
  }
  cat(sprintf("  Items per form : %d\n", x$n_items))
  cat(sprintf("  Number of forms: %d\n", x$n_forms))
  cat(sprintf("  Scale factor   : %.*f\n", digits, x$scale))
  cat(sprintf("  Centered       : %s\n", ifelse(x$centered, "Yes", "No")))

  cat("\nDifficulty (beta):\n")
  cat(sprintf("  Mean: %.*f, SD: %.*f, Range: [%.*f, %.*f]\n",
              digits, x$achieved$overall$beta_mean,
              digits, x$achieved$overall$beta_sd,
              digits, min(x$data$beta), digits, max(x$data$beta)))

  if (x$model %in% c("2pl", "3pl")) {
    cat("\nDiscrimination (lambda, scaled):\n")
    cat(sprintf("  Mean: %.*f, SD: %.*f, Range: [%.*f, %.*f]\n",
                digits, x$achieved$overall$lambda_mean,
                digits, x$achieved$overall$lambda_sd,
                digits, min(x$data$lambda), digits, max(x$data$lambda)))
    cat(sprintf("\nCorrelation (beta, log-lambda):\n"))
    target_rho <- .item_target_rho(x)
    if (!is.na(target_rho)) {
      cat(sprintf("  Target (rho): %.*f\n", digits, target_rho))
    }
    cat(sprintf("  Achieved Pearson : %.*f\n", digits, x$achieved$overall$cor_pearson_pooled))
    cat(sprintf("  Achieved Spearman: %.*f\n", digits, x$achieved$overall$cor_spearman_pooled))
  }

  if (x$model == "3pl") {
    cat("\nGuessing (lower asymptote):\n")
    cat(sprintf(
      "  Mean: %.*f, SD: %.*f, Range: [%.*f, %.*f]\n",
      digits, x$achieved$overall$guessing_mean,
      digits, x$achieved$overall$guessing_sd,
      digits, min(x$data$guessing), digits, max(x$data$guessing)
    ))
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
  if (object$model %in% c("2pl", "3pl")) {
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
  if (object$model %in% c("2pl", "3pl")) {
    achieved_cors <- list(
      target_rho       = .item_target_rho(object),
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
  if (object$model == "3pl") {
    guessing_vals <- object$data$guessing
    out$guessing_summary <- list(
      mean = mean(guessing_vals),
      sd = stats::sd(guessing_vals),
      min = min(guessing_vals),
      max = max(guessing_vals),
      quantiles = stats::quantile(
        guessing_vals,
        probs = c(0.25, 0.50, 0.75)
      )
    )
  }
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

  if (x$model %in% c("2pl", "3pl") && !is.null(x$lambda_summary)) {
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

  if (x$model %in% c("2pl", "3pl") && !is.null(x$achieved_cors)) {
    cat("\nCorrelation (beta, log-lambda):\n")
    if (!is.na(x$achieved_cors$target_rho)) {
      cat(sprintf("  Target (rho)     : %.*f\n", digits, x$achieved_cors$target_rho))
    }
    cat(sprintf("  Achieved Pearson : %.*f\n", digits, x$achieved_cors$pearson_pooled))
    cat(sprintf("  Achieved Spearman: %.*f\n", digits, x$achieved_cors$spearman_pooled))
  }

  if (x$model == "3pl" && !is.null(x$guessing_summary)) {
    cat("\nGuessing (lower asymptote):\n")
    cat(sprintf("  Mean     : %.*f\n", digits, x$guessing_summary$mean))
    cat(sprintf("  SD       : %.*f\n", digits, x$guessing_summary$sd))
    cat(sprintf("  Range    : [%.*f, %.*f]\n",
                digits, x$guessing_summary$min,
                digits, x$guessing_summary$max))
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
    if (x$model %in% c("2pl", "3pl")) {
      oldpar <- par(mfrow = c(1, if (x$model == "3pl") 3 else 2))
      on.exit(par(oldpar))
      hist(x$data$beta, main = "Difficulty Distribution", xlab = "beta", col = "steelblue")
      plot(x$data$beta, x$data$lambda, main = "Difficulty vs Discrimination",
           xlab = "beta", ylab = "lambda", pch = 19, col = rgb(0, 0, 0.5, 0.5))
      abline(lm(lambda ~ beta, data = x$data), col = "red", lty = 2)
      if (x$model == "3pl") {
        hist(
          x$data$guessing,
          main = "Guessing Distribution",
          xlab = "guessing",
          col = "darkseagreen3"
        )
      }
    } else {
      hist(x$data$beta, main = "Difficulty Distribution", xlab = "beta", col = "steelblue")
    }
    return(invisible(NULL))
  }

  df <- x$data

  # Scatter plot
  p_scatter <- NULL
  if (x$model %in% c("2pl", "3pl")) {
    rho_achieved <- x$achieved$overall$cor_spearman_pooled
    target_rho <- .item_target_rho(x)
    subtitle <- if (identical(x$source, "custom")) {
      "Source: custom"
    } else if (identical(x$source, "hierarchical")) {
      sprintf("Source: hierarchical | Target rho: %.2f", target_rho)
    } else {
      sprintf("Method: %s | Target rho: %.2f", x$method, target_rho)
    }
    p_scatter <- ggplot2::ggplot(df, ggplot2::aes(x = beta, y = lambda)) +
      ggplot2::geom_point(alpha = 0.6, color = "darkblue", size = 2) +
      ggplot2::geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
      ggplot2::labs(
        title = sprintf("Difficulty vs Discrimination (Spearman r = %.3f)", rho_achieved),
        subtitle = subtitle,
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
  if (x$model %in% c("2pl", "3pl")) {
    p_lambda <- ggplot2::ggplot(df, ggplot2::aes(x = lambda)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 30,
                     fill = "coral", alpha = 0.5, color = "white") +
      ggplot2::geom_density(color = "darkred", linewidth = 1) +
      ggplot2::labs(title = "Discrimination Distribution", x = expression(lambda), y = "Density") +
      ggplot2::theme_minimal()
  }

  p_guessing <- NULL
  if (x$model == "3pl") {
    p_guessing <- ggplot2::ggplot(df, ggplot2::aes(x = guessing)) +
      ggplot2::geom_histogram(
        ggplot2::aes(y = ggplot2::after_stat(density)),
        bins = 30,
        fill = "darkseagreen3",
        alpha = 0.6,
        color = "white"
      ) +
      ggplot2::geom_density(color = "darkgreen", linewidth = 1) +
      ggplot2::labs(
        title = "Guessing Distribution",
        x = "guessing",
        y = "Density"
      ) +
      ggplot2::theme_minimal()
  }

  # Return based on type
  if (type == "scatter") {
    if (x$model %in% c("2pl", "3pl")) return(p_scatter) else return(p_beta)
  } else if (type == "density") {
    if (x$model %in% c("2pl", "3pl") &&
        requireNamespace("patchwork", quietly = TRUE)) {
      density_plots <- if (x$model == "3pl") {
        list(p_beta, p_lambda, p_guessing)
      } else {
        list(p_beta, p_lambda)
      }
      return(patchwork::wrap_plots(density_plots, ncol = length(density_plots)))
    } else {
      return(p_beta)
    }
  } else {  # both
    if (x$model %in% c("2pl", "3pl") &&
        requireNamespace("patchwork", quietly = TRUE)) {
      density_plots <- if (x$model == "3pl") {
        list(p_beta, p_lambda, p_guessing)
      } else {
        list(p_beta, p_lambda)
      }
      combined_top <- patchwork::wrap_plots(
        density_plots,
        ncol = length(density_plots)
      )
      return(patchwork::wrap_plots(combined_top, p_scatter, ncol = 1))
    } else {
      return(p_beta)
    }
  }
}


#' Extract item parameters as data frame
#' @param x An `item_params` object.
#' @param row.names NULL or a character vector giving the row names.
#' @param optional Logical. Compatibility argument for \code{as.data.frame()};
#'   currently ignored.
#' @param ... Additional arguments (ignored).
#' @return A data frame of item parameters.
#' @export
as.data.frame.item_params <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- x$data
  if (!is.null(row.names)) {
    if (!is.character(row.names) || length(row.names) != nrow(out) ||
        any(is.na(row.names))) {
      stop("`row.names` must be NULL or a character vector matching the number of item rows.")
    }
    row.names(out) <- row.names
  }
  out
}
