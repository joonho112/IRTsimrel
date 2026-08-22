make_v03_comparison_pair <- function() {
  beta <- c(-0.8, 0.1, 0.7)
  lambda <- c(0.9, 1.1, 1.3)
  guessing <- c(0.10, 0.20, 0.25)
  scale_convention <- "global_discrimination_multiplier_D1"

  eqc <- list(
    c_star = 1.05,
    target_rho = 0.80,
    achieved_rho = 0.80,
    metric = "info",
    model = "3pl",
    n_items = 3L,
    theta_var = 1,
    theta_quad = c(-1, 0, 1),
    beta_vec = beta,
    lambda_base = lambda,
    lambda_scaled = lambda * 1.05,
    guessing_vec = guessing,
    schema_version = 1L,
    item_scope = "fixed_form",
    estimand_signature = list(
      metric = "info",
      theta_measure = "population",
      approximation_method = "fixed_mc",
      item_measure = "fixed_form"
    ),
    design_signature = list(
      schema_version = 1L,
      model = "3pl",
      metric = "info",
      n_items = 3L,
      item_scope = "fixed_form",
      latent_shape = "normal",
      latent_params = list(),
      scale_convention = scale_convention,
      beta = beta,
      lambda_base = lambda,
      guessing = guessing
    ),
    misc = list(root_status = "uniroot_success"),
    call = quote(eqc_calibrate(target_rho = 0.8, n_items = 3L))
  )
  class(eqc) <- c("eqc_result", "list")

  sac <- eqc
  sac$c_star <- 1.06
  sac$lambda_scaled <- lambda * 1.06
  sac$estimand_signature$approximation_method <- "stochastic_mc"
  sac$calibration_status <- "ok"
  sac$convergence <- list(status = "ok", status_flags = "ok")
  sac$misc <- NULL
  sac$call <- quote(sac_calibrate(
    target_rho = 0.8,
    n_items = 3L,
    resample_items = FALSE
  ))
  class(sac) <- c("sac_result", "list")

  list(eqc = eqc, sac = sac)
}

quiet_compare_v03 <- function(eqc, sac) {
  suppressWarnings(compare_eqc_sac(eqc, sac, verbose = FALSE))
}

test_that("same fixed 3PL bank is comparable and preserves legacy prefix", {
  pair <- make_v03_comparison_pair()
  out <- quiet_compare_v03(pair$eqc, pair$sac)

  legacy_prefix <- c(
    "c_eqc", "c_sac", "diff_abs", "diff_pct", "agreement",
    "target_rho", "achieved_eqc", "achieved_sac", "achieved_diff_abs",
    "metric_eqc", "metric_sac", "model_eqc", "model_sac",
    "n_items_eqc", "n_items_sac", "eqc_status", "sac_status",
    "sac_status_flags"
  )
  expect_identical(names(out)[seq_along(legacy_prefix)], legacy_prefix)
  expect_true(out$comparable)
  expect_identical(out$comparability_reasons, character())
  expect_identical(out$agreement_status, "evaluated")
  expect_true(out$agreement)
})


test_that("fixed-form comparability ignores irrelevant generator functions", {
  item_params <- list(custom_params = list(
    beta = function(n) seq(-1, 1, length.out = n),
    lambda = function(n) seq(0.8, 1.4, length.out = n)
  ))
  eqc <- suppressMessages(eqc_calibrate(
    target_rho = 0.4,
    n_items = 4L,
    model = "2pl",
    latent_shape = "normal",
    item_source = "custom",
    item_params = item_params,
    M = 128L,
    c_bounds = c(0.1, 3),
    seed = 9911L,
    root_controls = list(min_grid = 33L, max_evals = 129L)
  ))
  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.4,
    n_items = 4L,
    model = "2pl",
    latent_shape = "normal",
    item_source = "custom",
    item_params = item_params,
    reliability_metric = "info",
    c_init = eqc,
    M_per_iter = 32L,
    M_pre = 64L,
    n_iter = 4L,
    burn_in = 2L,
    c_bounds = c(0.1, 3),
    resample_items = FALSE,
    seed = 9912L,
    preflight_controls = list(
      M = 48L,
      split_check = FALSE,
      root_controls = list(min_grid = 33L, max_evals = 129L)
    ),
    evaluation_controls = list(n_forms = 2L, M = 48L)
  )))

  expect_true(sac$design_signature$verifiable)
  expect_null(sac$design_signature$item_params)
  comparison <- suppressWarnings(compare_eqc_sac(eqc, sac, verbose = FALSE))
  expect_true(comparison$comparable)
  expect_identical(comparison$comparability_reasons, character())
})


test_that("comparison canonicalizes reconstructed legacy shape names", {
  pair <- make_v03_comparison_pair()
  pair$eqc$design_signature$latent_shape <- "norm"

  comparison <- quiet_compare_v03(pair$eqc, pair$sac)
  expect_true(comparison$comparable)
  expect_false("latent_specification_mismatch" %in%
                 comparison$comparability_reasons)
})

test_that("target mismatch remains a hard comparison error", {
  pair <- make_v03_comparison_pair()
  pair$sac$target_rho <- 0.81
  expect_error(
    compare_eqc_sac(pair$eqc, pair$sac, verbose = FALSE),
    "target_rho differs"
  )
})

test_that("canonical estimand and design mismatches withhold agreement", {
  mutators <- list(
    model_mismatch = function(x) {
      x$model <- "2pl"
      x$design_signature$model <- "2pl"
      x
    },
    metric_mismatch = function(x) {
      x$metric <- "msem"
      x$estimand_signature$metric <- "msem"
      x$design_signature$metric <- "msem"
      x
    },
    n_items_mismatch = function(x) {
      x$n_items <- 4L
      x$design_signature$n_items <- 4L
      x
    },
    scale_convention_mismatch = function(x) {
      x$design_signature$scale_convention <- "difficulty_multiplier_D1"
      x
    },
    theta_measure_mismatch = function(x) {
      x$estimand_signature$theta_measure <- "sample"
      x
    },
    latent_specification_mismatch = function(x) {
      x$design_signature$latent_shape <- "uniform"
      x
    },
    item_scope_mismatch = function(x) {
      x$item_scope <- "item_superpopulation"
      x$estimand_signature$item_measure <- "item_superpopulation"
      x$design_signature$item_scope <- "item_superpopulation"
      x$design_signature$item_source <- "parametric"
      x$design_signature$item_params <- list()
      x
    },
    fixed_bank_beta_mismatch = function(x) {
      x$design_signature$beta[1L] <- x$design_signature$beta[1L] + 0.01
      x
    },
    fixed_bank_lambda_mismatch = function(x) {
      x$design_signature$lambda_base[1L] <-
        x$design_signature$lambda_base[1L] + 0.01
      x
    },
    fixed_bank_guessing_mismatch = function(x) {
      x$design_signature$guessing[1L] <-
        x$design_signature$guessing[1L] + 0.01
      x
    }
  )

  for (reason in names(mutators)) {
    pair <- make_v03_comparison_pair()
    pair$sac <- mutators[[reason]](pair$sac)
    out <- quiet_compare_v03(pair$eqc, pair$sac)
    expect_false(out$comparable, info = reason)
    expect_true(reason %in% out$comparability_reasons, info = reason)
    expect_true(is.na(out$agreement), info = reason)
    expect_identical(out$agreement_status, "not_comparable", info = reason)
  }
})

test_that("internally inconsistent and unstable signatures are diagnosed", {
  pair <- make_v03_comparison_pair()
  pair$sac$design_signature$model <- "2pl"
  inconsistent <- quiet_compare_v03(pair$eqc, pair$sac)
  expect_false(inconsistent$comparable)
  expect_true(
    "sac_design_signature_inconsistent" %in%
      inconsistent$comparability_reasons
  )

  pair <- make_v03_comparison_pair()
  for (object_name in c("eqc", "sac")) {
    x <- pair[[object_name]]
    x$item_scope <- "item_superpopulation"
    x$estimand_signature$item_measure <- "item_superpopulation"
    x$design_signature$item_scope <- "item_superpopulation"
    x$design_signature$item_source <- "custom"
    x$design_signature$item_params <- list(
      custom_params = list(beta = function(n) stats::rnorm(n))
    )
    pair[[object_name]] <- x
  }
  unstable <- quiet_compare_v03(pair$eqc, pair$sac)
  expect_false(unstable$comparable)
  expect_true("generator_signature_unstable" %in%
                unstable$comparability_reasons)
})

test_that("superpopulation generator specifications compare canonically", {
  pair <- make_v03_comparison_pair()
  for (object_name in c("eqc", "sac")) {
    x <- pair[[object_name]]
    x$item_scope <- "item_superpopulation"
    x$estimand_signature$item_measure <- "item_superpopulation"
    x$design_signature$item_scope <- "item_superpopulation"
    x$design_signature$item_source <- "parametric"
    x$design_signature$item_params <- list(
      guessing_params = list(value = 0.2, distribution = "fixed"),
      difficulty_params = list(sigma = 1, mu = 0)
    )
    pair[[object_name]] <- x
  }
  # Named-list order is cosmetic.
  pair$sac$design_signature$item_params <-
    pair$sac$design_signature$item_params[2:1]
  out <- quiet_compare_v03(pair$eqc, pair$sac)
  expect_true(out$comparable)

  pair$sac$design_signature$item_params$difficulty_params$sigma <- 1.1
  mismatch <- quiet_compare_v03(pair$eqc, pair$sac)
  expect_false(mismatch$comparable)
  expect_true("generator_signature_mismatch" %in%
                mismatch$comparability_reasons)
})

test_that("non-successful calibration keeps comparability but withholds agreement", {
  pair <- make_v03_comparison_pair()
  pair$sac$calibration_status <- "branch_lost"
  pair$sac$convergence$status <- "branch_lost"
  pair$sac$convergence$status_flags <- "branch_lost"
  out <- quiet_compare_v03(pair$eqc, pair$sac)

  expect_true(out$comparable)
  expect_true(is.na(out$agreement))
  expect_identical(out$agreement_status, "calibration_unsuccessful")
  expect_identical(out$agreement_reasons, "sac_calibration_unsuccessful")
})

test_that("safe legacy fixed-form provenance is inferred", {
  pair <- make_v03_comparison_pair()
  legacy <- pair$sac
  legacy$schema_version <- NULL
  legacy$item_scope <- NULL
  legacy$estimand_signature <- NULL
  legacy$design_signature <- NULL
  legacy$guessing_vec <- NULL
  legacy$model <- "2pl"
  legacy$item_design <- "fixed_iteration_items"
  legacy$call <- quote(sac_calibrate(
    target_rho = 0.8,
    n_items = 3L,
    model = "2pl",
    reliability_metric = "info",
    resample_items = FALSE
  ))
  legacy$lambda_scaled <- legacy$lambda_base * legacy$c_star

  eqc <- pair$eqc
  eqc$model <- "2pl"
  eqc$guessing_vec <- NULL
  eqc$design_signature$model <- "2pl"
  eqc$design_signature$guessing <- rep(0, 3L)

  out <- quiet_compare_v03(eqc, legacy)
  expect_true(out$comparable)
  expect_identical(out$comparability_reasons, character())
})

test_that("literal designless v0.1 SPC is returned as not comparable", {
  pair <- make_v03_comparison_pair()
  legacy <- pair$sac[c(
    "c_star", "target_rho", "achieved_rho", "metric", "model", "n_items",
    "calibration_status", "convergence"
  )]
  class(legacy) <- c("spc_result", "list")

  out <- quiet_compare_v03(pair$eqc, legacy)
  expect_false(out$comparable)
  expect_true("sac_design_signature_missing" %in%
                out$comparability_reasons)
  expect_true(is.na(out$agreement))
})

make_v03_sac_validator_fixture <- function(model = "3pl") {
  n_iter <- 4L
  beta <- c(-0.5, 0.5)
  lambda <- c(0.8, 1.2)
  guessing <- if (identical(model, "3pl")) c(0.15, 0.25) else c(0, 0)
  c_star <- 1.1
  item_data <- data.frame(
    form_id = 1L,
    item_id = 1:2,
    beta = beta,
    lambda = lambda,
    lambda_unscaled = lambda,
    guessing = guessing
  )
  calib_data <- item_data
  calib_data$lambda <- lambda * c_star

  out <- list(
    c_star = c_star,
    c_final = c_star,
    target_rho = 0.7,
    achieved_rho = 0.7,
    theta_var = 1,
    trajectory = rep(c_star, n_iter),
    rho_trajectory = rep(0.7, n_iter),
    rho_scale_trajectory = rep(c_star, n_iter),
    rho_update_trajectory = rep(0.7, n_iter),
    evaluation_trajectory = rep(c_star, n_iter),
    iteration_trace = data.frame(
      iteration = seq_len(n_iter),
      c = rep(c_star, n_iter),
      rho = rep(0.7, n_iter)
    ),
    raw_trajectory = rep(c_star, n_iter),
    step_size_trajectory = rep(0.1, n_iter),
    gradient_trajectory = rep(0, n_iter),
    projected = rep(FALSE, n_iter),
    projection_side = rep("none", n_iter),
    projection_count = 0L,
    projection_rate = 0,
    metric = "info",
    calibration_status = "branch_lost",
    model = model,
    n_items = 2L,
    n_iter = n_iter,
    burn_in = 2L,
    M_per_iter = 10L,
    M_pre = 20L,
    step_params = list(a = 1, A = 50, gamma = 0.67),
    c_bounds = c(0.1, 2),
    c_init = 1,
    init_method = "numeric",
    convergence = list(
      converged = FALSE,
      sd_post_burn = 0,
      hit_lower_bound = FALSE,
      hit_upper_bound = FALSE,
      achieved_gap_abs = 0,
      achieved_gap_tolerance = 0.05,
      large_achieved_gap = FALSE,
      status = "branch_lost",
      status_flags = "branch_lost"
    ),
    beta_vec = beta,
    lambda_base = lambda,
    lambda_scaled = lambda * c_star,
    guessing_vec = guessing,
    items_base = list(data = item_data, source = "custom"),
    items_calib = list(data = calib_data, source = "custom"),
    theta_quad = c(-1, 0, 1),
    call = quote(sac_calibrate(target_rho = 0.7, n_items = 2L)),
    schema_version = 1L,
    item_scope = "fixed_form",
    estimand_signature = list(
      metric = "info",
      theta_measure = "population",
      approximation_method = "stochastic_mc",
      item_measure = "fixed_form"
    ),
    design_signature = list(
      schema_version = 1L,
      model = model,
      metric = "info",
      n_items = 2L,
      item_scope = "fixed_form",
      latent_shape = "normal",
      latent_params = list(),
      scale_convention = "global_discrimination_multiplier_D1",
      beta = beta,
      lambda_base = lambda,
      guessing = guessing
    )
  )
  class(out) <- c("sac_result", "list")
  out
}

test_that("v0.3 SAC validator accepts 3PL branch and trace additions", {
  object <- make_v03_sac_validator_fixture()
  expect_invisible(
    .irtsimrel_validate_sac_result_object(object, "object", "test")
  )

  bad <- object
  bad$rho_update_trajectory <- bad$rho_update_trajectory[-1L]
  expect_error(
    .irtsimrel_validate_sac_result_object(bad, "object", "test"),
    "rho_update_trajectory"
  )

  bad <- object
  bad$iteration_trace <- bad$iteration_trace[-1L, ]
  expect_error(
    .irtsimrel_validate_sac_result_object(bad, "object", "test"),
    "iteration_trace"
  )
})

test_that("legacy 2PL SAC validator lazily supplies zero guessing", {
  object <- make_v03_sac_validator_fixture("2pl")
  object$schema_version <- NULL
  object$item_scope <- NULL
  object$estimand_signature <- NULL
  object$design_signature <- NULL
  object$guessing_vec <- NULL
  object$items_base$data$guessing <- NULL
  object$items_calib$data$guessing <- NULL
  object$calibration_status <- "ok"
  object$convergence$status <- "ok"
  object$convergence$status_flags <- "ok"
  object$convergence$converged <- TRUE

  expect_invisible(
    .irtsimrel_validate_sac_result_object(object, "object", "test")
  )
})
