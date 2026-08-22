# Independent Phase 6 acceptance tests
#
# This file verifies the public SAC contract from outside the implementation.
# Keep fixtures deliberately small: topology, alignment, provenance, and exact
# identities matter here; large stochastic recovery studies belong to Phase 7.

.phase6_legacy_sac_formals <- c(
  "target_rho", "n_items", "model", "latent_shape", "item_source",
  "latent_params", "item_params", "reliability_metric", "c_init",
  "M_per_iter", "M_pre", "n_iter", "burn_in", "step_params",
  "c_bounds", "resample_items", "seed", "verbose"
)

.phase6_beta <- c(-1.1, -0.55, -0.1, 0.35, 0.8, 1.25)
.phase6_lambda <- c(0.75, 0.9, 1.05, 1.2, 1.35, 1.5)
.phase6_guessing <- c(0.08, 0.12, 0.16, 0.2, 0.24, 0.28)

.phase6_item_params <- function(model = c("2pl", "3pl"),
                                guessing = .phase6_guessing) {
  model <- match.arg(model)
  custom <- list(beta = .phase6_beta, lambda = .phase6_lambda)
  if (identical(model, "3pl")) custom$guessing <- guessing
  list(
    custom_params = custom,
    center_difficulties = FALSE
  )
}

.phase6_bank <- function(guessing = rep(0, length(.phase6_beta)),
                         model = "2pl",
                         beta = .phase6_beta,
                         lambda = .phase6_lambda) {
  .irtsimrel_item_bank(
    beta = beta,
    lambda_base = lambda,
    guessing = guessing,
    model = model
  )
}

.phase6_quiet <- function(expr) {
  suppressWarnings(suppressMessages(force(expr)))
}

.phase6_fixed_sac <- function(model = c("2pl", "3pl"),
                              guessing = rep(0, length(.phase6_beta)),
                              seed = 7601L,
                              ...) {
  model <- match.arg(model)
  args <- list(
    target_rho = 0.4,
    n_items = length(.phase6_beta),
    model = model,
    latent_shape = "normal",
    item_source = "custom",
    item_params = .phase6_item_params(model, guessing),
    reliability_metric = "info",
    c_init = 0.7,
    M_per_iter = 64L,
    M_pre = 128L,
    n_iter = 8L,
    burn_in = 4L,
    c_bounds = c(0.1, 3),
    resample_items = FALSE,
    seed = seed,
    verbose = FALSE,
    preflight_controls = list(
      M = 96L,
      split_check = FALSE,
      root_controls = list(min_grid = 65L, max_evals = 257L)
    ),
    evaluation_controls = list(
      n_forms = 2L,
      M = 128L,
      probs = c(0.1, 0.5, 0.9)
    )
  )
  .phase6_quiet(do.call(sac_calibrate, utils::modifyList(args, list(...))))
}


test_that("SAC appends Phase 6 controls after the exact legacy formal prefix", {
  observed <- names(formals(sac_calibrate))

  expect_identical(
    observed[seq_along(.phase6_legacy_sac_formals)],
    .phase6_legacy_sac_formals
  )
  expect_identical(
    observed[-seq_along(.phase6_legacy_sac_formals)],
    c("root_policy", "preflight_controls", "evaluation_controls")
  )
})


test_that("SAC preflight rejects unresolved and inadmissible branches by class", {
  theta <- seq(-2.5, 2.5, length.out = 101L)
  bank <- .phase6_bank()
  common <- list(
    theta = theta,
    theta_var = stats::var(theta),
    item_banks = list(bank),
    metric_internal = "info",
    c_bounds = c(0.1, 10),
    split_check = FALSE
  )

  run_with_curve <- function(curve, target, args = common) {
    local_mocked_bindings(
      .compute_rho_generic = function(c, ...) curve(c),
      .package = "IRTsimrel",
      .env = parent.frame()
    )
    do.call(
      .irtsimrel_sac_preflight,
      c(args, list(target_rho = target))
    )
  }

  expect_error(
    run_with_curve(
      function(c) 0.2 + 0.5 * stats::plogis(log(c)),
      target = 0.9
    ),
    class = "irtsimrel_branch_unavailable"
  )
  expect_error(
    run_with_curve(function(c) stats::plogis(-log(c)), target = 0.5),
    class = "irtsimrel_branch_unavailable"
  )
  expect_error(
    run_with_curve(
      function(c) stats::plogis(log(c)),
      target = stats::plogis(log(0.1))
    ),
    class = "irtsimrel_branch_unavailable"
  )

  split_banks <- list(
    .phase6_bank(beta = .phase6_beta - 2),
    .phase6_bank(beta = .phase6_beta + 2),
    .phase6_bank(beta = .phase6_beta - 2),
    .phase6_bank(beta = .phase6_beta + 2)
  )
  split_args <- common
  split_args$item_banks <- split_banks
  split_args$split_check <- TRUE
  split_args$split_log_tolerance <- 0.5

  local_mocked_bindings(
    .compute_rho_generic = function(c, beta_vec, ...) {
      if (beta_vec[[1L]] < -2) {
        stats::plogis(log(c / 0.5))
      } else {
        stats::plogis(log(c / 2))
      }
    },
    .package = "IRTsimrel"
  )
  expect_error(
    do.call(
      .irtsimrel_sac_preflight,
      c(split_args, list(target_rho = 0.5))
    ),
    class = "irtsimrel_topology_uncertain"
  )
})


test_that("3PL guessing zero has the exact seeded 2PL SAC numeric path", {
  twopl <- .phase6_fixed_sac("2pl", seed = 7611L)
  threepl_zero <- .phase6_fixed_sac(
    "3pl",
    guessing = rep(0, length(.phase6_beta)),
    seed = 7611L
  )

  exact_fields <- c(
    "theta_var", "c_init", "c_final", "c_star", "achieved_rho",
    "trajectory", "rho_trajectory", "rho_update_trajectory",
    "rho_scale_trajectory", "raw_trajectory", "step_size_trajectory",
    "gradient_trajectory", "projected", "projection_side"
  )
  for (field in exact_fields) {
    expect_identical(threepl_zero[[field]], twopl[[field]], info = field)
  }
  expect_identical(
    threepl_zero$achieved_distribution$form_reliabilities,
    twopl$achieved_distribution$form_reliabilities
  )
  expect_identical(threepl_zero$guessing_vec, rep(0, length(.phase6_beta)))
  expect_identical(twopl$guessing_vec, rep(0, length(.phase6_beta)))
})


test_that("fixed-form 3PL preserves every guessing parameter", {
  out <- .phase6_fixed_sac(
    "3pl",
    guessing = .phase6_guessing,
    seed = 7621L
  )

  expect_identical(out$item_scope, "fixed_form")
  expect_equal(out$guessing_vec, .phase6_guessing, tolerance = 0)
  expect_equal(out$items_base$data$guessing, .phase6_guessing, tolerance = 0)
  expect_equal(out$items_calib$data$guessing, .phase6_guessing, tolerance = 0)
  expect_equal(out$design_signature$guessing, .phase6_guessing, tolerance = 0)
  expect_equal(out$lambda_scaled, out$lambda_base * out$c_star,
               tolerance = 1e-14)
})


test_that("seeded SAC is reproducible and restores its caller RNG state", {
  set.seed(7631L)
  caller_state <- .Random.seed
  first <- .phase6_fixed_sac(
    "3pl",
    guessing = .phase6_guessing,
    seed = 7632L
  )
  expect_identical(.Random.seed, caller_state)

  second <- .phase6_fixed_sac(
    "3pl",
    guessing = .phase6_guessing,
    seed = 7632L
  )
  expect_identical(serialize(second, NULL), serialize(first, NULL))
  expect_true(first$rng_provenance$caller_state_restored)
  expect_true(first$rng_provenance$caller_state_restored_when_seeded)
  expect_false(first$rng_provenance$stream_overlap)
  expect_false(identical(
    first$rng_provenance$preflight_stream_seed,
    first$rng_provenance$evaluation_stream_seed
  ))
})


test_that("superpopulation achieved reliability is the exact holdout summary", {
  out <- .phase6_quiet(sac_calibrate(
    target_rho = 0.35,
    n_items = 8L,
    model = "3pl",
    latent_shape = "normal",
    item_source = "parametric",
    item_params = list(
      guessing_params = list(
        distribution = "uniform",
        min = 0.1,
        max = 0.3
      )
    ),
    reliability_metric = "info",
    c_init = 0.8,
    M_per_iter = 64L,
    M_pre = 128L,
    n_iter = 8L,
    burn_in = 4L,
    c_bounds = c(0.1, 3),
    resample_items = TRUE,
    seed = 7641L,
    preflight_controls = list(
      n_forms = 4L,
      M = 96L,
      split_check = TRUE,
      split_log_tolerance = 2,
      root_controls = list(min_grid = 65L, max_evals = 257L)
    ),
    evaluation_controls = list(
      n_forms = 5L,
      M = 128L,
      conf_level = 0.9,
      probs = c(0.1, 0.5, 0.9)
    )
  ))
  distribution <- out$achieved_distribution
  values <- distribution$form_reliability

  expect_identical(out$item_scope, "item_superpopulation")
  expect_length(values, 5L)
  expect_equal(out$achieved_rho, mean(values), tolerance = 0)
  expect_equal(distribution$sd, stats::sd(values), tolerance = 0)
  expect_equal(
    distribution$se,
    stats::sd(values) / sqrt(length(values)),
    tolerance = 0
  )
  expect_equal(
    unname(distribution$quantiles),
    unname(stats::quantile(values, probs = out$evaluation_controls$probs)),
    tolerance = 0
  )
  expect_identical(distribution$n_forms, 5L)
  expect_identical(distribution$M_per_form, 128L)
  expect_identical(distribution$aggregation, "mean_of_form_reliabilities")
  expect_false(out$rng_provenance$stream_overlap)
  expect_false(identical(
    out$rng_provenance$preflight_stream_seed,
    out$rng_provenance$evaluation_stream_seed
  ))
})


test_that("SAC update and post-update reliability trajectories are aligned", {
  reliability_curve <- function(c) stats::plogis(log(c))
  local_mocked_bindings(
    .compute_rho_generic = function(c, ...) reliability_curve(c),
    .package = "IRTsimrel"
  )
  out <- .phase6_fixed_sac("2pl", seed = 7651L)

  expect_equal(
    out$rho_update_trajectory,
    reliability_curve(out$rho_scale_trajectory),
    tolerance = 1e-14
  )
  expect_equal(
    out$rho_trajectory,
    reliability_curve(out$trajectory),
    tolerance = 1e-14
  )
  expect_equal(
    out$gradient_trajectory,
    out$rho_update_trajectory - out$target_rho,
    tolerance = 0
  )
  expect_equal(
    out$raw_trajectory,
    out$rho_scale_trajectory -
      out$step_size_trajectory * out$gradient_trajectory,
    tolerance = 1e-14
  )
  expect_equal(out$rho_scale_trajectory[[1L]], out$c_init, tolerance = 0)
  expect_equal(
    out$rho_scale_trajectory[-1L],
    out$trajectory[-length(out$trajectory)],
    tolerance = 0
  )
  expect_equal(
    out$trajectory,
    pmin(pmax(out$raw_trajectory, out$c_bounds[[1L]]), out$c_bounds[[2L]]),
    tolerance = 0
  )
  expect_true(all(out$trajectory >= out$c_bounds[[1L]]))
  expect_true(all(out$trajectory <= out$c_bounds[[2L]]))
  expect_true(out$c_bounds[[1L]] >= out$requested_c_bounds[[1L]])
  expect_true(out$c_bounds[[2L]] <= out$requested_c_bounds[[2L]])
})


test_that("SAC rejects nonintegrable population MSEM with a package class", {
  condition <- tryCatch(
    sac_calibrate(
      target_rho = 0.3,
      n_items = 6L,
      model = "2pl",
      latent_shape = "heavy_tail",
      latent_params = list(df = 3),
      item_source = "custom",
      item_params = .phase6_item_params("2pl"),
      reliability_metric = "msem",
      M_per_iter = 16L,
      M_pre = 32L,
      n_iter = 2L,
      burn_in = 1L,
      resample_items = FALSE,
      seed = 7661L,
      preflight_controls = list(M = 32L, split_check = FALSE),
      evaluation_controls = list(n_forms = 2L, M = 32L)
    ),
    error = identity
  )

  expect_s3_class(condition, "error")
  expect_true(any(grepl("^irtsimrel_", class(condition))))
  expect_match(conditionMessage(condition), "heavy.tail|integrab", ignore.case = TRUE)
})


test_that("public SAC propagates infeasible preflight as a classed error", {
  expect_error(
    sac_calibrate(
      target_rho = 0.999,
      n_items = 6L,
      model = "2pl",
      item_source = "custom",
      item_params = .phase6_item_params("2pl"),
      reliability_metric = "info",
      c_init = 0.15,
      M_per_iter = 16L,
      M_pre = 64L,
      n_iter = 2L,
      burn_in = 1L,
      c_bounds = c(0.1, 0.2),
      resample_items = FALSE,
      seed = 7671L,
      preflight_controls = list(
        M = 48L,
        split_check = FALSE,
        root_controls = list(min_grid = 33L, max_evals = 129L)
      ),
      evaluation_controls = list(n_forms = 2L, M = 32L)
    ),
    class = "irtsimrel_branch_unavailable"
  )
})


test_that("fixed EQC warm start reuses its full 3PL bank and is comparable", {
  eqc <- .phase6_quiet(eqc_calibrate(
    target_rho = 0.4,
    n_items = length(.phase6_beta),
    model = "3pl",
    latent_shape = "normal",
    item_source = "custom",
    item_params = .phase6_item_params("3pl", .phase6_guessing),
    reliability_metric = "info",
    M = 128L,
    c_bounds = c(0.1, 3),
    seed = 7681L,
    root_controls = list(min_grid = 65L, max_evals = 257L)
  ))
  # Omitting item_source/item_params is the public opt-in for reusing the
  # exact fixed EQC form rather than drawing a nominally similar form.
  sac <- .phase6_quiet(sac_calibrate(
    target_rho = 0.4,
    n_items = length(.phase6_beta),
    model = "3pl",
    latent_shape = "normal",
    reliability_metric = "info",
    c_init = eqc,
    M_per_iter = 64L,
    M_pre = 128L,
    n_iter = 8L,
    burn_in = 4L,
    c_bounds = c(0.1, 3),
    resample_items = FALSE,
    seed = 7682L,
    preflight_controls = list(
      M = 96L,
      split_check = FALSE,
      root_controls = list(min_grid = 65L, max_evals = 257L)
    ),
    evaluation_controls = list(
      n_forms = 2L,
      M = 128L,
      probs = c(0.1, 0.5, 0.9)
    )
  ))

  expect_identical(sac$init_method, "eqc_warm_start")
  expect_true(sac$warm_start$form_reused)
  expect_true(sac$warm_start$comparable)
  expect_equal(sac$beta_vec, eqc$beta_vec, tolerance = 0)
  expect_equal(sac$lambda_base, eqc$lambda_base, tolerance = 0)
  expect_equal(sac$guessing_vec, eqc$guessing_vec, tolerance = 0)

  comparison <- suppressWarnings(compare_eqc_sac(eqc, sac, verbose = FALSE))
  expect_true(comparison$comparable)
  expect_identical(comparison$comparability_reasons, character())

  expect_error(
    .phase6_quiet(sac_calibrate(
      target_rho = 0.4,
      n_items = length(.phase6_beta),
      model = "3pl",
      reliability_metric = "info",
      c_init = eqc,
      M_per_iter = 16L,
      M_pre = 48L,
      n_iter = 2L,
      burn_in = 1L,
      c_bounds = c(0.1, 3),
      resample_items = TRUE,
      seed = 7684L,
      preflight_controls = list(
        n_forms = 4L,
        M = 32L,
        split_check = FALSE
      ),
      evaluation_controls = list(n_forms = 2L, M = 32L)
    )),
    class = "irtsimrel_sac_warm_start_mismatch"
  )
  expect_error(
    .phase6_quiet(sac_calibrate(
      target_rho = 0.4,
      n_items = length(.phase6_beta),
      model = "3pl",
      reliability_metric = "msem",
      c_init = eqc,
      M_per_iter = 16L,
      M_pre = 48L,
      n_iter = 2L,
      burn_in = 1L,
      c_bounds = c(0.1, 3),
      resample_items = FALSE,
      seed = 7685L,
      preflight_controls = list(M = 32L, split_check = FALSE),
      evaluation_controls = list(n_forms = 2L, M = 32L)
    )),
    class = "irtsimrel_sac_warm_start_mismatch"
  )
  shifted_params <- .phase6_item_params("3pl", .phase6_guessing)
  shifted_params$custom_params$beta[[1L]] <-
    shifted_params$custom_params$beta[[1L]] + 0.05
  expect_error(
    .phase6_quiet(sac_calibrate(
      target_rho = 0.4,
      n_items = length(.phase6_beta),
      model = "3pl",
      item_source = "custom",
      item_params = shifted_params,
      reliability_metric = "info",
      c_init = eqc,
      M_per_iter = 16L,
      M_pre = 48L,
      n_iter = 2L,
      burn_in = 1L,
      c_bounds = c(0.1, 3),
      resample_items = FALSE,
      seed = 7686L,
      preflight_controls = list(M = 32L, split_check = FALSE),
      evaluation_controls = list(n_forms = 2L, M = 32L)
    )),
    class = "irtsimrel_sac_warm_start_mismatch"
  )

  mismatches <- list(
    metric = function(x) {
      x$metric <- "msem"
      x$estimand_signature$metric <- "msem"
      x$design_signature$metric <- "msem"
      x
    },
    scope = function(x) {
      x$item_scope <- "item_superpopulation"
      x$estimand_signature$item_measure <- "item_superpopulation"
      x$design_signature$item_scope <- "item_superpopulation"
      x$design_signature$item_source <- "custom"
      x$design_signature$item_params <- .phase6_item_params(
        "3pl", .phase6_guessing
      )
      x
    },
    bank = function(x) {
      x$design_signature$beta[[1L]] <-
        x$design_signature$beta[[1L]] + 0.01
      x
    },
    guessing = function(x) {
      x$design_signature$guessing[[1L]] <-
        x$design_signature$guessing[[1L]] + 0.01
      x
    },
    signature = function(x) {
      x$design_signature$latent_shape <- "uniform"
      x
    }
  )
  for (label in names(mismatches)) {
    changed <- mismatches[[label]](sac)
    result <- suppressWarnings(compare_eqc_sac(eqc, changed, verbose = FALSE))
    expect_false(result$comparable, info = label)
    expect_true(is.na(result$agreement), info = label)
  }

  model_mismatch <- .phase6_fixed_sac("2pl", seed = 7683L)
  result <- suppressWarnings(
    compare_eqc_sac(eqc, model_mismatch, verbose = FALSE)
  )
  expect_false(result$comparable)
  expect_true(is.na(result$agreement))

  nonconverged <- sac
  nonconverged$calibration_status <- "not_converged"
  nonconverged$convergence$status <- "not_converged"
  nonconverged$convergence$converged <- FALSE
  nonconverged$convergence$status_flags <- "not_converged"
  result <- suppressWarnings(
    compare_eqc_sac(eqc, nonconverged, verbose = FALSE)
  )
  expect_true(result$comparable)
  expect_true(is.na(result$agreement))
  expect_identical(result$agreement_status, "calibration_unsuccessful")
})


test_that("3PL SAC coef and predict preserve and use guessing", {
  out <- .phase6_fixed_sac(
    "3pl",
    guessing = .phase6_guessing,
    seed = 7691L
  )
  coefficients <- coef(out)
  legacy_prefix <- c(
    "item_id", "beta", "lambda_base", "lambda_scaled", "c_star"
  )

  expect_identical(names(coefficients)[seq_along(legacy_prefix)], legacy_prefix)
  expect_identical(names(coefficients)[[length(legacy_prefix) + 1L]], "guessing")
  expect_equal(coefficients$guessing, .phase6_guessing, tolerance = 0)

  theta <- seq(-2.5, 2.5, length.out = 151L)
  scales <- c(0.5, 0.9, 1.3)
  observed <- predict(out, newdata = scales, theta_vec = theta)
  expected <- vapply(
    scales,
    function(scale) {
      .compute_rho_generic(
        c = scale,
        theta_vec = theta,
        beta_vec = out$beta_vec,
        lambda_base = out$lambda_base,
        theta_var = stats::var(theta),
        metric_internal = out$metric,
        guessing = out$guessing_vec
      )
    },
    numeric(1)
  )
  expect_equal(unname(observed), expected, tolerance = 1e-13)
  expect_equal(
    unname(predict(out, newdata = out$c_star)),
    out$achieved_rho,
    tolerance = 1e-13
  )
})


test_that("legacy 2PL SAC is lazily zero-guessing and v0.1 loss is actionable", {
  legacy <- .phase6_fixed_sac("2pl", seed = 7701L)
  legacy$guessing_vec <- NULL
  legacy$schema_version <- NULL
  legacy$item_scope <- NULL
  legacy$estimand_signature <- NULL
  legacy$design_signature <- NULL
  if ("guessing" %in% names(legacy$items_base$data)) {
    legacy$items_base$data$guessing <- NULL
  }
  if ("guessing" %in% names(legacy$items_calib$data)) {
    legacy$items_calib$data$guessing <- NULL
  }
  before <- serialize(legacy, NULL)
  theta <- seq(-2, 2, length.out = 101L)
  observed <- predict(legacy, newdata = c(0.6, 1), theta_vec = theta)
  expected <- vapply(
    c(0.6, 1),
    function(scale) {
      .compute_rho_generic(
        c = scale,
        theta_vec = theta,
        beta_vec = legacy$beta_vec,
        lambda_base = legacy$lambda_base,
        theta_var = stats::var(theta),
        metric_internal = legacy$metric,
        guessing = rep(0, legacy$n_items)
      )
    },
    numeric(1)
  )
  expect_equal(unname(observed), expected, tolerance = 1e-13)
  expect_identical(serialize(legacy, NULL), before)

  literal_v01 <- structure(
    list(
      c_star = 1,
      c_final = 1,
      target_rho = 0.4,
      achieved_rho = 0.4,
      theta_var = 1,
      trajectory = rep(1, 4L),
      rho_trajectory = rep(0.4, 4L),
      metric = "info",
      model = "2pl",
      n_items = 6L,
      n_iter = 4L,
      burn_in = 2L,
      M_per_iter = 10L,
      M_pre = 20L,
      step_params = list(a = 1, A = 50, gamma = 0.67),
      c_bounds = c(0.1, 3),
      c_init = 1,
      init_method = "apc_warm_start",
      convergence = list(
        mean_first_half = 1,
        mean_second_half = 1,
        sd_post_burn = 0,
        range_post_burn = c(1, 1),
        final_gradient = 0,
        hit_lower_bound = FALSE,
        hit_upper_bound = FALSE,
        converged = TRUE
      ),
      call = quote(spc_calibrate(target_rho = 0.4, n_items = 6L))
    ),
    class = c("spc_result", "list")
  )
  literal_before <- serialize(literal_v01, NULL)
  expect_error(
    simulate_response_data(literal_v01, n_persons = 10L, seed = 7702L),
    "[Rr]e-?run.*calibration|[Rr]erun.*calibration"
  )
  expect_identical(serialize(literal_v01, NULL), literal_before)
})
