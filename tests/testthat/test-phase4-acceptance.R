# Independent Phase 4 acceptance tests
#
# These tests lock the public item/response contract without relying on the
# calibration algorithms that are introduced in later phases.

.phase4_item_data <- function(model,
                              beta,
                              lambda,
                              guessing = NULL) {
  out <- data.frame(
    form_id = 1L,
    item_id = seq_along(beta),
    beta = beta,
    lambda = lambda,
    lambda_unscaled = lambda
  )
  if (identical(model, "3pl")) {
    out$guessing <- guessing
  }
  out
}

.phase4_result_fixture <- function(model = "3pl",
                                   beta = 0,
                                   lambda = 1,
                                   guessing = 0.2,
                                   include_top_level_guessing = TRUE) {
  data <- .phase4_item_data(model, beta, lambda, guessing)
  items <- structure(
    list(
      data = data,
      model = model,
      source = "custom",
      method = NA_character_,
      n_items = as.integer(length(beta)),
      n_forms = 1L,
      scale = 1,
      centered = FALSE,
      params = list(
        difficulty = list(),
        discrimination = list(),
        hierarchical = list(),
        custom = list(beta = beta, lambda = lambda, guessing = guessing),
        guessing = list(distribution = "custom")
      ),
      achieved = list()
    ),
    class = c("item_params", "list")
  )

  out <- list(
    c_star = 1,
    target_rho = 0.5,
    achieved_rho = 0.5,
    metric = "info",
    model = model,
    n_items = as.integer(length(beta)),
    M = 100L,
    theta_quad = c(-1, 1),
    theta_var = 1,
    beta_vec = beta,
    lambda_base = lambda,
    lambda_scaled = lambda,
    items_base = items,
    items_calib = items,
    call = quote(eqc_calibrate(target_rho = 0.5, n_items = 1)),
    misc = list(root_status = "ok"),
    schema_version = 1L,
    item_scope = "fixed_form",
    estimand_signature = list(metric = "info", theta_measure = "population"),
    design_signature = list(
      version = 1L,
      model = model,
      beta = beta,
      lambda_base = lambda,
      guessing = guessing
    )
  )
  if (include_top_level_guessing) {
    out$guessing_vec <- if (identical(model, "3pl")) {
      guessing
    } else {
      rep(0, length(beta))
    }
  }
  structure(out, class = c("eqc_result", "list"))
}

.phase4_literal_v01_spc <- function() {
  convergence <- list(
    mean_first_half = 0.9,
    mean_second_half = 0.9,
    sd_post_burn = 0.01,
    range_post_burn = c(0.88, 0.92),
    final_gradient = 0,
    hit_lower_bound = FALSE,
    hit_upper_bound = FALSE,
    converged = TRUE
  )

  structure(
    list(
      c_star = 0.9,
      c_final = 0.9,
      target_rho = 0.8,
      achieved_rho = 0.8,
      theta_var = 1,
      trajectory = rep(0.9, 20L),
      rho_trajectory = rep(0.8, 20L),
      metric = "info",
      model = "rasch",
      n_items = 10L,
      n_iter = 20L,
      burn_in = 10L,
      M_per_iter = 100L,
      M_pre = 1000L,
      step_params = list(a = 1, A = 50, gamma = 0.67),
      c_bounds = c(0.1, 5),
      c_init = 0.9,
      init_method = "apc_warm_start",
      convergence = convergence,
      call = quote(spc_calibrate(target_rho = 0.8, n_items = 10))
    ),
    class = c("spc_result", "list")
  )
}


test_that("sim_item_params appends guessing_params after the v0.2 formals", {
  v02_formals <- c(
    "n_items", "model", "source", "method", "n_forms",
    "difficulty_params", "discrimination_params", "hierarchical_params",
    "custom_params", "scale", "center_difficulties", "seed"
  )
  observed <- names(formals(sim_item_params))

  expect_identical(observed[seq_along(v02_formals)], v02_formals)
  expect_identical(observed[length(observed)], "guessing_params")
  expect_length(observed, length(v02_formals) + 1L)
})


test_that("Rasch and 2PL retain their exact five-column item schema", {
  expected <- c(
    "form_id", "item_id", "beta", "lambda", "lambda_unscaled"
  )
  rasch <- sim_item_params(5, model = "rasch", seed = 4101)
  twopl <- sim_item_params(5, model = "2pl", seed = 4102)

  expect_identical(names(rasch$data), expected)
  expect_identical(names(twopl$data), expected)
  expect_false("guessing" %in% names(rasch$data))
  expect_false("guessing" %in% names(twopl$data))
})


test_that("3PL fixed, beta, and uniform guessing generators follow the contract", {
  expected_columns <- c(
    "form_id", "item_id", "beta", "lambda", "lambda_unscaled", "guessing"
  )

  default_fixed <- sim_item_params(7, model = "3pl", seed = 4201)
  fixed_vector <- sim_item_params(
    4,
    model = "3pl",
    guessing_params = list(
      distribution = "fixed",
      value = c(0.05, 0.10, 0.15, 0.20)
    ),
    seed = 4202
  )
  beta_generated <- suppressWarnings(sim_item_params(
    500,
    model = "3pl",
    guessing_params = list(distribution = "beta", shape1 = 5, shape2 = 17),
    seed = 4203
  ))
  uniform_generated <- sim_item_params(
    500,
    model = "3pl",
    guessing_params = list(distribution = "uniform", min = 0.1, max = 0.3),
    seed = 4204
  )

  expect_identical(names(default_fixed$data), expected_columns)
  expect_equal(default_fixed$data$guessing, rep(0.2, 7))
  expect_equal(
    fixed_vector$data$guessing,
    c(0.05, 0.10, 0.15, 0.20)
  )
  expect_true(all(is.finite(beta_generated$data$guessing)))
  expect_true(all(beta_generated$data$guessing > 0))
  expect_true(all(beta_generated$data$guessing < 1))
  expect_equal(mean(beta_generated$data$guessing), 5 / 22, tolerance = 0.035)
  expect_true(all(uniform_generated$data$guessing >= 0.1))
  expect_true(all(uniform_generated$data$guessing <= 0.3))
  expect_equal(
    mean(uniform_generated$data$guessing),
    0.2,
    tolerance = 0.02
  )

  expect_identical(default_fixed$params$guessing$distribution, "fixed")
  expect_identical(beta_generated$params$guessing$distribution, "beta")
  expect_identical(uniform_generated$params$guessing$distribution, "uniform")
})


test_that("guessing generators reject malformed or irrelevant specifications", {
  make_3pl <- function(spec) {
    sim_item_params(4, model = "3pl", guessing_params = spec)
  }

  expect_error(
    make_3pl(list(distribution = "fixed", value = -0.01)),
    "guessing|0.*1"
  )
  expect_error(
    make_3pl(list(distribution = "fixed", value = 1)),
    "guessing|0.*1"
  )
  expect_error(
    make_3pl(list(distribution = "fixed", value = c(0.1, 0.2))),
    "length|n_items"
  )
  expect_error(
    make_3pl(list(distribution = "beta", shape1 = 0, shape2 = 2)),
    "shape1|positive"
  )
  expect_error(
    make_3pl(list(distribution = "beta", shape1 = 2, shape2 = Inf)),
    "shape2|finite"
  )
  expect_error(
    make_3pl(list(distribution = "uniform", min = 0.2, max = 0.2)),
    "min|max"
  )
  expect_error(
    make_3pl(list(distribution = "uniform", min = -0.1, max = 0.2)),
    "min|0"
  )
  expect_error(
    make_3pl(list(distribution = "uniform", min = 0.1, max = 1)),
    "max|1"
  )
  expect_error(
    make_3pl(list(distribution = "logit-normal")),
    "distribution|fixed|beta|uniform"
  )
  expect_error(
    make_3pl(list(distribution = "fixed", value = 0.2, shape1 = 2)),
    "shape1|unknown|unused|irrelevant"
  )

  for (legacy_model in c("rasch", "2pl")) {
    expect_error(
      sim_item_params(
        4,
        model = legacy_model,
        guessing_params = list(distribution = "fixed", value = 0.2)
      ),
      "guessing_params|3pl"
    )
  }
})


test_that("high guessing warns and seeded guessing generation restores RNG", {
  expect_warning(
    sim_item_params(
      4,
      model = "3pl",
      guessing_params = list(distribution = "fixed", value = 0.5)
    ),
    "0\\.5|practical|high"
  )

  set.seed(4310)
  seed_before <- .Random.seed
  invisible(sim_item_params(
    20,
    model = "3pl",
    guessing_params = list(distribution = "beta", shape1 = 4, shape2 = 16),
    seed = 4311
  ))
  expect_identical(.Random.seed, seed_before)
})


test_that("custom 3PL parameters are retained as custom provenance", {
  beta <- c(-1, -0.25, 0.5, 1.25)
  lambda <- c(0.8, 1.0, 1.2, 1.4)
  guessing <- c(0.05, 0.10, 0.15, 0.20)
  items <- sim_item_params(
    4,
    model = "3pl",
    source = "custom",
    custom_params = list(
      beta = beta,
      lambda = lambda,
      guessing = guessing
    ),
    center_difficulties = FALSE
  )

  expect_true(is.na(items$method))
  expect_equal(items$data$beta, beta)
  expect_equal(items$data$lambda_unscaled, lambda)
  expect_equal(items$data$guessing, guessing)
  expect_true("custom" %in% names(items$params))
  expect_equal(items$params$custom$beta, beta)
  expect_equal(items$params$custom$lambda, lambda)
  expect_equal(items$params$custom$guessing, guessing)

  generated <- sim_item_params(
    3,
    model = "3pl",
    source = "custom",
    custom_params = list(
      beta = function(n) seq(-1, 1, length.out = n),
      lambda = function(n) seq(0.8, 1.2, length.out = n),
      guessing = function(n) seq(0.1, 0.2, length.out = n)
    ),
    center_difficulties = FALSE
  )
  expect_equal(generated$data$guessing, c(0.1, 0.15, 0.2))

  scalar <- sim_item_params(
    3,
    model = "3pl",
    source = "custom",
    custom_params = list(beta = -1:1, lambda = rep(1, 3), guessing = 0.15),
    center_difficulties = FALSE
  )
  expect_equal(scalar$data$guessing, rep(0.15, 3))

  expect_error(
    sim_item_params(
      3,
      model = "3pl",
      source = "custom",
      custom_params = list(beta = -1:1, lambda = rep(1, 3), guessing = 0.2),
      guessing_params = list(distribution = "fixed", value = 0.2)
    ),
    "ambig|guessing_params|only one"
  )
  expect_error(
    sim_item_params(
      3,
      model = "3pl",
      source = "custom",
      custom_params = list(beta = -1:1, lambda = rep(1, 3))
    ),
    "custom_params.*guessing|guessing.*required"
  )
  expect_error(
    sim_item_params(
      3,
      model = "2pl",
      source = "custom",
      custom_params = list(
        beta = -1:1,
        lambda = rep(1, 3),
        guessing = 0.2
      )
    ),
    "guessing|3pl"
  )
})


test_that("hierarchical sources handle n_items one and preserve Rasch slopes", {
  twopl_one <- sim_item_params(
    1,
    model = "2pl",
    source = "hierarchical",
    seed = 4401
  )
  threepl_one <- sim_item_params(
    1,
    model = "3pl",
    source = "hierarchical",
    seed = 4402
  )
  rasch <- sim_item_params(
    8,
    model = "rasch",
    source = "hierarchical",
    scale = 2.5,
    seed = 4403
  )

  expect_equal(nrow(twopl_one$data), 1L)
  expect_true(is.finite(twopl_one$data$beta))
  expect_gt(twopl_one$data$lambda_unscaled, 0)
  expect_equal(nrow(threepl_one$data), 1L)
  expect_equal(threepl_one$data$guessing, 0.2)
  expect_equal(rasch$data$lambda_unscaled, rep(1, 8))
  expect_equal(rasch$data$lambda, rep(2.5, 8))
})


test_that("global scale changes discrimination but not generated guessing", {
  args <- list(
    n_items = 40,
    model = "3pl",
    source = "parametric",
    guessing_params = list(distribution = "beta", shape1 = 4, shape2 = 16),
    seed = 4501
  )
  base <- do.call(sim_item_params, c(args, list(scale = 1)))
  scaled <- do.call(sim_item_params, c(args, list(scale = 3)))

  expect_identical(scaled$data$beta, base$data$beta)
  expect_identical(scaled$data$lambda_unscaled, base$data$lambda_unscaled)
  expect_identical(scaled$data$guessing, base$data$guessing)
  expect_equal(scaled$data$lambda, 3 * base$data$lambda)

  adapted <- IRTsimrel:::.irtsimrel_apply_item_scale(
    base,
    base$data$lambda_unscaled,
    scale = 2
  )
  expect_identical(adapted$data$guessing, base$data$guessing)
})


test_that("3PL item S3 methods expose discrimination and guessing", {
  items <- sim_item_params(
    40,
    model = "3pl",
    guessing_params = list(distribution = "uniform", min = 0.1, max = 0.3),
    seed = 4601
  )

  printed <- capture.output(returned <- print(items))
  expect_identical(returned, items)
  expect_true(any(grepl("3PL", printed, fixed = TRUE)))
  expect_true(any(grepl("Guess", printed, ignore.case = TRUE)))

  summary_items <- summary(items)
  expect_s3_class(summary_items, "summary.item_params")
  expect_true("lambda_summary" %in% names(summary_items))
  expect_true("guessing_summary" %in% names(summary_items))
  expect_true(is.list(summary_items$lambda_summary))
  expect_true(is.list(summary_items$guessing_summary))
  summary_printed <- capture.output(summary_returned <- print(summary_items))
  expect_identical(summary_returned, summary_items)
  expect_true(any(grepl("Guess", summary_printed, ignore.case = TRUE)))

  expect_identical(
    names(as.data.frame(items)),
    c("form_id", "item_id", "beta", "lambda", "lambda_unscaled", "guessing")
  )

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(plot(items, type = "scatter"), "gg")
    density_plot <- plot(items, type = "density")
    expect_true(inherits(density_plot, "gg") || inherits(density_plot, "patchwork"))
  } else {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_no_error(plot(items, type = "scatter"))
    expect_no_error(plot(items, type = "density"))
  }
})


test_that("3PL response simulation agrees with its empirical success probability", {
  result <- .phase4_result_fixture(
    model = "3pl",
    beta = 0,
    lambda = 1.3,
    guessing = 0.25
  )
  simulated <- simulate_response_data(
    result,
    n_persons = 50000,
    latent_shape = "normal",
    latent_params = list(mu = 0.35, sigma = 0.8),
    seed = 4701
  )

  expected_probability <- mean(
    0.25 + 0.75 * plogis(1.3 * (simulated$theta - 0))
  )
  expect_equal(
    mean(simulated$response_matrix),
    expected_probability,
    tolerance = 0.008
  )
  expect_equal(simulated$guessing, 0.25)
  expect_identical(
    names(simulated),
    c("response_matrix", "theta", "beta", "lambda", "provenance", "guessing")
  )
  expect_true("item_parameter_provenance" %in% names(simulated$provenance))
  expect_equal(simulated$provenance$model, "3pl")
})


test_that("g zero gives seeded response equivalence with the 2PL", {
  beta <- c(-0.75, 0.25, 1)
  lambda <- c(0.8, 1.1, 1.4)
  twopl <- .phase4_result_fixture(
    model = "2pl",
    beta = beta,
    lambda = lambda,
    guessing = rep(0, 3)
  )
  threepl_zero <- .phase4_result_fixture(
    model = "3pl",
    beta = beta,
    lambda = lambda,
    guessing = rep(0, 3)
  )

  simulated_2pl <- simulate_response_data(twopl, 1000, seed = 4801)
  simulated_3pl <- simulate_response_data(threepl_zero, 1000, seed = 4801)

  expect_identical(simulated_3pl$theta, simulated_2pl$theta)
  expect_identical(
    simulated_3pl$response_matrix,
    simulated_2pl$response_matrix
  )
  expect_identical(simulated_3pl$guessing, simulated_2pl$guessing)
})


test_that("legacy Rasch and 2PL results lazily receive zero guessing", {
  legacy <- .phase4_result_fixture(
    model = "2pl",
    beta = c(-0.5, 0.5),
    lambda = c(0.9, 1.1),
    guessing = c(0, 0),
    include_top_level_guessing = FALSE
  )
  legacy$items_base$params$custom$guessing <- NULL
  legacy$items_base$params$guessing <- NULL
  legacy$items_calib$params$custom$guessing <- NULL
  legacy$items_calib$params$guessing <- NULL
  legacy$items_base$params$custom <- NULL
  legacy$items_calib$params$custom <- NULL
  legacy$schema_version <- NULL
  legacy$item_scope <- NULL
  legacy$estimand_signature <- NULL
  legacy$design_signature <- NULL
  before <- serialize(legacy, NULL)

  simulated <- simulate_response_data(legacy, 100, seed = 4901)

  expect_identical(simulated$guessing, c(0, 0))
  expect_identical(serialize(legacy, NULL), before)
  expect_false("guessing_vec" %in% names(legacy))
  expect_false("guessing" %in% names(legacy$items_calib$data))
})


test_that("response validation detects beta and guessing provenance mismatches", {
  result <- .phase4_result_fixture(
    model = "3pl",
    beta = c(-0.5, 0.5),
    lambda = c(0.9, 1.1),
    guessing = c(0.1, 0.2)
  )

  bad_beta <- result
  bad_beta$items_calib$data$beta[1] <- bad_beta$items_calib$data$beta[1] + 0.1
  expect_error(
    simulate_response_data(bad_beta, 10),
    "beta|items_calib"
  )

  bad_guessing <- result
  bad_guessing$items_calib$data$guessing[1] <- 0.3
  expect_error(
    simulate_response_data(bad_guessing, 10),
    "guessing|items_calib"
  )

  bad_length <- result
  bad_length$guessing_vec <- 0.1
  expect_error(
    simulate_response_data(bad_length, 10),
    "guessing|length|n_items"
  )

  missing_guessing <- result
  missing_guessing$guessing_vec <- NULL
  missing_guessing$items_base$data$guessing <- NULL
  missing_guessing$items_calib$data$guessing <- NULL
  expect_error(
    simulate_response_data(missing_guessing, 10),
    "3PL|3pl|guessing"
  )
})


test_that("literal v0.1 SPC design loss gives an actionable immutable error", {
  literal_spc <- .phase4_literal_v01_spc()
  before <- serialize(literal_spc, NULL)

  expect_error(
    simulate_response_data(literal_spc, n_persons = 10, seed = 5001),
    "[Rr]e-?run.*calibration|[Rr]erun.*calibration"
  )
  expect_identical(serialize(literal_spc, NULL), before)
  expect_false("beta_vec" %in% names(literal_spc))
  expect_false("lambda_scaled" %in% names(literal_spc))
})
