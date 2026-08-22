# Phase 6 SAC S3 contracts ---------------------------------------------------

.sac_s3_beta <- c(-1, -0.3, 0.4, 1)
.sac_s3_lambda <- c(0.8, 1, 1.2, 1.4)
.sac_s3_guessing <- c(0.1, 0.15, 0.2, 0.25)

.sac_s3_params <- function(model = c("2pl", "3pl")) {
  model <- match.arg(model)
  custom <- list(beta = .sac_s3_beta, lambda = .sac_s3_lambda)
  if (identical(model, "3pl")) custom$guessing <- .sac_s3_guessing
  list(custom_params = custom, center_difficulties = FALSE)
}

.sac_s3_fit <- function(model = c("2pl", "3pl"),
                        resample_items = FALSE,
                        seed = 9901L) {
  model <- match.arg(model)
  suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.3,
    n_items = length(.sac_s3_beta),
    model = model,
    item_source = "custom",
    item_params = .sac_s3_params(model),
    reliability_metric = "info",
    c_init = 0.8,
    M_per_iter = 24L,
    M_pre = 48L,
    n_iter = 4L,
    burn_in = 2L,
    c_bounds = c(0.1, 4),
    resample_items = resample_items,
    seed = seed,
    preflight_controls = list(
      n_forms = if (resample_items) 4L else 1L,
      M = 48L,
      split_check = resample_items,
      split_log_tolerance = 1,
      root_controls = list(min_grid = 33L, max_evals = 129L)
    ),
    evaluation_controls = list(n_forms = 3L, M = 48L)
  )))
}


test_that("SAC method formals retain their public contracts", {
  expect_identical(names(formals(print.sac_result)), c("x", "digits", "..."))
  expect_identical(names(formals(summary.sac_result)), c("object", "..."))
  expect_identical(
    names(formals(print.summary.sac_result)),
    c("x", "digits", "...")
  )
  expect_identical(names(formals(coef.sac_result)), c("object", "..."))
  expect_identical(
    names(formals(predict.sac_result)),
    c("object", "newdata", "theta_vec", "...")
  )
})


test_that("fixed-form 3PL coef and predict preserve and use guessing", {
  fit <- .sac_s3_fit("3pl", resample_items = FALSE, seed = 9902L)
  coefficients <- coef(fit)

  expect_identical(
    names(coefficients),
    c(
      "item_id", "beta", "lambda_base", "lambda_scaled", "c_star",
      "guessing"
    )
  )
  expect_equal(coefficients$guessing, .sac_s3_guessing, tolerance = 0)

  theta <- seq(-2.5, 2.5, length.out = 101L)
  scales <- c(0.55, 0.9, 1.25)
  observed <- predict(fit, newdata = scales, theta_vec = theta)
  expected <- vapply(
    scales,
    function(scale) {
      .compute_rho_generic(
        c = scale,
        theta_vec = theta,
        beta_vec = fit$beta_vec,
        lambda_base = fit$lambda_base,
        theta_var = stats::var(theta),
        metric_internal = fit$metric,
        guessing = fit$guessing_vec
      )
    },
    numeric(1)
  )
  expect_equal(as.numeric(observed), expected, tolerance = 1e-13)
  expect_equal(
    unname(predict(fit, newdata = fit$c_star)),
    fit$achieved_rho,
    tolerance = 1e-13
  )
  expect_identical(predict(fit), fit$achieved_rho)
})


test_that("superpopulation prediction averages every stored holdout form", {
  fit <- .sac_s3_fit("3pl", resample_items = TRUE, seed = 9903L)
  scales <- c(0.6, fit$c_star, 1.1)
  observed <- predict(fit, newdata = scales)

  banks <- fit$evaluation_design$item_banks
  theta_blocks <- fit$evaluation_design$theta_blocks
  expected <- vapply(
    scales,
    function(scale) {
      mean(vapply(
        seq_along(banks),
        function(index) {
          bank <- banks[[index]]
          .compute_rho_generic(
            c = scale,
            theta_vec = theta_blocks[[index]],
            beta_vec = bank$beta,
            lambda_base = bank$lambda_base,
            theta_var = fit$theta_var,
            metric_internal = fit$metric,
            guessing = bank$guessing
          )
        },
        numeric(1)
      ))
    },
    numeric(1)
  )

  expect_equal(as.numeric(observed), expected, tolerance = 1e-13)
  expect_identical(attr(observed, "prediction_scope"), "item_superpopulation")
  expect_equal(unname(observed[[2L]]), fit$achieved_rho, tolerance = 1e-13)

  common_theta <- seq(-2, 2, length.out = 81L)
  common <- predict(fit, newdata = scales[[1L]], theta_vec = common_theta)
  common_expected <- mean(vapply(
    banks,
    function(bank) {
      .compute_rho_generic(
        c = scales[[1L]],
        theta_vec = common_theta,
        beta_vec = bank$beta,
        lambda_base = bank$lambda_base,
        theta_var = stats::var(common_theta),
        metric_internal = fit$metric,
        guessing = bank$guessing
      )
    },
    numeric(1)
  ))
  expect_equal(as.numeric(common), common_expected, tolerance = 1e-13)
})


test_that("legacy 2PL methods use lazy zero guessing without mutation", {
  legacy <- .sac_s3_fit("2pl", resample_items = FALSE, seed = 9904L)
  legacy$guessing_vec <- NULL
  legacy$schema_version <- NULL
  legacy$item_scope <- NULL
  legacy$estimand_signature <- NULL
  legacy$design_signature <- NULL
  legacy$items_base$data$guessing <- NULL
  legacy$items_calib$data$guessing <- NULL
  before <- serialize(legacy, NULL)

  coefficients <- coef(legacy)
  expect_identical(
    names(coefficients),
    c("item_id", "beta", "lambda_base", "lambda_scaled", "c_star")
  )

  theta <- seq(-2, 2, length.out = 81L)
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
  expect_false("guessing_summary" %in% names(summary(legacy)))
  expect_identical(serialize(legacy, NULL), before)
})


test_that("SAC print and summary expose additive v0.3 diagnostics", {
  fit <- .sac_s3_fit("3pl", resample_items = FALSE, seed = 9905L)
  summary_object <- summary(fit)
  legacy_prefix <- c(
    "c_star", "target_rho", "achieved_rho", "metric", "model",
    "n_items", "n_iter", "M_per_iter", "M_pre", "theta_var",
    "init_method", "convergence", "burn_in"
  )

  expect_identical(
    names(summary_object)[seq_along(legacy_prefix)],
    legacy_prefix
  )
  expect_identical(summary_object$item_scope, "fixed_form")
  expect_equal(summary_object$achieved_se, fit$achieved_se, tolerance = 0)
  expect_identical(
    summary_object$achieved_distribution_n,
    fit$achieved_distribution$n_forms
  )
  expect_equal(summary_object$selected_root, fit$branch$selected_root)
  expect_identical(summary_object$selected_branch, fit$branch$selected_branch)
  expect_identical(summary_object$calibration_status, fit$calibration_status)
  expect_identical(
    summary_object$guessing_summary,
    c(
      min = min(.sac_s3_guessing),
      mean = mean(.sac_s3_guessing),
      max = max(.sac_s3_guessing)
    )
  )

  output <- capture.output(returned <- print(fit))
  expect_identical(returned, fit)
  expect_true(any(grepl("Item scope", output, fixed = TRUE)))
  expect_true(any(grepl("Achieved reliability SE", output, fixed = TRUE)))
  expect_true(any(grepl("Selected preflight root", output, fixed = TRUE)))
  expect_true(any(grepl("Guessing", output, fixed = TRUE)))

  summary_output <- capture.output(returned_summary <- print(summary_object))
  expect_identical(returned_summary, summary_object)
  expect_true(any(grepl("Evaluation forms", summary_output, fixed = TRUE)))
  expect_true(any(grepl("Status", summary_output, fixed = TRUE)))
})
