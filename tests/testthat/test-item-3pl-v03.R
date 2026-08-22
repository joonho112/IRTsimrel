# Version 0.3 item-generation tests: 3PL guessing, source contracts, and S3 schema

test_that("sim_item_params appends guessing_params without moving legacy formals", {
  legacy_formals <- c(
    "n_items", "model", "source", "method", "n_forms",
    "difficulty_params", "discrimination_params", "hierarchical_params",
    "custom_params", "scale", "center_difficulties", "seed"
  )
  expect_identical(names(formals(sim_item_params))[seq_along(legacy_formals)],
                   legacy_formals)
  expect_identical(tail(names(formals(sim_item_params)), 1L), "guessing_params")
})

test_that("3PL fixed guessing defaults to .20 and accepts scalar or vector", {
  default_items <- sim_item_params(7L, model = "3pl", seed = 1)
  expect_identical(
    names(default_items$data),
    c("form_id", "item_id", "beta", "lambda", "lambda_unscaled", "guessing")
  )
  expect_identical(default_items$data$guessing, rep(0.2, 7L))
  expect_identical(default_items$params$guessing$distribution, "fixed")

  guessing <- seq(0, 0.3, length.out = 7L)
  vector_items <- sim_item_params(
    7L,
    model = "3pl",
    guessing_params = list(distribution = "fixed", value = guessing),
    seed = 1
  )
  expect_identical(vector_items$data$guessing, guessing)
})

test_that("beta and uniform guessing generators use documented defaults and bounds", {
  beta_items <- sim_item_params(
    200L,
    model = "3pl",
    guessing_params = list(distribution = "beta"),
    seed = 23
  )
  expect_identical(
    beta_items$params$guessing,
    list(distribution = "beta", shape1 = 5, shape2 = 17)
  )
  expect_true(all(beta_items$data$guessing >= 0 & beta_items$data$guessing < 1))

  uniform_items <- sim_item_params(
    100L,
    model = "3pl",
    guessing_params = list(distribution = "uniform"),
    seed = 23
  )
  expect_identical(
    uniform_items$params$guessing,
    list(distribution = "uniform", min = 0.1, max = 0.3)
  )
  expect_true(all(uniform_items$data$guessing >= 0.1))
  expect_true(all(uniform_items$data$guessing <= 0.3))
})

test_that("guessing validation rejects malformed values without clipping", {
  expect_error(
    sim_item_params(4L, model = "3pl", guessing_params = list(0.2)),
    "named"
  )
  expect_error(
    sim_item_params(
      4L, model = "3pl",
      guessing_params = list(distribution = "fixed", value = c(0.1, 0.2))
    ),
    "length 1 or n_items"
  )
  expect_error(
    sim_item_params(
      4L, model = "3pl",
      guessing_params = list(distribution = "fixed", value = 1)
    ),
    "0 <= guessing < 1"
  )
  expect_error(
    sim_item_params(
      4L, model = "3pl",
      guessing_params = list(distribution = "beta", shape1 = 0)
    ),
    "shape1"
  )
  expect_error(
    sim_item_params(
      4L, model = "3pl",
      guessing_params = list(distribution = "uniform", min = 0.4, max = 0.2)
    ),
    "0 <= min < max < 1"
  )
  expect_error(
    sim_item_params(
      4L, model = "3pl",
      guessing_params = list(distribution = "fixed", typo = 0.2)
    ),
    "Unknown or inapplicable"
  )
  expect_warning(
    high <- sim_item_params(
      4L, model = "3pl",
      guessing_params = list(distribution = "fixed", value = 0.6),
      seed = 1
    ),
    "at least 0.5"
  )
  expect_identical(high$data$guessing, rep(0.6, 4L))
})

test_that("Rasch and 2PL reject nonempty guessing inputs", {
  for (model in c("rasch", "2pl")) {
    expect_error(
      sim_item_params(
        5L, model = model,
        guessing_params = list(distribution = "fixed", value = 0.2)
      ),
      "only supported.*3pl"
    )
    expect_error(
      sim_item_params(
        5L, model = model, source = "custom",
        custom_params = list(
          beta = rep(0, 5L),
          lambda = rep(1, 5L),
          guessing = 0.2
        )
      ),
      "only supported.*3pl"
    )
  }
})

test_that("custom 3PL requires and retains beta, lambda, and guessing", {
  beta <- seq(-1, 1, length.out = 5L)
  lambda <- seq(0.8, 1.2, length.out = 5L)
  custom <- list(beta = beta, lambda = lambda, guessing = 0.15)
  items <- sim_item_params(
    5L,
    model = "3pl",
    source = "custom",
    custom_params = custom,
    center_difficulties = FALSE,
    scale = 2
  )

  expect_identical(items$data$beta, beta)
  expect_identical(items$data$lambda_unscaled, lambda)
  expect_equal(items$data$lambda, 2 * lambda)
  expect_identical(items$data$guessing, rep(0.15, 5L))
  expect_identical(items$params$custom, custom)
  expect_identical(items$params$guessing, list(distribution = "custom"))
  expect_true(is.na(items$method))

  function_items <- sim_item_params(
    5L,
    model = "3pl",
    source = "custom",
    custom_params = list(
      beta = function(n) seq_len(n),
      lambda = function(n) rep(1.1, n),
      guessing = function(n) seq(0.1, 0.2, length.out = n)
    ),
    center_difficulties = FALSE
  )
  expect_equal(function_items$data$guessing, seq(0.1, 0.2, length.out = 5L))

  expect_error(
    sim_item_params(
      5L, model = "3pl", source = "custom",
      custom_params = list(beta = beta, lambda = lambda)
    ),
    "custom_params\\$guessing.*required"
  )
  expect_error(
    sim_item_params(
      5L, model = "3pl", source = "custom",
      custom_params = custom,
      guessing_params = list(distribution = "fixed", value = 0.2)
    ),
    "ambiguous"
  )
})

test_that("3PL and 2PL share beta and discrimination generation", {
  twopl <- sim_item_params(30L, model = "2pl", source = "parametric", seed = 92)
  threepl <- sim_item_params(30L, model = "3pl", source = "parametric", seed = 92)
  expect_identical(threepl$data$beta, twopl$data$beta)
  expect_identical(threepl$data$lambda, twopl$data$lambda)

  twopl_h <- sim_item_params(30L, model = "2pl", source = "hierarchical", seed = 92)
  threepl_h <- sim_item_params(30L, model = "3pl", source = "hierarchical", seed = 92)
  expect_identical(threepl_h$data$beta, twopl_h$data$beta)
  expect_identical(threepl_h$data$lambda, twopl_h$data$lambda)
})

test_that("hierarchical source handles one item and Rasch keeps unit base lambda", {
  one <- sim_item_params(
    1L,
    model = "3pl",
    source = "hierarchical",
    guessing_params = list(distribution = "fixed", value = 0.25),
    seed = 9
  )
  expect_equal(nrow(one$data), 1L)
  expect_true(all(is.finite(one$data$beta)))
  expect_true(all(one$data$lambda > 0))

  rasch <- sim_item_params(15L, model = "rasch", source = "hierarchical", seed = 9)
  expect_identical(rasch$data$lambda_unscaled, rep(1, 15L))
  expect_identical(rasch$data$lambda, rep(1, 15L))

  hierarchical <- sim_item_params(
    15L,
    model = "2pl",
    source = "hierarchical",
    hierarchical_params = list(rho = -0.72),
    seed = 9
  )
  output <- capture.output(print(hierarchical))
  expect_true(any(grepl("Target \\(rho\\): -0.7200", output)))
  expect_true(is.na(hierarchical$method))
})

test_that("scale changes only lambda for generated 3PL items", {
  base <- sim_item_params(
    25L,
    model = "3pl",
    guessing_params = list(distribution = "beta"),
    scale = 1,
    seed = 892
  )
  scaled <- sim_item_params(
    25L,
    model = "3pl",
    guessing_params = list(distribution = "beta"),
    scale = 2.4,
    seed = 892
  )
  expect_identical(scaled$data$beta, base$data$beta)
  expect_identical(scaled$data$lambda_unscaled, base$data$lambda_unscaled)
  expect_equal(scaled$data$lambda, 2.4 * base$data$lambda_unscaled)
  expect_identical(scaled$data$guessing, base$data$guessing)

  adapted <- IRTsimrel:::.irtsimrel_apply_item_scale(
    base,
    lambda_base = base$data$lambda_unscaled,
    scale = 3
  )
  expect_identical(adapted$data$guessing, base$data$guessing)
})

test_that("seeded guessing generation restores caller RNG state", {
  set.seed(191)
  state_before <- .Random.seed
  invisible(sim_item_params(
    20L,
    model = "3pl",
    guessing_params = list(distribution = "uniform"),
    seed = 300
  ))
  expect_identical(.Random.seed, state_before)

  first <- sim_item_params(
    20L,
    model = "3pl",
    guessing_params = list(distribution = "beta"),
    seed = 300
  )
  second <- sim_item_params(
    20L,
    model = "3pl",
    guessing_params = list(distribution = "beta"),
    seed = 300
  )
  expect_identical(first$data, second$data)
})

test_that("legacy item schemas and S3 output remain model-specific", {
  expected_legacy <- c("form_id", "item_id", "beta", "lambda", "lambda_unscaled")
  rasch <- sim_item_params(8L, model = "rasch", seed = 6)
  twopl <- sim_item_params(8L, model = "2pl", seed = 6)
  threepl <- sim_item_params(8L, model = "3pl", seed = 6)

  expect_identical(names(rasch$data), expected_legacy)
  expect_identical(names(twopl$data), expected_legacy)
  expect_identical(names(threepl$data), c(expected_legacy, "guessing"))
  expect_false("guessing_summary" %in% names(summary(rasch)))
  expect_false("guessing_summary" %in% names(summary(twopl)))
  expect_true("guessing_summary" %in% names(summary(threepl)))
  expect_false(any(grepl("Guessing", capture.output(print(twopl)))))
  expect_true(any(grepl("Guessing", capture.output(print(threepl)))))
  expect_identical(as.data.frame(threepl), threepl$data)

  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(threepl, type = "scatter"), "gg")
})
