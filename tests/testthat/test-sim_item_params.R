# =============================================================================
# test-sim_item_params.R
# =============================================================================
# Comprehensive tests for sim_item_params() and related functions
# =============================================================================

# ---------------------------------------------------------------------------
# Rasch parametric: all discriminations equal 1
# ---------------------------------------------------------------------------

test_that("Rasch parametric: all lambda_unscaled = 1 and lambda = 1", {
  items <- sim_item_params(
    n_items = 25, model = "rasch", source = "parametric", seed = 42
  )
  expect_true(all(items$data$lambda_unscaled == 1))
  expect_true(all(items$data$lambda == 1))
})


# ---------------------------------------------------------------------------
# 2PL parametric with different methods
# ---------------------------------------------------------------------------

describe("2PL parametric methods for discrimination generation", {

  test_that("copula method: negative rho produces negative Spearman correlation", {
    items <- sim_item_params(
      n_items = 100, model = "2pl", source = "parametric",
      method = "copula",
      discrimination_params = list(rho = -0.3),
      seed = 42
    )
    # Spearman correlation between beta and log(lambda) should be negative
    spearman_r <- cor(items$data$beta, log(items$data$lambda_unscaled),
                      method = "spearman")
    expect_true(spearman_r < 0,
                info = paste("Spearman r =", round(spearman_r, 3)))
  })

  test_that("copula method: positive rho produces positive Spearman correlation", {
    items <- sim_item_params(
      n_items = 100, model = "2pl", source = "parametric",
      method = "copula",
      discrimination_params = list(rho = 0.5),
      seed = 42
    )
    spearman_r <- cor(items$data$beta, log(items$data$lambda_unscaled),
                      method = "spearman")
    expect_true(spearman_r > 0,
                info = paste("Spearman r =", round(spearman_r, 3)))
  })

  test_that("conditional method: negative rho produces negative correlation", {
    items <- sim_item_params(
      n_items = 100, model = "2pl", source = "parametric",
      method = "conditional",
      discrimination_params = list(rho = -0.5),
      seed = 42
    )
    obs_cor <- cor(items$data$beta, log(items$data$lambda_unscaled))
    expect_true(obs_cor < 0,
                info = paste("Pearson r =", round(obs_cor, 3)))
  })

  test_that("independent method: no systematic correlation", {
    items <- sim_item_params(
      n_items = 200, model = "2pl", source = "parametric",
      method = "independent",
      seed = 42
    )
    obs_cor <- cor(items$data$beta, log(items$data$lambda_unscaled))
    # With 200 independent items, correlation should be near zero
    expect_true(abs(obs_cor) < 0.25,
                info = paste("Pearson r =", round(obs_cor, 3)))
  })
})


# ---------------------------------------------------------------------------
# Hierarchical 2PL source
# ---------------------------------------------------------------------------

test_that("hierarchical 2PL: all lambda > 0 and reasonable values", {
  items <- sim_item_params(
    n_items = 30, model = "2pl", source = "hierarchical",
    hierarchical_params = list(mu = c(0, 0), tau = c(0.25, 1), rho = -0.3),
    seed = 42
  )
  expect_true(all(items$data$lambda > 0))
  expect_true(all(is.finite(items$data$beta)))
  expect_equal(nrow(items$data), 30)
})


# ---------------------------------------------------------------------------
# Custom source: vectors and functions
# ---------------------------------------------------------------------------

describe("Custom source for item parameters", {

  test_that("custom with exact beta and lambda vectors", {
    beta_custom <- seq(-2, 2, length.out = 10)
    lambda_custom <- rep(1.5, 10)
    items <- sim_item_params(
      n_items = 10, model = "2pl", source = "custom",
      custom_params = list(beta = beta_custom, lambda = lambda_custom),
      center_difficulties = FALSE
    )
    expect_equal(items$data$lambda_unscaled, lambda_custom)
    expect_equal(items$data$beta, beta_custom)
  })

  test_that("custom with beta as a function returns correct n_items output", {
    beta_fn <- function(n) seq(-1, 1, length.out = n)
    items <- sim_item_params(
      n_items = 15, model = "rasch", source = "custom",
      custom_params = list(beta = beta_fn),
      center_difficulties = FALSE
    )
    expected_beta <- seq(-1, 1, length.out = 15)
    expect_equal(items$data$beta, expected_beta, tolerance = 1e-10)
    expect_equal(nrow(items$data), 15)
  })

  test_that("custom with lambda as a function", {
    lambda_fn <- function(n) rep(2.0, n)
    items <- sim_item_params(
      n_items = 10, model = "2pl", source = "custom",
      custom_params = list(
        beta = rep(0, 10),
        lambda = lambda_fn
      ),
      center_difficulties = FALSE
    )
    expect_equal(items$data$lambda_unscaled, rep(2.0, 10))
  })
})


# ---------------------------------------------------------------------------
# center_difficulties
# ---------------------------------------------------------------------------

describe("Centering of difficulties", {

  test_that("center_difficulties=TRUE centers beta to mean ~ 0", {
    items <- sim_item_params(
      n_items = 50, model = "rasch", source = "parametric",
      center_difficulties = TRUE, seed = 42
    )
    expect_true(abs(mean(items$data$beta)) < 1e-10)
  })

  test_that("center_difficulties=FALSE preserves original mean", {
    items <- sim_item_params(
      n_items = 10, model = "rasch", source = "custom",
      custom_params = list(beta = rep(5, 10)),
      center_difficulties = FALSE
    )
    expect_equal(mean(items$data$beta), 5)
  })
})


# ---------------------------------------------------------------------------
# Scale factor
# ---------------------------------------------------------------------------

test_that("scale=2.5 multiplies lambda_unscaled by 2.5", {
  items <- sim_item_params(
    n_items = 20, model = "2pl", source = "parametric",
    scale = 2.5, seed = 42
  )
  ratio <- items$data$lambda / items$data$lambda_unscaled
  expect_true(all(abs(ratio - 2.5) < 1e-10))
})

test_that("Rasch with scale=2: lambda = 2 * lambda_unscaled = 2", {
  items <- sim_item_params(
    n_items = 20, model = "rasch", source = "parametric",
    scale = 2, seed = 42
  )
  expect_true(all(items$data$lambda == 2))
  expect_true(all(items$data$lambda_unscaled == 1))
})


# ---------------------------------------------------------------------------
# Multiple forms
# ---------------------------------------------------------------------------

test_that("n_forms=5 produces 5 * n_items rows with 5 unique form_ids", {
  items <- sim_item_params(
    n_items = 10, model = "rasch", source = "parametric",
    n_forms = 5, seed = 42
  )
  expect_equal(nrow(items$data), 50)  # 5 * 10
  expect_equal(length(unique(items$data$form_id)), 5)
  # Each form should have exactly 10 items
  form_counts <- table(items$data$form_id)
  expect_true(all(form_counts == 10))
})


# ---------------------------------------------------------------------------
# Seed reproducibility
# ---------------------------------------------------------------------------

test_that("same seed produces identical data frames", {
  items1 <- sim_item_params(
    n_items = 20, model = "2pl", source = "parametric",
    method = "copula", seed = 123
  )
  items2 <- sim_item_params(
    n_items = 20, model = "2pl", source = "parametric",
    method = "copula", seed = 123
  )
  expect_identical(items1$data, items2$data)
})


# ---------------------------------------------------------------------------
# S3 methods: print, summary, as.data.frame
# ---------------------------------------------------------------------------

describe("S3 methods for item_params", {

  test_that("print.item_params runs without error and returns invisibly", {
    items <- sim_item_params(
      n_items = 10, model = "2pl", source = "parametric", seed = 1
    )
    out <- expect_invisible(print(items))
    expect_identical(out, items)
  })

  test_that("summary.item_params returns summary.item_params object", {
    items <- sim_item_params(
      n_items = 10, model = "rasch", source = "parametric", seed = 1
    )
    s <- summary(items)
    expect_s3_class(s, "summary.item_params")
    expect_equal(s$model, "rasch")
    expect_equal(s$n_items, 10)
  })

  test_that("as.data.frame.item_params returns data frame with correct columns", {
    items <- sim_item_params(
      n_items = 15, model = "2pl", source = "parametric", seed = 1
    )
    df <- as.data.frame(items)
    expect_s3_class(df, "data.frame")
    expect_true(all(c("form_id", "item_id", "beta", "lambda",
                       "lambda_unscaled") %in% names(df)))
    expect_equal(nrow(df), 15)
  })
})


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

describe("Input validation for sim_item_params", {

  test_that("n_items=0 errors", {
    expect_error(
      sim_item_params(n_items = 0, model = "rasch", source = "parametric"),
      "`n_items` must be a positive integer"
    )
  })

  test_that("n_items=-5 errors", {
    expect_error(
      sim_item_params(n_items = -5, model = "rasch", source = "parametric"),
      "`n_items` must be a positive integer"
    )
  })

  test_that("invalid model errors", {
    expect_error(
      sim_item_params(n_items = 10, model = "3pl", source = "parametric"),
      "'arg' should be one of"
    )
  })

  test_that("scale=0 errors", {
    expect_error(
      sim_item_params(n_items = 10, model = "rasch", source = "parametric", scale = 0),
      "`scale` must be a positive scalar"
    )
  })

  test_that("scale=-1 errors", {
    expect_error(
      sim_item_params(n_items = 10, model = "rasch", source = "parametric", scale = -1),
      "`scale` must be a positive scalar"
    )
  })

  test_that("custom 2pl without lambda errors", {
    expect_error(
      sim_item_params(
        n_items = 5, model = "2pl", source = "custom",
        custom_params = list(beta = rep(0, 5))
      ),
      "lambda.*required"
    )
  })

  test_that("custom with wrong-length beta vector errors", {
    expect_error(
      sim_item_params(
        n_items = 5, model = "rasch", source = "custom",
        custom_params = list(beta = rep(0, 3))
      ),
      "must have length"
    )
  })
})


# ---------------------------------------------------------------------------
# Return structure
# ---------------------------------------------------------------------------

describe("Return structure of sim_item_params", {

  test_that("return object has class 'item_params'", {
    items <- sim_item_params(
      n_items = 10, model = "rasch", source = "parametric", seed = 1
    )
    expect_s3_class(items, "item_params")
  })

  test_that("data has correct columns", {
    items <- sim_item_params(
      n_items = 10, model = "2pl", source = "parametric", seed = 1
    )
    expected_cols <- c("form_id", "item_id", "beta", "lambda", "lambda_unscaled")
    expect_true(all(expected_cols %in% names(items$data)))
  })

  test_that("output contains all expected metadata fields", {
    items <- sim_item_params(
      n_items = 10, model = "rasch", source = "parametric", seed = 1
    )
    for (field in c("model", "source", "n_items", "n_forms", "scale",
                     "centered", "params", "achieved")) {
      expect_true(field %in% names(items),
                  info = paste("Missing field:", field))
    }
  })

  test_that("all lambda values are positive for 2PL", {
    items <- sim_item_params(
      n_items = 50, model = "2pl", source = "parametric",
      method = "copula", seed = 42
    )
    expect_true(all(items$data$lambda > 0))
    expect_true(all(items$data$lambda_unscaled > 0))
  })

  test_that("achieved statistics include correlations for 2PL", {
    items <- sim_item_params(
      n_items = 30, model = "2pl", source = "parametric", seed = 42
    )
    expect_false(is.na(items$achieved$overall$cor_pearson_pooled))
    expect_false(is.na(items$achieved$overall$cor_spearman_pooled))
  })

  test_that("achieved correlations are NA for Rasch model", {
    items <- sim_item_params(
      n_items = 20, model = "rasch", source = "parametric", seed = 42
    )
    expect_true(is.na(items$achieved$overall$cor_pearson_pooled))
  })
})


# ---------------------------------------------------------------------------
# IRW source (skip if irw not installed)
# ---------------------------------------------------------------------------

test_that("IRW source requires irw package", {
  skip_if_not_installed("irw")
  items <- sim_item_params(
    n_items = 10, model = "rasch", source = "irw", seed = 42
  )
  expect_s3_class(items, "item_params")
  expect_equal(nrow(items$data), 10)
})
