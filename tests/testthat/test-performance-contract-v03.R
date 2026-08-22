test_that("reliability reducer derives a bounded automatic row chunk", {
  chunk <- IRTsimrel:::.irtsimrel_reducer_chunk_size

  expect_identical(chunk(200000L, 50L), 5000L)
  expect_identical(chunk(100L, 50L), 100L)
  expect_identical(chunk(10L, 500000L), 1L)
  expect_identical(chunk(1000L, 20L, max_cells = 1000L), 50L)

  expect_error(chunk(0L, 20L), "positive integer scalars")
  expect_error(chunk(20L, 0L), "positive integer scalars")
  expect_error(chunk(20L, 2L, max_cells = 1.5), "positive integer scalars")
})


test_that("automatic multi-chunk reduction is exactly unchunked-equivalent", {
  set.seed(77101)
  n_nodes <- 10000L
  n_items <- 30L
  theta <- stats::rnorm(n_nodes)
  bank <- IRTsimrel:::.irtsimrel_item_bank(
    beta = stats::rnorm(n_items),
    lambda_base = exp(stats::rnorm(n_items, 0, 0.2)),
    guessing = seq(0.1, 0.3, length.out = n_items),
    model = "3pl"
  )

  automatic <- IRTsimrel:::.irtsimrel_log_test_information(
    theta,
    bank,
    scale = 1.4
  )
  unchunked <- IRTsimrel:::.irtsimrel_log_test_information(
    theta,
    bank,
    scale = 1.4,
    chunk_size = n_nodes
  )
  explicit_small <- IRTsimrel:::.irtsimrel_log_test_information(
    theta,
    bank,
    scale = 1.4,
    chunk_size = 317L
  )

  expect_identical(automatic, unchunked)
  expect_identical(explicit_small, unchunked)
  expect_identical(
    IRTsimrel:::.irtsimrel_reducer_chunk_size(n_nodes, n_items),
    8333L
  )
})


test_that("public reliability values retain exact chunked results", {
  set.seed(77102)
  n_nodes <- 9000L
  n_items <- 31L
  theta <- stats::rnorm(n_nodes)
  beta <- stats::rnorm(n_items)
  lambda <- exp(stats::rnorm(n_items, 0, 0.25))
  guessing <- seq(0.08, 0.28, length.out = n_items)

  public <- compute_rho_both(
    c = 1.2,
    theta_vec = theta,
    beta_vec = beta,
    lambda_base = lambda,
    theta_var = 1,
    guessing = guessing,
    return_diagnostics = TRUE
  )
  internal_info <- IRTsimrel:::.compute_rho_generic(
    c = 1.2,
    theta_vec = theta,
    beta_vec = beta,
    lambda_base = lambda,
    theta_var = 1,
    metric_internal = "info",
    guessing = guessing,
    chunk_size = n_nodes
  )
  internal_msem <- IRTsimrel:::.compute_rho_generic(
    c = 1.2,
    theta_vec = theta,
    beta_vec = beta,
    lambda_base = lambda,
    theta_var = 1,
    metric_internal = "msem",
    guessing = guessing,
    chunk_size = n_nodes
  )

  expect_identical(public$rho_tilde, internal_info)
  expect_identical(public$rho_bar, internal_msem)
  expect_gte(public$rho_tilde, public$rho_bar)
})
