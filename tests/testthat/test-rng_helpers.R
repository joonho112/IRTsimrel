# =============================================================================
# test-rng_helpers.R
# =============================================================================
# Internal RNG helpers and public seed side effects
# =============================================================================

capture_rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    list(exists = TRUE, value = get(".Random.seed", envir = globalenv()))
  } else {
    list(exists = FALSE, value = NULL)
  }
}

restore_rng_state <- function(state) {
  if (state$exists) {
    assign(".Random.seed", state$value, envir = globalenv())
  } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
}

expect_preserves_rng <- function(expr) {
  set.seed(24601)
  state_before <- .Random.seed
  force(expr)
  expect_identical(.Random.seed, state_before)
}

expect_preserves_rng_on_error <- function(expr, regexp) {
  set.seed(13579)
  state_before <- .Random.seed
  expect_error(force(expr), regexp)
  expect_identical(.Random.seed, state_before)
}

test_that("internal seed helper restores an existing RNG state", {
  original_state <- capture_rng_state()
  on.exit(restore_rng_state(original_state), add = TRUE)

  set.seed(1234)
  state_before <- .Random.seed
  restore_seed <- IRTsimrel:::.irtsimrel_set_seed(99)

  expect_true(is.function(restore_seed))
  invisible(runif(5))
  restore_seed()

  expect_identical(.Random.seed, state_before)
})

test_that("internal seed helper restores missing RNG state", {
  original_state <- capture_rng_state()
  on.exit(restore_rng_state(original_state), add = TRUE)

  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))

  restore_seed <- IRTsimrel:::.irtsimrel_set_seed(99)
  expect_true(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  invisible(runif(5))
  restore_seed()

  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("internal seed helper accepts signed 32-bit boundary seeds", {
  original_state <- capture_rng_state()
  on.exit(restore_rng_state(original_state), add = TRUE)

  set.seed(4321)
  state_before <- .Random.seed

  for (seed in c(-.Machine$integer.max, 0, .Machine$integer.max)) {
    restore_seed <- IRTsimrel:::.irtsimrel_set_seed(seed)
    expect_true(is.function(restore_seed))
    invisible(runif(1))
    restore_seed()
    expect_identical(.Random.seed, state_before)
  }
})

test_that("internal seed helper validates seed before changing RNG state", {
  original_state <- capture_rng_state()
  on.exit(restore_rng_state(original_state), add = TRUE)

  set.seed(4321)
  state_before <- .Random.seed

  expect_null(IRTsimrel:::.irtsimrel_set_seed(NULL))
  expect_identical(.Random.seed, state_before)

  expect_error(IRTsimrel:::.irtsimrel_set_seed("99"), "`seed` must be a single integer value")
  expect_error(IRTsimrel:::.irtsimrel_set_seed(c(1, 2)), "`seed` must be a single integer value")
  expect_error(IRTsimrel:::.irtsimrel_set_seed(NA_real_), "`seed` must be a single integer value")
  expect_error(IRTsimrel:::.irtsimrel_set_seed(1.5), "`seed` must be a single integer value")
  expect_error(IRTsimrel:::.irtsimrel_set_seed(.Machine$integer.max + 1),
               "signed 32-bit integer range")
  expect_error(IRTsimrel:::.irtsimrel_set_seed(-.Machine$integer.max - 1),
               "signed 32-bit integer range")
  expect_identical(.Random.seed, state_before)
})

test_that("seeded public generators preserve existing RNG state", {
  expect_preserves_rng(
    sim_latentG(n = 30L, shape = "normal", seed = 1)
  )

  expect_preserves_rng(
    sim_item_params(n_items = 8L, model = "rasch", source = "parametric", seed = 1)
  )

  result <- suppressWarnings(suppressMessages(eqc_calibrate(
    target_rho = 0.70,
    n_items = 6L,
    item_source = "parametric",
    M = 100L,
    seed = 1
  )))

  expect_preserves_rng(
    simulate_response_data(result, n_persons = 20L, seed = 1)
  )
})

test_that("seeded public generators restore missing RNG state", {
  original_state <- capture_rng_state()
  on.exit(restore_rng_state(original_state), add = TRUE)

  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv(), inherits = FALSE)
  }

  sim_item_params(n_items = 8L, model = "rasch", source = "parametric", seed = 1)

  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("seeded public calls restore RNG state after downstream errors", {
  expect_preserves_rng_on_error(
    suppressMessages(eqc_calibrate(
      target_rho = 0.70,
      n_items = 6L,
      item_source = "parametric",
      latent_shape = "custom",
      M = 100L,
      seed = 1,
      verbose = FALSE
    )),
    "mixture_spec"
  )
})

test_that("unseeded public generators consume caller RNG state", {
  set.seed(777)
  state_before <- .Random.seed
  sim_latentG(n = 30L, shape = "normal")

  expect_false(identical(.Random.seed, state_before))
})

test_that("seeded calibration utilities preserve existing RNG state", {
  expect_preserves_rng(
    suppressMessages(eqc_calibrate(
      target_rho = 0.70,
      n_items = 6L,
      item_source = "parametric",
      M = 500L,
      seed = 1,
      verbose = FALSE
    ))
  )

  expect_preserves_rng(
    suppressWarnings(suppressMessages(sac_calibrate(
      target_rho = 0.70,
      n_items = 6L,
      item_source = "parametric",
      M_per_iter = 60L,
      M_pre = 120L,
      n_iter = 8L,
      burn_in = 4L,
      seed = 1,
      verbose = FALSE
    )))
  )

  expect_preserves_rng(
    suppressMessages(check_feasibility(
      n_items = 6L,
      M = 300L,
      seed = 1,
      verbose = FALSE
    ))
  )

  expect_preserves_rng(
    rho_curve(
      c_values = c(0.5, 1, 1.5),
      n_items = 6L,
      M = 300L,
      seed = 1,
      plot = FALSE
    )
  )
})

test_that("SAC with same seed produces identical trajectories", {
  sac1 <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.70,
    n_items = 8L,
    item_source = "parametric",
    M_per_iter = 80L,
    M_pre = 200L,
    n_iter = 12L,
    burn_in = 6L,
    seed = 42,
    verbose = FALSE
  )))

  sac2 <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.70,
    n_items = 8L,
    item_source = "parametric",
    M_per_iter = 80L,
    M_pre = 200L,
    n_iter = 12L,
    burn_in = 6L,
    seed = 42,
    verbose = FALSE
  )))

  expect_identical(sac1$trajectory, sac2$trajectory)
  expect_identical(sac1$rho_trajectory, sac2$rho_trajectory)
  expect_equal(sac1$c_star, sac2$c_star, tolerance = 1e-12)
  expect_equal(sac1$achieved_rho, sac2$achieved_rho, tolerance = 1e-12)
})
