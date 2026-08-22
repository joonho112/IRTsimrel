canonical_v03_item_params <- function() {
  list(custom_params = list(
    beta = c(-1, -0.3, 0.3, 1),
    lambda = c(0.8, 1, 1.2, 1.4)
  ))
}


test_that("latent-shape partial matches are canonicalized once", {
  expect_identical(.irtsimrel_match_latent_shape("norm"), "normal")
  expect_identical(.irtsimrel_match_latent_shape("heavy"), "heavy_tail")
  expect_identical(.irtsimrel_match_latent_shape("light"), "light_tail")
  expect_error(.irtsimrel_match_latent_shape("skew"), "arg")
  expect_error(.irtsimrel_match_latent_shape(NA_character_), "non-empty")
})


test_that("partial heavy-tail names cannot bypass the MSEM guard", {
  condition <- tryCatch(
    sac_calibrate(
      target_rho = 0.3,
      n_items = 4L,
      model = "2pl",
      latent_shape = "heavy",
      item_source = "custom",
      item_params = canonical_v03_item_params(),
      reliability_metric = "msem",
      M_per_iter = 16L,
      M_pre = 32L,
      n_iter = 2L,
      burn_in = 1L,
      resample_items = FALSE,
      seed = 9903L,
      preflight_controls = list(M = 32L, split_check = FALSE),
      evaluation_controls = list(n_forms = 2L, M = 32L)
    ),
    error = identity
  )

  expect_s3_class(condition, "irtsimrel_sac_nonintegrable")
  expect_match(conditionMessage(condition), "heavy-tail|not finite")
})


test_that("public latent-shape orchestrators store canonical contracts", {
  item_params <- canonical_v03_item_params()

  eqc <- suppressMessages(eqc_calibrate(
    target_rho = 0.4,
    n_items = 4L,
    model = "2pl",
    latent_shape = "norm",
    item_source = "custom",
    item_params = item_params,
    M = 128L,
    c_bounds = c(0.1, 3),
    seed = 9901L,
    root_controls = list(min_grid = 33L, max_evals = 129L)
  ))
  expect_identical(eqc$design_signature$latent_shape, "normal")

  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.4,
    n_items = 4L,
    model = "2pl",
    latent_shape = "norm",
    reliability_metric = "info",
    c_init = eqc,
    M_per_iter = 32L,
    M_pre = 64L,
    n_iter = 4L,
    burn_in = 2L,
    c_bounds = c(0.1, 3),
    resample_items = FALSE,
    seed = 9902L,
    preflight_controls = list(
      M = 48L,
      split_check = FALSE,
      root_controls = list(min_grid = 33L, max_evals = 129L)
    ),
    evaluation_controls = list(n_forms = 2L, M = 48L)
  )))
  expect_identical(sac$design_signature$latent_shape, "normal")
  expect_identical(sac$calibration_design$latent_shape, "normal")
  expect_true(sac$warm_start$comparable)

  feasibility <- check_feasibility(
    n_items = 4L,
    model = "2pl",
    latent_shape = "norm",
    item_source = "custom",
    item_params = item_params,
    M = 64L,
    c_bounds = c(0.2, 2),
    seed = 9904L,
    verbose = FALSE
  )
  expect_identical(feasibility$latent_shape, "normal")

  curve <- rho_curve(
    c_values = c(0.5, 1, 1.5),
    n_items = 4L,
    model = "2pl",
    latent_shape = "norm",
    item_source = "custom",
    item_params = item_params,
    M = 64L,
    seed = 9905L,
    plot = FALSE
  )
  expect_identical(attr(curve, "latent_shape"), "normal")

  response <- simulate_response_data(
    eqc,
    n_persons = 10L,
    latent_shape = "norm",
    seed = 9906L
  )
  expect_identical(response$provenance$latent_shape, "normal")
})
