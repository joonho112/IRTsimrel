# =============================================================================
# Adaptive feasibility and rho-curve topology (v0.3)
# =============================================================================

.topology_fixture <- function(guessing = NULL) {
  list(
    theta = c(
      -2.8945639829, -2.1638196461, -1.1142421522, -0.4808900280,
      -0.4195555953, 0.0231953129, 1.6704952657
    ),
    beta = c(-3.4556089868, -1.7336287633, 0.8669555536),
    lambda = c(4.8618527891, 1.4673837406, 0.0537872782),
    guessing = guessing
  )
}

.mock_topology_design <- function(fixture, model) {
  local_mocked_bindings(
    sim_latentG = function(...) list(theta = fixture$theta),
    sim_item_params = function(...) {
      item_data <- data.frame(
        form_id = 1L,
        item_id = seq_along(fixture$beta),
        beta = fixture$beta,
        lambda = fixture$lambda,
        lambda_unscaled = fixture$lambda
      )
      if (identical(model, "3pl")) {
        item_data$guessing <- fixture$guessing
      }
      list(data = item_data)
    },
    .package = "IRTsimrel",
    .env = parent.frame()
  )
}

test_that("check_feasibility exposes four-root 2PL topology", {
  fixture <- .topology_fixture()
  .mock_topology_design(fixture, "2pl")

  result <- check_feasibility(
    n_items = 3L,
    model = "2pl",
    c_bounds = c(0.005, 500),
    M = length(fixture$theta),
    target_rho = 0.18,
    verbose = FALSE
  )

  expect_s3_class(result$topology_info, "irtsimrel_root_topology")
  expect_equal(result$root_count_info, 4L)
  expect_equal(result$admissible_root_count_info, 2L)
  expect_equal(result$target_status_info_canonical, "feasible")
  expect_equal(result$target_status_info, "feasible")
  expect_equal(result$rho_bounds_info, result$topology_info$rho_bounds)
  expect_equal(result$rho_range_info, result$topology_info$rho_range)
  expect_gt(result$rho_range_info[["upper"]], max(result$rho_bounds_info))
  expect_true(all(c(
    "rho_range_info", "rho_range_msem", "rho_bounds_info",
    "rho_bounds_msem", "rho_info_max_c", "rho_msem_max_c", "target_rho",
    "target_status_info", "target_status_msem", "n_items", "model",
    "latent_shape", "c_bounds", "M", "theta_var"
  ) %in% names(result)))
  expect_named(
    result$best_achievable_info,
    c("c", "rho", "residual", "absolute_error", "location")
  )
})

test_that("check_feasibility passes 3PL guessing into both metric scans", {
  fixture <- .topology_fixture(c(0.20, 0.25, 0.10))
  .mock_topology_design(fixture, "3pl")

  result <- check_feasibility(
    n_items = 3L,
    model = "3pl",
    c_bounds = c(0.005, 500),
    M = length(fixture$theta),
    target_rho = 0.09,
    verbose = FALSE
  )

  expect_equal(result$model, "3pl")
  expect_equal(result$root_count_info, 4L)
  expect_equal(result$target_status_info_canonical, "feasible")
  expect_true(all(result$topology_info$roots$direction %in%
                    c("increasing", "decreasing")))
  expect_equal(
    result$rho_bounds_info,
    c(
      lower = compute_rho_tilde(
        0.005, fixture$theta, fixture$beta, fixture$lambda,
        guessing = fixture$guessing
      ),
      upper = compute_rho_tilde(
        500, fixture$theta, fixture$beta, fixture$lambda,
        guessing = fixture$guessing
      )
    ),
    tolerance = 1e-12
  )
})

test_that("canonical infeasibility statuses coexist with legacy aliases", {
  fixture <- .topology_fixture()
  .mock_topology_design(fixture, "2pl")
  high <- check_feasibility(
    n_items = 3L, model = "2pl", c_bounds = c(0.005, 500),
    M = length(fixture$theta), target_rho = 0.99, verbose = FALSE
  )

  expect_equal(high$target_status_info_canonical, "infeasible_above_range")
  expect_equal(high$target_status_info, "above_upper")
  expect_equal(high$root_count_info, 0L)
  expect_true(is.finite(high$best_achievable_info$c))
  expect_true(is.finite(high$best_achievable_info$rho))
  expect_gt(high$best_achievable_info$absolute_error, 0)
})

test_that("rho_curve g=0 is the exact 2PL special case", {
  beta <- c(-1, -0.2, 0.7)
  lambda <- c(0.7, 1.1, 1.5)
  common <- list(
    c_values = c(0.3, 0.8, 1.4, 2.2),
    n_items = 3L,
    latent_shape = "normal",
    item_source = "custom",
    metric = "both",
    M = 400L,
    seed = 73L,
    plot = FALSE
  )
  curve_2pl <- do.call(
    rho_curve,
    c(common, list(
      model = "2pl",
      item_params = list(
        custom_params = list(beta = beta, lambda = lambda),
        center_difficulties = FALSE
      )
    ))
  )
  curve_3pl <- do.call(
    rho_curve,
    c(common, list(
      model = "3pl",
      item_params = list(
        custom_params = list(
          beta = beta, lambda = lambda, guessing = rep(0, 3L)
        ),
        center_difficulties = FALSE
      )
    ))
  )

  expect_named(curve_2pl, c("c", "rho_tilde", "rho_bar"))
  expect_named(curve_3pl, c("c", "rho_tilde", "rho_bar"))
  expect_equal(curve_3pl, curve_2pl, tolerance = 0, ignore_attr = TRUE)
})

test_that("rho_curve topology is adaptive rather than inferred from display points", {
  fixture <- .topology_fixture()
  .mock_topology_design(fixture, "2pl")

  curve <- rho_curve(
    c_values = c(0.005, 500),
    n_items = 3L,
    model = "2pl",
    metric = "info",
    M = length(fixture$theta),
    plot = FALSE
  )
  topology <- attr(curve, "topology_info")

  expect_s3_class(topology, "irtsimrel_root_topology")
  expect_gt(nrow(topology$extrema), 0L)
  expect_gt(topology$rho_range[["upper"]], max(curve$rho_tilde))
})
