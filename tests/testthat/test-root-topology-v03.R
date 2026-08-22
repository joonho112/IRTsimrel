# =============================================================================
# Shared root-topology engine (v0.3)
# =============================================================================

test_that("monotone topology has one increasing selectable root", {
  topology <- .irtsimrel_scan_topology(
    function(c) plogis(log(c)),
    c_bounds = c(0.01, 100),
    target = 0.5
  )

  expect_s3_class(topology, "irtsimrel_root_topology")
  expect_identical(topology$topology_status, "resolved")
  expect_identical(topology$target_status, "feasible")
  expect_equal(topology$root_count, 1L)
  expect_equal(topology$admissible_root_count, 1L)
  expect_equal(topology$roots$c, 1, tolerance = 1e-7)
  expect_identical(topology$roots$direction, "increasing")
  expect_named(
    topology$best_achievable,
    c("c", "rho", "residual", "absolute_error", "location")
  )

  selected <- .irtsimrel_select_root(topology)
  expect_equal(selected$c, 1, tolerance = 1e-7)
  expect_true(selected$topology$roots$selected[[1L]])
})

test_that("single-hump topology enumerates both crossing directions", {
  topology <- .irtsimrel_scan_topology(
    function(c) exp(-(log(c))^2),
    c_bounds = c(0.01, 100),
    target = 0.5
  )

  expected <- exp(c(-sqrt(log(2)), sqrt(log(2))))
  expect_equal(topology$root_count, 2L)
  expect_equal(topology$roots$c, expected, tolerance = 1e-5)
  expect_identical(
    topology$roots$direction,
    c("increasing", "decreasing")
  )
  expect_equal(.irtsimrel_select_root(topology)$c, expected[[1L]],
               tolerance = 1e-5)
  expect_equal(
    .irtsimrel_select_root(topology, "highest_any")$c,
    expected[[2L]],
    tolerance = 1e-5
  )
})

test_that("boundary, tangent, and plateau candidates remain distinct", {
  boundary <- .irtsimrel_scan_topology(
    function(c) plogis(log(c)), c(0.1, 10), target = plogis(log(0.1))
  )
  expect_identical(boundary$roots$kind[[1L]], "boundary")
  expect_equal(.irtsimrel_select_root(boundary)$c, 0.1, tolerance = 1e-8)

  tangent <- .irtsimrel_scan_topology(
    function(c) 1 - log(c)^2, c(0.1, 10), target = 1
  )
  expect_equal(tangent$root_count, 1L)
  expect_identical(tangent$roots$kind, "tangent")
  expect_identical(tangent$roots$direction, "flat")
  expect_error(
    .irtsimrel_select_root(tangent),
    "No admissible root",
    class = "irtsimrel_selection_unavailable"
  )
  expect_equal(
    .irtsimrel_select_root(
      tangent, root_policy = "lowest_any", allow_tangent = TRUE
    )$c,
    1,
    tolerance = 1e-6
  )

  plateau <- .irtsimrel_scan_topology(
    function(c) 0.6 - pmax(abs(log(c)) - 0.25, 0)^2,
    c(0.1, 10),
    target = 0.6
  )
  expect_equal(plateau$root_count, 1L)
  expect_identical(plateau$roots$kind, "plateau")
  expect_identical(plateau$roots$direction, "flat")
  expect_error(
    .irtsimrel_select_root(plateau),
    "No admissible root",
    class = "irtsimrel_selection_unavailable"
  )
  expect_equal(
    .irtsimrel_select_root(
      plateau, root_policy = "lowest_any", allow_plateau = TRUE
    )$c,
    1,
    tolerance = 1e-6
  )
})

test_that("infeasible and unresolved statuses are orthogonal", {
  above <- .irtsimrel_scan_topology(
    function(c) plogis(log(c)), c(0.1, 10), target = 0.99
  )
  expect_identical(above$topology_status, "resolved")
  expect_identical(above$target_status, "infeasible_above_range")
  expect_equal(above$root_count, 0L)
  expect_true(is.finite(above$best_achievable$c))

  nonfinite <- .irtsimrel_scan_topology(
    function(c) {
      if (c > 0.8 && c < 1.2) return(NA_real_)
      plogis(log(c))
    },
    c(0.1, 10),
    target = 0.5
  )
  expect_identical(nonfinite$topology_status, "nonfinite")
  expect_identical(nonfinite$target_status, "uncertain")

  all_nonfinite <- .irtsimrel_scan_topology(
    function(c) NA_real_, c(0.1, 10), target = 0.5
  )
  expect_identical(all_nonfinite$topology_status, "nonfinite")
  expect_true(is.na(all_nonfinite$best_achievable$c))
  expect_identical(all_nonfinite$best_achievable$location, "unavailable")

  limited <- .irtsimrel_scan_topology(
    function(c) plogis(log(c)),
    c(0.1, 10),
    target = 0.5,
    controls = list(min_grid = 65L, max_evals = 10L)
  )
  expect_identical(limited$topology_status, "resolution_limited")
  expect_identical(limited$target_status, "uncertain")
  expect_error(
    .irtsimrel_select_root(limited),
    "status",
    class = "irtsimrel_topology_unresolved"
  )
})

test_that("capped curvature refinement detects a narrow hidden bump", {
  initial_spacing <- diff(log(c(0.1, 10))) / 64
  center <- initial_spacing / 2
  narrow_bump <- function(c) {
    0.2 + 0.5 * exp(-((log(c) - center) / 0.012)^2)
  }

  resolved <- .irtsimrel_scan_topology(
    narrow_bump, c(0.1, 10), target = 0.4
  )
  expect_identical(resolved$topology_status, "resolved")
  expect_equal(resolved$root_count, 2L)
  expect_identical(
    resolved$roots$direction,
    c("increasing", "decreasing")
  )

  limited <- .irtsimrel_scan_topology(
    narrow_bump,
    c(0.1, 10),
    target = 0.4,
    controls = list(max_evals = 65L)
  )
  expect_identical(limited$topology_status, "resolution_limited")
  expect_identical(limited$target_status, "uncertain")
  expect_equal(limited$root_count, 0L)
})

test_that("empirical 2PL fixture exposes all four roots", {
  theta <- c(
    -2.8945639829, -2.1638196461, -1.1142421522, -0.4808900280,
    -0.4195555953, 0.0231953129, 1.6704952657
  )
  beta <- c(-3.4556089868, -1.7336287633, 0.8669555536)
  lambda <- c(4.8618527891, 1.4673837406, 0.0537872782)
  objective <- function(c) compute_rho_tilde(c, theta, beta, lambda)

  topology <- .irtsimrel_scan_topology(
    objective, c(0.005, 500), target = 0.18
  )

  expect_identical(topology$topology_status, "resolved")
  expect_equal(topology$root_count, 4L)
  expect_equal(topology$admissible_root_count, 2L)
  expect_equal(
    topology$roots$c,
    c(0.1890451, 11.592282, 13.402881, 102.05913),
    tolerance = 2e-5
  )
  expect_identical(
    topology$roots$direction,
    c("increasing", "decreasing", "increasing", "decreasing")
  )
  expect_equal(.irtsimrel_select_root(topology)$c, 0.1890451,
               tolerance = 2e-5)
})

test_that("empirical 3PL fixture exposes all four roots", {
  theta <- c(
    -2.8945639829, -2.1638196461, -1.1142421522, -0.4808900280,
    -0.4195555953, 0.0231953129, 1.6704952657
  )
  beta <- c(-3.4556089868, -1.7336287633, 0.8669555536)
  lambda <- c(4.8618527891, 1.4673837406, 0.0537872782)
  guessing <- c(0.20, 0.25, 0.10)
  objective <- function(c) {
    compute_rho_tilde(c, theta, beta, lambda, guessing = guessing)
  }

  topology <- .irtsimrel_scan_topology(
    objective, c(0.005, 500), target = 0.09
  )

  expect_equal(topology$root_count, 4L)
  expect_equal(topology$admissible_root_count, 2L)
  expect_equal(
    topology$roots$c,
    c(0.1256152, 7.4568044, 10.2013650, 106.04184),
    tolerance = 2e-5
  )
  expect_equal(.irtsimrel_select_root(topology)$c, 0.1256152,
               tolerance = 2e-5)
})

test_that("topology controls and selection policies validate inputs", {
  expect_error(
    .irtsimrel_scan_topology(identity, c(0, 1)),
    "c_bounds"
  )
  expect_error(
    .irtsimrel_scan_topology(identity, c(0.1, 1), controls = list(foo = 1)),
    "Unknown"
  )

  topology <- .irtsimrel_scan_topology(
    function(c) exp(-(log(c))^2), c(0.1, 10), target = 0.5
  )
  expect_error(
    .irtsimrel_select_root(topology, "nearest_increasing"),
    "reference"
  )
  expect_equal(
    .irtsimrel_select_root(
      topology, "nearest_increasing", reference = 0.5
    )$root_id,
    1L
  )
})
