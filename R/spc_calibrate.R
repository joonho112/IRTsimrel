# =============================================================================
# spc_calibrate.R
# =============================================================================
# Backward-compatibility entry point plus the SAC plotting method.  The single
# production definition of sac_calibrate() lives in stochastic_calibration.R;
# the remaining SAC S3 methods live in stochastic_methods.R.
# =============================================================================


#' Deprecated Alias for SAC Calibration
#'
#' `spc_calibrate()` is the deprecated v0.1.x name for
#' [sac_calibrate()]. It forwards all arguments unchanged and returns the
#' resulting `sac_result` after issuing R's standard deprecation warning.
#'
#' @param ... Arguments passed to [sac_calibrate()].
#'
#' @return The result of [sac_calibrate()].
#'
#' @seealso [sac_calibrate()]
#'
#' @export
spc_calibrate <- function(...) {
  .Deprecated("sac_calibrate")
  sac_calibrate(...)
}


#' Plot an SAC Convergence Trajectory
#'
#' Visualizes the aligned Robbins--Monro scale and reliability trajectories
#' stored in a `sac_result`.
#'
#' @param x A `sac_result` object from [sac_calibrate()].
#' @param type Character. One of `"both"`, `"trajectory"`, `"c"`, or
#'   `"rho"`. `"trajectory"` is an alias for `"c"`.
#' @param ... Additional arguments (currently unused).
#'
#' @return A `ggplot` object if \pkg{ggplot2} is available for
#'   `type = "c"` or `type = "rho"`. For `type = "both"`, a patchwork
#'   object is returned when \pkg{patchwork} is available; otherwise both
#'   plots are printed and returned invisibly in a list. Without \pkg{ggplot2},
#'   base R graphics are drawn and `NULL` is returned invisibly.
#'
#' @export
plot.sac_result <- function(x,
                            type = c("both", "trajectory", "c", "rho"),
                            ...) {
  .irtsimrel_validate_sac_result_object(x, "x", "SAC plot")

  type <- match.arg(type)
  if (type == "trajectory") type <- "c"

  n_iter <- x$n_iter
  burn_in <- x$burn_in

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    if (type == "both") {
      oldpar <- par(mfrow = c(2, 1))
      on.exit(par(oldpar))
    }

    if (type %in% c("both", "c")) {
      plot(
        seq_len(n_iter), x$trajectory,
        type = "l", col = "steelblue",
        xlab = "Iteration", ylab = "c",
        main = "SAC Convergence Trajectory"
      )
      abline(h = x$c_star, col = "red", lty = 2, lwd = 2)
      abline(v = burn_in, col = "gray", lty = 3)
      legend(
        "topright",
        legend = c("Trajectory", "c* (average)", "Burn-in"),
        col = c("steelblue", "red", "gray"),
        lty = c(1, 2, 3),
        lwd = c(1, 2, 1)
      )
    }

    if (type %in% c("both", "rho")) {
      plot(
        seq_len(n_iter), x$rho_trajectory,
        type = "l", col = "darkorange",
        xlab = "Iteration", ylab = expression(hat(rho)),
        main = "Reliability Estimates"
      )
      abline(h = x$target_rho, col = "red", lty = 2, lwd = 2)
      abline(v = burn_in, col = "gray", lty = 3)
    }

    return(invisible(NULL))
  }

  df <- data.frame(
    iteration = seq_len(n_iter),
    c = x$trajectory,
    rho = x$rho_trajectory,
    phase = factor(
      ifelse(seq_len(n_iter) <= burn_in, "Burn-in", "Averaging"),
      levels = c("Burn-in", "Averaging")
    )
  )

  p1 <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = iteration, y = c, color = phase)
  ) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_hline(
      yintercept = x$c_star,
      color = "red",
      linetype = "dashed",
      linewidth = 1
    ) +
    ggplot2::geom_vline(
      xintercept = burn_in,
      color = "gray50",
      linetype = "dotted"
    ) +
    ggplot2::scale_color_manual(
      values = c("Burn-in" = "steelblue", "Averaging" = "darkblue")
    ) +
    ggplot2::labs(
      title = "SAC Convergence Trajectory",
      subtitle = sprintf(
        "Target rho* = %.3f | c* = %.4f | Init: %s",
        x$target_rho, x$c_star, x$init_method
      ),
      x = "Iteration", y = "Scaling Factor c", color = "Phase"
    ) +
    ggplot2::annotate(
      "text",
      x = n_iter * 0.95,
      y = x$c_star,
      label = sprintf("c* = %.3f", x$c_star),
      vjust = -0.5,
      color = "red",
      size = 3.5
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  p2 <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = iteration, y = rho, color = phase)
  ) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_hline(
      yintercept = x$target_rho,
      color = "red",
      linetype = "dashed",
      linewidth = 1
    ) +
    ggplot2::geom_vline(
      xintercept = burn_in,
      color = "gray50",
      linetype = "dotted"
    ) +
    ggplot2::scale_color_manual(
      values = c("Burn-in" = "darkorange", "Averaging" = "darkred")
    ) +
    ggplot2::labs(
      title = "Reliability Estimates Across Iterations",
      subtitle = sprintf(
        "Target rho* = %.3f | Achieved = %.3f",
        x$target_rho, x$achieved_rho
      ),
      x = "Iteration", y = expression(hat(rho)), color = "Phase"
    ) +
    ggplot2::annotate(
      "text",
      x = n_iter * 0.95,
      y = x$target_rho,
      label = sprintf("rho* = %.3f", x$target_rho),
      vjust = -0.5,
      color = "red",
      size = 3.5
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  if (type == "c") {
    p1
  } else if (type == "rho") {
    p2
  } else if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(p1, p2, ncol = 1)
  } else {
    print(p1)
    print(p2)
    invisible(list(p1 = p1, p2 = p2))
  }
}
