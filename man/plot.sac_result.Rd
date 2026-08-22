% Please edit documentation in R/spc_calibrate.R
\name{plot.sac_result}
\alias{plot.sac_result}
\title{Plot an SAC Convergence Trajectory}
\usage{
\method{plot}{sac_result}(x, type = c("both", "trajectory", "c", "rho"), ...)
}
\arguments{
\item{x}{A \code{sac_result} object from \code{\link[=sac_calibrate]{sac_calibrate()}}.}

\item{type}{Character. One of \code{"both"}, \code{"trajectory"}, \code{"c"}, or
\code{"rho"}. \code{"trajectory"} is an alias for \code{"c"}.}

\item{...}{Additional arguments (currently unused).}
}
\value{
A \code{ggplot} object if \pkg{ggplot2} is available for
\code{type = "c"} or \code{type = "rho"}. For \code{type = "both"}, a patchwork
object is returned when \pkg{patchwork} is available; otherwise both
plots are printed and returned invisibly in a list. Without \pkg{ggplot2},
base R graphics are drawn and \code{NULL} is returned invisibly.
}
\description{
Visualizes the aligned Robbins--Monro scale and reliability trajectories
stored in a \code{sac_result}.
}
