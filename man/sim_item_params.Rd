% Please edit documentation in R/sim_item_params.R
\name{sim_item_params}
\alias{sim_item_params}
\alias{print.item_params}
\alias{summary.item_params}
\title{Simulate Item Parameters for IRT Studies}
\usage{
sim_item_params(
  n_items,
  model = c("rasch", "2pl"),
  source = c("parametric", "irw", "hierarchical", "custom"),
  method = c("copula", "conditional", "independent"),
  n_forms = 1L,
  difficulty_params = list(),
  discrimination_params = list(),
  hierarchical_params = list(),
  custom_params = list(),
  scale = 1,
  center_difficulties = TRUE,
  seed = NULL
)

\method{print}{item_params}(x, digits = 4, ...)

\method{summary}{item_params}(object, ...)
}
\arguments{
\item{n_items}{Integer. Number of items to generate per form.}

\item{model}{Character. The data-generating model: "rasch" or "2pl".}

\item{source}{Character. Source for generating difficulties:
\describe{
\item{\code{"parametric"}}{Generate from a parametric difficulty distribution (default).}
\item{\code{"irw"}}{Use an IRW difficulty pool when the external IRW package is installed.}
\item{\code{"hierarchical"}}{Joint MVN for both parameters (Glas & van der Linden)}
\item{\code{"custom"}}{User-supplied parameters or function}
}}

\item{method}{Character. Method for generating discriminations (when model = "2pl"):
\describe{
\item{\code{"copula"}}{Gaussian copula - preserves marginals exactly (RECOMMENDED)}
\item{\code{"conditional"}}{Conditional normal regression on difficulty}
\item{\code{"independent"}}{Independent generation (no correlation)}
}}

\item{n_forms}{Integer. Number of test forms to generate. Default is 1.
When > 1, returns a data frame with form_id column.}

\item{difficulty_params}{List. Parameters for difficulty generation:
\describe{
\item{For \code{source = "irw"}:}{\code{pool} - optional difficulty pool data frame.}
\item{For \code{source = "parametric"}:}{\code{mu} - finite mean
(default 0), \code{sigma} - positive finite SD (default 1), and
\code{distribution} - \code{"normal"} or \code{"uniform"} (default
\code{"normal"}).}
}}

\item{discrimination_params}{List. Parameters for discrimination generation:
\describe{
\item{\code{mu_log}}{Finite mean of log-discrimination (default: 0).}
\item{\code{sigma_log}}{Positive finite SD of log-discrimination (default: 0.3).}
\item{\code{rho}}{Finite target correlation between \eqn{\beta} and
\eqn{\log(\lambda)} in the closed interval from -1 to 1
(default: -0.3).}
}}

\item{hierarchical_params}{List. For source = "hierarchical":
\describe{
\item{\code{mu}}{Finite 2-vector: means of \eqn{(\log\lambda, \beta)}.}
\item{\code{tau}}{Positive finite 2-vector of SDs.}
\item{\code{rho}}{Finite correlation in (-1, 1).}
}}

\item{custom_params}{List. For source = "custom":
\describe{
\item{\code{beta}}{Finite numeric vector of length \code{n_items}, or
function returning one.}
\item{\code{lambda}}{Positive finite numeric vector of length
\code{n_items}, or function returning one; required for
\code{model = "2pl"}.}
}}

\item{scale}{Numeric. Global discrimination scaling factor for reliability targeting.
Final discriminations are \eqn{\lambda_i^* = c \cdot \lambda_i}. Default is 1.}

\item{center_difficulties}{Logical. If TRUE, center difficulties to sum to zero
for identification. Default is TRUE.}

\item{seed}{Integer. Random seed for reproducibility.}

\item{x}{An object of class \code{"item_params"}.}

\item{digits}{Integer. Number of decimal places for printing.}

\item{...}{Additional arguments passed to or from other methods.}

\item{object}{An object of class \code{"item_params"}.}
}
\value{
An object of class \code{"item_params"} containing:
\describe{
\item{\code{data}}{Data frame with columns: form_id, item_id, beta, lambda, lambda_unscaled}
\item{\code{model}}{Model type used}
\item{\code{source}}{Source used for generation}
\item{\code{method}}{Method used for discrimination generation}
\item{\code{n_items}}{Number of items per form}
\item{\code{n_forms}}{Number of forms generated}
\item{\code{scale}}{Scale factor applied}
\item{\code{centered}}{Whether difficulties were centered}
\item{\code{params}}{Parameters used for generation}
\item{\code{achieved}}{Achieved statistics (correlations, moments)}
}

The input object, invisibly.

An object of class \code{"summary.item_params"} containing key
parameter summaries.
}
\description{
\code{sim_item_params()} generates item parameters (difficulty \eqn{\beta} and
discrimination \eqn{\lambda}) for Item Response Theory (IRT) simulation studies.
It provides parametric, hierarchical, custom, and optional Item Response
Warehouse (IRW) sources, plus multiple methods for generating correlated
discriminations.

The function is designed with five key principles:
\enumerate{
\item \strong{Portable defaults:} Parametric item generation works without
external item-pool data.
\item \strong{Empirical integration:} Optional IRW integration is available
when the Item Response Warehouse package is installed.
\item \strong{Correlated parameters:} Support for the empirically observed negative
correlation between difficulty and discrimination (Sweeney et al., 2022).
\item \strong{Marginal preservation:} Copula method preserves exact marginal
distributions while achieving target correlation.
\item \strong{Reliability targeting:} Scale factor for subsequent calibration.
}
}
\details{
\subsection{Why the Copula Method is Recommended}{

When difficulties come from realistic or otherwise non-normal marginal
distributions, the conditional normal method can distort the achieved
correlation because it assumes linearity. The Gaussian copula method:

\enumerate{
\item Transforms difficulties to uniform scale via empirical CDF
\item Generates correlated uniforms through Gaussian copula
\item Transforms back to desired marginals (log-normal for discrimination)
}

This guarantees:
\itemize{
\item Exact preservation of the sampled difficulty marginal
\item Exact log-normal marginal for discriminations
\item Spearman correlation \eqn{\approx \rho} (rank-based, robust to non-normality)
}
}

\subsection{Connection to Reliability-Targeted Framework}{

The \code{scale} parameter implements "separation of structure and scale":
\itemize{
\item \strong{Structure}: Item characteristics and optional difficulty-discrimination correlation
\item \strong{Scale}: Global factor \eqn{c} calibrated for target reliability
}
}
}
\examples{
# Example 1: Rasch with parametric difficulties
items1 <- sim_item_params(n_items = 25, model = "rasch",
                          source = "parametric", seed = 42)

# Example 2: 2PL with copula method (recommended)
items2 <- sim_item_params(
  n_items = 30, model = "2pl", source = "parametric",
  method = "copula",
  discrimination_params = list(rho = -0.3),
  seed = 42
)

# Example 3: Hierarchical 2PL
items3 <- sim_item_params(
  n_items = 25, model = "2pl", source = "hierarchical",
  hierarchical_params = list(mu = c(0, 0), tau = c(0.25, 1), rho = -0.3),
  seed = 42
)

\dontrun{
# Example 4: Using IRW difficulty pool (requires irw package)
if (requireNamespace("irw", quietly = TRUE)) {
  items4 <- sim_item_params(n_items = 25, model = "rasch", source = "irw")

  # Example 5: Multiple forms with IRW
  items5 <- sim_item_params(
    n_items = 20, model = "2pl", n_forms = 5,
    source = "irw", method = "copula"
  )
}
}

}
\references{
Glas, C. A. W., & van der Linden, W. J. (2003). Computerized adaptive testing
with item cloning. \emph{Applied Psychological Measurement, 27}(4), 247-261.

Sweeney, S. M., et al. (2022). An investigation of the nature and consequence
of the relationship between IRT difficulty and discrimination.
\emph{EM:IP, 41}(4), 50-67.

Zhang, L., et al. (2025). Realistic simulation of item difficulties. \emph{PsyArXiv}.
}
