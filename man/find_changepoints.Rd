% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_changepoints.R
\name{find_changepoints}
\alias{find_changepoints}
\title{Find changepoints.}
\usage{
find_changepoints(
  y,
  method,
  threshold = NULL,
  maxiter = NULL,
  num_rand_ints = NULL,
  rand_ints = NULL,
  seeded = FALSE,
  decay = NULL
)
}
\arguments{
\item{y}{Numeric vector of data.}

\item{method}{Character string: \code{"bs"} for binary segmentation; \code{"wbs"} for wild binary segmentation;
\code{"not"} for narrowest over threshold.}

\item{threshold}{Numeric; changepoint detection threshold for CUSUM statistic.}

\item{maxiter}{Integer; maximum number of changepoints to detect.}

\item{num_rand_ints}{Integer; number of random intervals. Ignored if \code{rand_ints} is specified, or if
\code{method = "bs"}.}

\item{rand_ints}{Matrix containing random intervals for wild binary segmentation or narrowest over threshold.
Ignored if \code{method = "bs"}. Optional.}

\item{seeded}{Logical; if \code{TRUE} and \code{method = "wbs"}, then seeded binary segmentation is implemented.}

\item{decay}{Decay parameter for seeded binary segmentation; only used if \code{method = "wbs"} and \code{seeded = TRUE}.}
}
\value{
A list:
\itemize{
\item \code{results} Dataframe containing results
\item \code{changepoints} Vector of changepoints detected
\item \code{rand_ints} (Not if \code{method = "bs"}) Matrix containing random intervals used
\item \code{threshold} Value of \code{threshold}
\item \code{maxiter} Value of \code{maxiter}
}
}
\description{
Detect changepoints in the mean using binary segmentation, wild binary segmentation, or narrowest over threshold.
}
\examples{
set.seed(100)
y <- rnorm(100) + c(rep(1,45), rep(-1,10), rep(1,45))
results_bs <- find_cps(y, "bs", threshold=4)
print(results_bs$results)

results_wbs <- find_cps(y, "wbs", threshold=4, num_rand_ints=100)
print(results_wbs$results)

results_not <- find_cps(y, "not", threshold=4, num_rand_ints=100)
print(results_not$results)

}
