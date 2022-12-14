% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_interval_bs.R
\name{calculate_interval_bs}
\alias{calculate_interval_bs}
\title{Calculate interval for binary segmentation}
\usage{
calculate_interval_bs(
  y,
  nu,
  b,
  d,
  nu2 = NULL,
  nuTy = NULL,
  threshold = NULL,
  n.cp = NULL
)
}
\arguments{
\item{y}{Numeric vector of data.}

\item{nu}{Numeric vector.}

\item{b}{Vector of changepoints detected by binary segmentation algorithm.}

\item{d}{Directions of changepoints detected by binary segmentation algorithm. Vector whose entries are all equal to \code{1}
or \code{-1}.}

\item{nu2}{Value of \eqn{||\nu||_2^2}; optional.}

\item{nuTy}{Value of \eqn{\nu^T y}; optional.}

\item{threshold}{Threshold used in binary segmentation algorithm.}

\item{n.cp}{Maximum number of changepoints to detect of binary segmentation algorithm.}
}
\value{
A 2-dimensional vector.
}
\description{
Find values of phi which satisfy the required inequalities so that applying binary segmentation to \eqn{y'(\phi)}
returns \code{(b, d)}.
Option to only consider part of \code{(b, d)}, e.g. if we want \code{b[1]} to be in the set of detected changepoints but are not concerned
about other changepoints.
}
\details{
Used inside \code{calculate_S} if \code{method = "bs"}.
}
\examples{
set.seed(100)
y <- rnorm(100) + c(rep(1,50), rep(-1,50))
results <- binary_segmentation(y, threshold=4)
b <- results$results$b[ results$results$cp==1 ]
d <- results$results$d[ results$results$cp==1 ]
h <- 10
nu <- c(rep(0, b[1]-h), rep(1/h, h), rep(-1/h, h), rep(0, length(y)-b[1]-h))
calculate_interval_bs(y, nu, b, d, threshold=4)

}
