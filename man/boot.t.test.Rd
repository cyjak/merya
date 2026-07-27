\name{boot.t.test}
\alias{boot.t.test}
\title{Student's t-Test with BCa Bootstrap Confidence Interval}
\description{
Performs one-sample, paired, and two independent samples t-tests, using the
same calling convention and output structure as \code{\link[stats]{t.test}},
but the confidence interval is obtained from a bias-corrected and
accelerated (BCa) bootstrap of the relevant mean / mean-difference, and
the p-value is obtained by inverting that same BCa confidence interval at
the hypothesised value \code{mu} (confidence-interval-inversion testing),
rather than from the classical t-distribution.
}
\usage{
boot.t.test(x, y = NULL,
            alternative = c("two.sided", "less", "greater"),
            mu = 0, paired = FALSE, var.equal = FALSE,
            conf.level = 0.95, R = 10000, seed = 123, \dots)
}
\arguments{
  \item{x}{a (non-empty) numeric vector of data values.}
  \item{y}{an optional (non-empty) numeric vector of data values.}
  \item{alternative}{character string, one of \code{"two.sided"} (default),
    \code{"greater"} or \code{"less"}.}
  \item{mu}{a number indicating the true value of the mean (or difference
    in means if performing a two-sample test).}
  \item{paired}{logical; if \code{TRUE}, a paired test is performed.}
  \item{var.equal}{currently ignored; kept for interface compatibility
    with \code{\link[stats]{t.test}}. The bootstrap does not assume equal
    variances.}
  \item{conf.level}{confidence level of the returned interval.}
  \item{R}{number of bootstrap replicates (default \code{10000}).}
  \item{seed}{integer seed used to make the bootstrap resampling
    reproducible; set via \code{\link[base]{set.seed}} at the start of
    the function. Defaults to \code{123}; pass \code{NULL} to use
    whatever random state is currently active (no reseeding), or any
    other integer for a different reproducible draw.}
  \item{\dots}{further arguments (ignored; present for compatibility).}
}
\details{
The bootstrap statistic is the sample mean (one-sample / paired case) or
the difference in sample means (two-sample case). Resampling is fully
vectorized: a single \code{n x R} index matrix is drawn and all replicate
means are computed with \code{colMeans()}. The BCa acceleration constant
is obtained from a closed-form jackknife on the mean, so no
per-replicate loop is required. The p-value is obtained by analytically
inverting the BCa transformation at \code{mu} -- no bisection or search
is performed, which keeps the procedure fast even for large \code{R}.
}
\value{
A list with class \code{"htest"} containing the same components as
\code{\link[stats]{t.test}}: \code{statistic}, \code{parameter},
\code{p.value}, \code{conf.int}, \code{estimate}, \code{null.value},
\code{stderr}, \code{alternative}, \code{method}, \code{data.name}, plus
\code{R} (the number of bootstrap replicates used).
}
\seealso{
\code{\link[stats]{t.test}}, \code{\link{boot.cor.test}},
\code{\link{boot.lm}}, \code{\link{boot.glm}}
}
\examples{
x <- rnorm(30, mean = 1)
boot.t.test(x)

y <- rnorm(30, mean = 0)
boot.t.test(x, y)

boot.t.test(x, y, paired = TRUE)
}
