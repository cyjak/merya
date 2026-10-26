\name{boot.cor.test}
\alias{boot.cor.test}
\title{Correlation Test with BCa Bootstrap Confidence Interval}
\description{
Test for association between paired samples, using the same calling
convention and output structure as \code{\link[stats]{cor.test}}, but
with a bias-corrected and accelerated (BCa) bootstrap confidence interval
for the correlation coefficient, and a p-value obtained by inverting that
interval at zero correlation (confidence-interval inversion testing).
}
\usage{
boot.cor.test(x, y, method = c("pearson", "kendall", "spearman"),
              alternative = c("two.sided", "less", "greater"),
              conf.level = 0.95, R = 10000,
              ci.type = c("bca", "percentile"),
              seed = 123, \dots)
}
\arguments{
  \item{x, y}{numeric vectors of the same length (at least 3 complete
    observations required).}
  \item{method}{character string: \code{"pearson"} (default),
    \code{"kendall"}, or \code{"spearman"}.}
  \item{alternative}{character string: \code{"two.sided"} (default),
    \code{"greater"} or \code{"less"}.}
  \item{conf.level}{confidence level for the returned interval.}
  \item{R}{number of bootstrap replicates (default \code{10000}).}
  \item{ci.type}{how the confidence interval (and, by CI inversion, the
    p-value) is obtained from the bootstrap distribution: \code{"bca"}
    (default) or \code{"percentile"}; see \code{\link{boot.lm}}.}
  \item{seed}{integer seed used to make the bootstrap resampling
    reproducible; set via \code{\link[base]{set.seed}} at the start of
    the function. Defaults to \code{123}; pass \code{NULL} to use
    whatever random state is currently active (no reseeding), or any
    other integer for a different reproducible draw.}
  \item{\dots}{further arguments (ignored; present for compatibility).}
}
\details{
For \code{method = "pearson"} and \code{"spearman"}, the bootstrap
replicates are computed in a fully vectorized way: a single \code{n x R}
resampling-index matrix is drawn, and all \code{R} replicate correlations
are obtained from columnwise sufficient statistics (sums, sums of
squares/products) rather than a per-replicate loop. \code{method =
"kendall"} does not admit the same closed-form vectorization and falls
back to one compiled \code{cor()} call per replicate (slower). The
acceleration constant of the BCa interval is obtained from a closed-form
leave-one-out (jackknife) formula for Pearson/Spearman, evaluated in
O(n).
}
\value{
A list with class \code{"htest"} with the same components as
\code{\link[stats]{cor.test}}: \code{statistic}, \code{parameter},
\code{p.value}, \code{estimate}, \code{null.value}, \code{alternative},
\code{method}, \code{data.name}, \code{conf.int}, plus \code{R}.
}
\seealso{
\code{\link[stats]{cor.test}}, \code{\link{boot.t.test}},
\code{\link{boot.lm}}, \code{\link{boot.glm}}
}
\examples{
x <- rnorm(40); y <- x * 0.5 + rnorm(40)
boot.cor.test(x, y)
boot.cor.test(x, y, method = "spearman")
}
