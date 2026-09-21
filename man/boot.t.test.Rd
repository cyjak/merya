\name{boot.t.test}
\alias{boot.t.test}
\title{Student's t-Test with BCa Bootstrap Confidence Intervals}
\description{
Performs one-sample, paired, and two independent samples t-tests, using
the same arguments as \code{\link[stats]{t.test}} and returning an
object of the same \code{"htest"} class, but with a bootstrap
confidence interval and CI-inversion p-value in place of the classical
asymptotic ones.

Internally, every case is expressed as a linear regression and handed
to the exact same engines \code{\link{boot.lm}} uses (an intercept-only
regression of \code{x} for the one-sample case, of the paired
differences \code{x - y} for the paired case, or of the combined sample
on a 0/1 group indicator for the two-sample case, testing that
indicator's coefficient -- which is algebraically identical to the
difference in means). This means \code{boot.t.test(x, y)} and
\code{boot.lm(c(x, y) ~ group)}'s group coefficient, for the same
\code{boot.method}/\code{wild.dist}/\code{ci.type}/\code{seed}, are not
just conceptually comparable but numerically identical.
}
\usage{
boot.t.test(x, y = NULL,
            alternative = c("two.sided", "less", "greater"),
            mu = 0, paired = FALSE, var.equal = FALSE,
            conf.level = 0.95, R = 10000,
            boot.method = c("wild", "case"),
            wild.dist = c("rademacher", "mammen", "normal"),
            ci.type = c("bca", "percentile"),
            seed = 123, \dots)
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
    with \code{\link[stats]{t.test}}. Neither bootstrap method assumes
    equal variances in the classical (pooled-variance) sense; see
    Details.}
  \item{conf.level}{confidence level of the returned interval.}
  \item{R}{number of bootstrap replicates (default \code{10000}).}
  \item{boot.method}{resampling method: \code{"wild"} (default) or
    \code{"case"}; see \code{\link{boot.lm}}, whose engines this
    function reuses directly (see Details).}
  \item{wild.dist}{distribution of the wild-bootstrap multipliers (only
    used when \code{boot.method = "wild"}); see \code{\link{boot.lm}}.}
  \item{ci.type}{how the confidence interval (and CI-inversion p-value)
    is obtained from the bootstrap distribution: \code{"bca"} (default)
    or \code{"percentile"}; see \code{\link{boot.lm}}.}
  \item{seed}{integer seed used to make the bootstrap resampling
    reproducible; set via \code{\link[base]{set.seed}} at the start of
    the function. Defaults to \code{123}; pass \code{NULL} to use
    whatever random state is currently active (no reseeding), or any
    other integer for a different reproducible draw.}
  \item{\dots}{further arguments (ignored; present for compatibility).}
}
\value{
A list with class \code{"htest"} containing the same components as
\code{\link[stats]{t.test}} (\code{statistic}, \code{parameter},
\code{p.value}, \code{conf.int}, \code{estimate}, \code{null.value},
\code{stderr}, \code{alternative}, \code{method}, \code{data.name}),
plus \code{R}, \code{boot.method}, and \code{ci.type}.
}
\details{
The bootstrap statistic is the sample mean (one-sample / paired case) or
the difference in sample means (two-sample case), obtained as the
relevant coefficient of an equivalent OLS regression via
\code{\link{boot.lm}}'s own internal engines -- the bootstrap
distribution, and a closed-form leave-one-out jackknife for the BCa
acceleration constant, both closed-form/vectorized with no explicit
loop over observations or (for \code{"wild"}) replicates. Because this
reduces the two-sample case to a *combined*, unstratified regression on
a group indicator, both \code{boot.method}s here resample/perturb the
pooled sample rather than resampling each group separately at its own
fixed size; this is what guarantees consistency with
\code{\link{boot.lm}} for the same data, and is standard practice for
the wild bootstrap in particular (which has no natural "resample within
group" analogue, since it does not resample rows at all). The p-value
is obtained by analytically inverting the BCa (or, for \code{ci.type =
"percentile"}, empirical-quantile) transformation at \code{mu}, in the
same way as every other function in this package.
}
\examples{
x <- rnorm(30, mean = 1)
boot.t.test(x)

y <- rnorm(30, mean = 0)
boot.t.test(x, y)

boot.t.test(x, y, paired = TRUE)

## the classical case-resampling bootstrap, instead of the wild-
## bootstrap default
boot.t.test(x, y, boot.method = "case")

## consistency with boot.lm(): the "group" coefficient below matches
## boot.t.test(x, y)'s estimate/CI/p-value exactly (same seed, method)
d <- data.frame(z = c(x, y), group = c(rep(1, length(x)), rep(0, length(y))))
boot.lm(z ~ group, data = d)
}
\seealso{
\code{\link[stats]{t.test}}, \code{\link{boot.lm}}
}
