\name{merya-package}
\alias{merya-package}
\alias{merya}
\docType{package}
\title{Fast BCa Bootstrap Confidence Intervals and CI-Inversion P-Values}
\description{
\pkg{merya} provides drop-in replacements for the most commonly used
R statistical tests and models -- \code{\link[stats]{t.test}},
\code{\link[stats]{cor.test}}, \code{\link[stats]{lm}}, and
\code{\link[stats]{glm}} -- that report bias-corrected and accelerated
(BCa) bootstrap confidence intervals and confidence-interval-inversion
p-values, instead of relying on classical asymptotic distributions. All
bootstrap-based functions default to \code{R = 10000} replicates and to
a fixed \code{seed = 123} for reproducibility (set \code{seed = NULL}
to use whatever random state is currently active instead). The entire
package depends only on base R (\code{stats} and \code{utils}), with no
other package dependencies.
}
\details{
The exported functions mirror the argument lists and return values of
their original counterparts as closely as possible:
\itemize{
  \item \code{\link{boot.t.test}} -- like \code{\link[stats]{t.test}},
    returns an \code{"htest"} object. Performs one-sample, paired, and
    two independent samples t-tests.
  \item \code{\link{boot.cor.test}} -- like \code{\link[stats]{cor.test}},
    returns an \code{"htest"} object.
  \item \code{\link{boot.lm}} -- like \code{\link[stats]{lm}}, returns an
    \code{"lm"} object (subclassed \code{"boot.lm"}) with a bootstrap
    inference table available through \code{summary()}. By default
    reports the slope coefficients; can instead report, per predictor,
    the partial correlation, partial eta-squared, or eta-squared (see
    \code{effect} in \code{\link{boot.lm}}).
  \item \code{\link{boot.glm}} -- like \code{\link[stats]{glm}}, returns
    a \code{"glm"} object (subclassed \code{"boot.glm"}) with a bootstrap
    inference table available through \code{summary()}. For a binomial
    model, can instead report the marginal risk ratio or risk difference
    for a given binary predictor via g-computation (see \code{effect}
    and \code{exposure} in \code{\link{boot.glm}}).
  \item \code{\link{boot.ci}} -- a generic that takes an object already
    produced by \code{t.test()}, \code{cor.test()}, \code{lm()}, or
    \code{glm()} and returns it with BCa/CI-inversion inference
    substituted in, using the exact same internal machinery as the
    matching function above.
}

The bootstrap inference tables shown by \code{summary()} for
\code{\link{boot.lm}} and \code{\link{boot.glm}} report \code{Estimate},
\code{CI lower}, \code{CI upper}, and \code{Pr(>|z|)} (no separate
bootstrap standard error column).

Every function is written for speed:
\itemize{
  \item resampling indices are drawn once as a single \code{n x R}
    matrix rather than in an R-length loop;
  \item statistics that admit a closed form (means, mean differences,
    Pearson/Spearman correlations, OLS/GLM coefficients, partial and
    semi-partial correlations, g-computed risk ratios/differences) are
    computed with vectorized matrix/columnwise operations, from each
    replicate's own single refit -- no separate reduced-model refits are
    needed for the effect-size or g-computation options;
  \item the BCa acceleration constant is obtained from closed-form or
    one-step (infinitesimal jackknife) leave-one-out formulas, avoiding
    n full refits;
  \item the CI-inversion p-value is obtained by analytically inverting
    the BCa transformation, with no bisection or numerical root-finding.
}
}
\author{
Your Name <you@example.com>
}
\seealso{
\code{\link{boot.t.test}}, \code{\link{boot.cor.test}},
\code{\link{boot.lm}}, \code{\link{boot.glm}}, \code{\link{boot.ci}}
}
\keyword{package}
