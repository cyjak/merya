\name{boot.glm}
\alias{boot.glm}
\alias{summary.boot.glm}
\alias{print.summary.boot.glm}
\title{Generalized Linear Regression with BCa Bootstrap Confidence Intervals}
\description{
Fits a generalized linear model exactly as \code{\link[stats]{glm}} does
(same formula interface, same fitted \code{"glm"} object), but
additionally performs a case-resampling bootstrap and attaches
bias-corrected and accelerated (BCa) confidence intervals and
CI-inversion p-values, retrievable via \code{summary()}. By default the
bootstrapped quantity is the model coefficients (as before); for a
binomial model, set \code{effect} to instead report a marginal
(population-averaged) risk ratio or risk difference for a given binary
predictor, obtained by g-computation.
}
\usage{
boot.glm(formula, family = stats::gaussian(), data, subset,
         weights, na.action, offset,
         conf.level = 0.95, R = 10000,
         irls.maxit = 25L, irls.tol = 1e-8,
         effect = c("coef", "rr", "rd"), exposure = NULL,
         seed = 123, \dots)

\method{summary}{boot.glm}(object, \dots)

\method{print}{summary.boot.glm}(x, digits = max(3L, getOption("digits") - 3L), \dots)
}
\arguments{
  \item{formula}{an object of class \code{"formula"}, as in \code{glm}.}
  \item{family}{as in \code{\link[stats]{glm}} (e.g. \code{binomial()},
    \code{poisson()}, \code{gaussian()}).}
  \item{data}{a data frame containing the variables in the model.}
  \item{subset, weights, na.action, offset}{as in \code{\link[stats]{glm}}.}
  \item{conf.level}{confidence level for the bootstrap intervals.}
  \item{R}{number of bootstrap replicates (default \code{10000}).}
  \item{irls.maxit}{maximum number of IRLS (Newton) iterations per
    bootstrap replicate (default \code{25}).}
  \item{irls.tol}{relative convergence tolerance on the coefficient
    update between IRLS iterations (default \code{1e-8}).}
  \item{effect}{which quantity to bootstrap and report: \code{"coef"}
    (default) for the ordinary model coefficients; \code{"rr"} for the
    marginal (population-averaged) risk ratio associated with
    \code{exposure}, obtained by g-computation; or \code{"rd"} for the
    corresponding marginal risk difference. \code{"rr"}/\code{"rd"} are
    only available when \code{family} is binomial.}
  \item{exposure}{the name (character string) of a single binary
    predictor in \code{formula} for which the marginal risk ratio/
    difference is computed; required (and only used) when \code{effect}
    is \code{"rr"} or \code{"rd"}. Must correspond to a two-level factor
    or a 0/1-coded main-effect term (no interactions); see Details.}
  \item{seed}{integer seed used to make the bootstrap resampling
    reproducible; set via \code{\link[base]{set.seed}} at the start of
    the function. Defaults to \code{123}; pass \code{NULL} to use
    whatever random state is currently active (no reseeding), or any
    other integer for a different reproducible draw.}
  \item{object}{a \code{boot.glm} object, for \code{summary}.}
  \item{x}{a \code{summary.boot.glm} object, for \code{print}.}
  \item{digits}{number of significant digits for printing.}
  \item{\dots}{further arguments passed to \code{\link[stats]{glm}} for
    the single full-data fit (e.g. \code{contrasts}); does not affect
    the bootstrap replicates.}
}
\details{
\strong{\code{effect = "coef"}} (default): each bootstrap replicate is
refit with a lean, hand-rolled IRLS solver, warm-started from the
full-data MLE coefficients, using a Cholesky factorization and two
backsolves per iteration, and converges typically in 2-4 iterations
thanks to the warm start. The special case
\code{family = gaussian(link = "identity")} skips iteration entirely and
solves the (weighted) normal equations exactly in one step. The
acceleration constant for the BCa interval is approximated from a
one-step Newton (infinitesimal jackknife) update using the converged
IRLS working weights, avoiding n full leave-one-out refits.

\strong{\code{effect = "rr"} / \code{"rd"}} (binomial only):
g-computation ("standardization") estimates the marginal effect of
\code{exposure} by predicting, for every row of the data, the response
probability the model would give if that row's exposure were set to 1,
and again as if it were set to 0 -- leaving every other covariate at its
observed value -- and averaging each set of predictions over the whole
sample: \eqn{R_1 = \mathrm{mean}(\hat p_i \mid \mathrm{exposure}_i = 1)},
\eqn{R_0 = \mathrm{mean}(\hat p_i \mid \mathrm{exposure}_i = 0)}. The
marginal risk ratio is \eqn{R_1 / R_0} and the marginal risk difference
is \eqn{R_1 - R_0}. Each bootstrap replicate resamples whole rows
(outcome and covariates together), refits the same lean IRLS solver used
for \code{"coef"}, and then repeats the standardization step within that
replicate's own resampled covariate distribution -- the standard
nonparametric bootstrap for a g-computed effect. The BCa acceleration
constant uses a closed-form (no-refit) leave-one-out approximation built
from the same one-step Newton update used for \code{"coef"}. The p-value
is obtained by CI inversion against the natural null value for each
scale (1 for the risk ratio, 0 for the risk difference).
\code{exposure} must resolve to exactly one design-matrix column (a
plain 0/1 numeric/logical variable, or a two-level factor with the usual
treatment contrasts) -- multi-level factors, interactions involving
\code{exposure}, or non-0/1-coded columns are not supported and raise an
error.
}
\value{
An object of class \code{c("boot.glm", "glm", "lm")}: identical to what
\code{\link[stats]{glm}} returns, with an additional \code{boot} element
(the bootstrap replicates of whichever quantity \code{effect} selected,
its BCa CI, CI-inversion p-value, \code{R}, and \code{effect} itself).
All standard \code{glm} methods keep working unchanged, regardless of
\code{effect}; use \code{summary()} for the bootstrap inference table
(columns \code{Estimate}, \code{CI lower}, \code{CI upper}, and
\code{Pr(>|z|)}).
}
\examples{
fit <- boot.glm(am ~ wt + hp, data = mtcars, family = binomial(), R = 500)
summary(fit)

## marginal risk ratio / risk difference for a binary predictor,
## via g-computation
d <- mtcars; d$vs <- factor(d$vs)
rr_fit <- boot.glm(am ~ vs + wt, data = d, family = binomial(),
                    R = 500, effect = "rr", exposure = "vs")
summary(rr_fit)
}
\seealso{
\code{\link[stats]{glm}}
}
