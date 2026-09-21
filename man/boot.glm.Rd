\name{boot.glm}
\alias{boot.glm}
\alias{summary.boot.glm}
\alias{print.boot.glm}
\alias{print.summary.boot.glm}
\title{Generalized Linear Regression with BCa Bootstrap Confidence Intervals}
\description{
Fits a generalized linear model exactly as \code{\link[stats]{glm}} does
(same formula interface, same fitted \code{"glm"} object), but
additionally performs a bootstrap and attaches bias-corrected and
accelerated (BCa) confidence intervals and CI-inversion p-values,
retrievable via \code{summary()}. By default the bootstrapped quantity
is the model coefficients (as before); \code{effect} can instead
request odds ratios, risk ratios/differences, or the population
attributable fraction -- see Details.
}
\usage{
boot.glm(formula, family = stats::gaussian(), data, subset,
         weights, na.action, offset,
         conf.level = 0.95, R = 10000,
         boot.method = c("wild", "case"),
         wild.dist = c("rademacher", "mammen", "normal"),
         ci.type = c("bca", "percentile"),
         seed = 123,
         irls.maxit = 25L, irls.tol = 1e-8,
         effect = c("coef", "OR", "RR", "RD", "PAF"),
         exposure = NULL, \dots)

\method{summary}{boot.glm}(object, \dots)

\method{print}{boot.glm}(x, \dots)

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
  \item{boot.method}{resampling method: \code{"wild"} (default) or
    \code{"case"}. \code{"wild"} keeps the design matrix fixed and
    perturbs the residual of the IRLS working response at convergence by
    an independent random multiplier (see \code{wild.dist}), then solves
    one weighted-least-squares step per replicate using the converged
    IRLS weights -- a linearization around the full-data fit rather than
    a fresh nonlinear IRLS refit -- which is well-defined for any
    family/link and needs no explicit loop over replicates at all,
    making it both the default and the faster option in most cases.
    \code{"case"} is the classical resample-the-rows-with-replacement
    bootstrap (a full nonlinear IRLS refit per replicate); it remains
    available whenever an exact nonlinear refit, or resampling the
    covariate distribution itself, is preferred.}
  \item{wild.dist}{distribution of the wild-bootstrap multipliers
    \eqn{v_i} (only used when \code{boot.method = "wild"}); see
    \code{\link{boot.lm}}.}
  \item{ci.type}{how the confidence interval (and, by CI inversion, the
    p-value) is obtained from the bootstrap distribution: \code{"bca"}
    (default) or \code{"percentile"}; see \code{\link{boot.lm}}. Applies
    uniformly regardless of \code{boot.method} or \code{effect}.}
  \item{seed}{integer seed used to make the bootstrap resampling
    reproducible; set via \code{\link[base]{set.seed}} at the start of
    the function. Defaults to \code{123}; pass \code{NULL} to use
    whatever random state is currently active (no reseeding), or any
    other integer for a different reproducible draw.}
  \item{irls.maxit}{maximum number of IRLS (Newton) iterations per
    bootstrap replicate when \code{boot.method = "case"} (default
    \code{25}). Unused when \code{boot.method = "wild"}.}
  \item{irls.tol}{relative convergence tolerance on the coefficient
    update between IRLS iterations when \code{boot.method = "case"}
    (default \code{1e-8}). Unused when \code{boot.method = "wild"}.}
  \item{effect}{which quantity to bootstrap and report: \code{"coef"}
    (default) for the ordinary model coefficients; \code{"OR"} for odds
    ratios (the same bootstrap/jackknife coefficient replicates as
    \code{"coef"}, exponentiated; requires
    \code{family = binomial(link = "logit")}); \code{"RR"} for the
    risk ratio -- if \code{exposure} is \code{NULL}, the *conditional*
    risk ratio for every predictor (exponentiated coefficient, requires
    a log link, any family), or if \code{exposure} is given, the
    *marginal* risk ratio for that predictor via g-computation (requires
    \code{family = binomial()}); \code{"RD"} for the marginal risk
    difference via g-computation (requires \code{family = binomial()}
    and a non-\code{NULL} \code{exposure}); or \code{"PAF"} for the
    population attributable fraction for every eligible predictor term
    via g-computation (requires \code{family = binomial()}).}
  \item{exposure}{the name (character string) of a single binary
    predictor in \code{formula}; required for \code{effect = "RD"},
    optional for \code{effect = "RR"} (see above), and unused otherwise.
    Must correspond to a two-level factor or a 0/1-coded main-effect
    term (no interactions); see Details.}
  \item{object}{a \code{boot.glm} object, for \code{summary}.}
  \item{x}{a \code{boot.glm} object (for \code{print}, which delegates
    to \code{summary()} so the console default always reflects
    \code{effect}), or a \code{summary.boot.glm} object (for
    \code{print.summary.boot.glm}).}
  \item{digits}{number of significant digits for printing.}
  \item{\dots}{further arguments passed to \code{\link[stats]{glm}} for
    the single full-data fit (e.g. \code{contrasts}); does not affect
    the bootstrap replicates.}
}
\details{
\strong{\code{boot.method = "wild"}} (default): perturbs the residual of
the IRLS working response at convergence,
\eqn{z_i = \eta_i + (y_i-\mu_i)/g'(\mu_i)}, by an independent random
multiplier \eqn{v_i}, and solves one weighted-least-squares step using
the fixed, converged IRLS weights. Because the weighted design is fixed
across replicates, the entire \code{R}-replicate bootstrap for the
coefficients reduces to a handful of full-matrix operations with no
explicit loop over replicates, exactly as for \code{\link{boot.lm}}'s
wild bootstrap. For the g-computation effects (\code{"RR"}/\code{"RD"}
with \code{exposure}, and \code{"PAF"}), the standardization step uses
the fixed, full covariate distribution together with each replicate's
own coefficients, again fully vectorized across all \code{R} replicates
at once. For \code{"PAF"} specifically, the observed prevalence used as
the denominator does not vary under wild bootstrap (there is no
resampled response for it to be computed from) and is held fixed at the
full-sample value for every replicate; only the model-based
counterfactual prevalence varies per replicate.

\strong{\code{boot.method = "case"}}: each bootstrap replicate is refit
with a lean, hand-rolled IRLS solver, warm-started from the full-data
MLE coefficients. For the g-computation effects, the standardization
step is repeated within that replicate's own resampled covariate
distribution.

For \strong{both} \code{boot.method}s, the BCa acceleration constant is
approximated from a one-step Newton (infinitesimal jackknife) update on
the original, unperturbed sample, using the converged IRLS working
weights.

\strong{\code{effect = "OR"}} and the conditional-risk-ratio form of
\strong{\code{effect = "RR"}} (\code{exposure = NULL}): reuse the
identical bootstrap/jackknife coefficient replicates as \code{"coef"}
and simply report \code{exp(estimate)}/\code{exp(CI)} instead of the
raw coefficient and its CI (valid because the interval endpoints are
specific quantiles and \code{exp()} is monotonic, so it commutes with
taking quantiles). The p-value is left as computed on the coefficient
scale.

\strong{The g-computation form of \code{effect = "RR"}} (\code{exposure}
given) \strong{, \code{effect = "RD"}, and \code{effect = "PAF"}}: all
three use g-computation ("standardization"). For \code{"RR"}/\code{"RD"},
every row's \code{exposure} is set to 1 for everyone, then to 0 for
everyone, giving population-averaged risks \eqn{R_1} and \eqn{R_0}; the
marginal risk ratio is \eqn{R_1/R_0} and the marginal risk difference is
\eqn{R_1-R_0}. For \code{"PAF"}, the observed prevalence \eqn{P} is
compared, for every eligible predictor term in turn, against the
counterfactual prevalence obtained by setting that term to its
reference category (categorical) or sample mean (continuous) for
everyone; \eqn{\mathrm{PAF}_i = (P-P_i)/P}. Multi-level categorical
predictors contribute a single row, and only simple, untransformed,
non-interaction predictor terms are supported for \code{"PAF"}
(interactions or transformations such as \code{poly()}/\code{log()}
raise an error). \code{"PAF"} works whether or not \code{data} is
supplied. \code{exposure} must resolve to exactly one design-matrix
column.
}
\value{
An object of class \code{c("boot.glm", "glm", "lm")}: identical to what
\code{\link[stats]{glm}} returns, with an additional \code{boot} element
(the bootstrap replicates of whichever quantity \code{effect} selected,
its confidence interval(s), CI-inversion p-value(s), \code{R},
\code{boot.method}, and \code{effect} itself). All standard \code{glm}
methods keep working unchanged, regardless of \code{effect} or
\code{boot.method}; use \code{summary()} for the bootstrap inference
table (columns \code{Estimate}, \code{CI lower}, \code{CI upper}, and
\code{Pr(>|z|)}).
}
\examples{
fit <- boot.glm(am ~ wt + hp, data = mtcars, family = binomial(), R = 500)
summary(fit)

## the classical case-resampling bootstrap, instead of the wild-
## bootstrap default
summary(boot.glm(am ~ wt + hp, data = mtcars, family = binomial(),
                  R = 500, boot.method = "case"))

## odds ratios (requires a logit link)
summary(boot.glm(am ~ wt + hp, data = mtcars, family = binomial(),
                  R = 500, effect = "OR"))

## conditional risk ratios for every predictor (requires a log link,
## but works for any family -- e.g. Poisson, not just binomial)
summary(boot.glm(carb ~ wt, data = mtcars, family = poisson(link = "log"),
                  R = 500, effect = "RR"))

## marginal risk ratio / risk difference for a binary predictor,
## via g-computation
d <- mtcars; d$vs <- factor(d$vs)
rr_fit <- boot.glm(am ~ vs + wt, data = d, family = binomial(),
                    R = 500, effect = "RR", exposure = "vs")
summary(rr_fit)

## population attributable fraction for every eligible predictor
summary(boot.glm(am ~ vs + wt, data = d, family = binomial(),
                  R = 500, effect = "PAF"))
}
\seealso{
\code{\link[stats]{glm}}
}
