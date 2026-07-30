\name{boot.glm}
\alias{boot.glm}
\alias{summary.boot.glm}
\alias{print.boot.glm}
\alias{print.summary.boot.glm}
\title{Generalized Linear Regression with BCa Bootstrap Confidence Intervals}
\description{
Fits a generalized linear model exactly as \code{\link[stats]{glm}} does
(same formula interface, same fitted \code{"glm"} object), but
additionally performs a case-resampling bootstrap and attaches
bias-corrected and accelerated (BCa) confidence intervals and
CI-inversion p-values, retrievable via \code{summary()}. By default the
bootstrapped quantity is the model coefficients (as before); \code{effect}
can instead request odds ratios, risk ratios/differences, or the
population attributable fraction.
}
\usage{
boot.glm(formula, family = stats::gaussian(), data, subset,
         weights, na.action, offset,
         conf.level = 0.95, R = 10000,
         irls.maxit = 25L, irls.tol = 1e-8,
         effect = c("coef", "OR", "RR", "RD", "PAF"),
         exposure = NULL, seed = 123, \dots)

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
  \item{irls.maxit}{maximum number of IRLS (Newton) iterations per
    bootstrap replicate (default \code{25}).}
  \item{irls.tol}{relative convergence tolerance on the coefficient
    update between IRLS iterations (default \code{1e-8}).}
  \item{effect}{which quantity to bootstrap and report:
    \describe{
      \item{\code{"coef"}}{(default) the ordinary model coefficients.}
      \item{\code{"OR"}}{odds ratios: the exact same bootstrap/jackknife
        coefficient replicates as \code{"coef"}, simply exponentiated
        (estimate and CI); the p-value is untouched. Requires
        \code{family = binomial(link = "logit")}.}
      \item{\code{"RR"}}{risk ratio. If \code{exposure} is \code{NULL},
        reports the conditional risk ratio for every predictor as the
        exponentiated coefficient (same replicates as \code{"coef"},
        just exponentiated, p-value untouched); requires a log link
        (any family). If \code{exposure} is given, instead reports the
        marginal risk ratio for that one binary predictor via
        g-computation; requires \code{family = binomial()} (any link).}
      \item{\code{"RD"}}{the marginal risk difference for
        \code{exposure} via g-computation; requires \code{family =
        binomial()} and a non-\code{NULL} \code{exposure}.}
      \item{\code{"PAF"}}{the population attributable fraction for every
        eligible predictor term via g-computation; requires \code{family
        = binomial()}. \code{exposure} is not used.}
    }
  }
  \item{exposure}{the name (character string) of a single binary
    predictor in \code{formula}; required for \code{effect = "RD"},
    optional for \code{effect = "RR"} (see above), and unused otherwise.
    Must correspond to a two-level factor or a 0/1-coded main-effect
    term (no interactions); see Details.}
  \item{seed}{integer seed used to make the bootstrap resampling
    reproducible; set via \code{\link[base]{set.seed}} at the start of
    the function. Defaults to \code{123}; pass \code{NULL} to use
    whatever random state is currently active (no reseeding), or any
    other integer for a different reproducible draw.}
  \item{object}{a \code{boot.glm} object, for \code{summary}.}
  \item{x}{a \code{boot.glm} object (for \code{print}, which delegates to
    \code{summary()} so the console default always reflects \code{effect}),
    or a \code{summary.boot.glm} object (for \code{print.summary.boot.glm}).}
  \item{digits}{number of significant digits for printing.}
  \item{\dots}{further arguments passed to \code{\link[stats]{glm}} for
    the single full-data fit (e.g. \code{contrasts}); does not affect
    the bootstrap replicates.}
}
\details{
\strong{\code{effect = "coef"}} (default): each bootstrap replicate is
refit with a lean, hand-rolled IRLS solver, warm-started from the
full-data MLE coefficients, converging in a handful of iterations. The
acceleration constant for the BCa interval is approximated from a
one-step Newton (infinitesimal jackknife) update using the converged
IRLS working weights, avoiding n full leave-one-out refits.

\strong{\code{effect = "OR"}} and the conditional-risk-ratio form of
\strong{\code{effect = "RR"}} (\code{exposure = NULL}): these reuse the
identical bootstrap/jackknife coefficient replicates as \code{"coef"} --
no separate resampling is done -- and simply report
\code{exp(estimate)}/\code{exp(CI)} instead of the raw coefficient and
its CI. Because the BCa interval's endpoints are specific
(bias-and-acceleration-corrected) quantiles of the bootstrap
distribution, and \code{exp()} is monotonic, exponentiating those
endpoints gives the exact corresponding BCa interval on the OR/RR
scale. The p-value is left as computed on the coefficient scale (testing
coefficient = 0, equivalently OR/RR = 1), since exponentiation does not
change which side of the null a replicate falls on.

\strong{The g-computation form of \code{effect = "RR"}} (\code{exposure}
given) \strong{, \code{effect = "RD"}, and \code{effect = "PAF"}}: all
three use g-computation ("standardization"). For \code{"RR"}/\code{"RD"},
every row's \code{exposure} is set to 1 for everyone, then to 0 for
everyone (leaving every other covariate at its observed value), giving
population-averaged risks \eqn{R_1} and \eqn{R_0}; the marginal risk
ratio is \eqn{R_1 / R_0} and the marginal risk difference is
\eqn{R_1 - R_0}. For \code{"PAF"}, the true (observed) prevalence
\eqn{P} is compared, for every eligible predictor term \eqn{i} in turn,
against the counterfactual prevalence \eqn{P_i} obtained by setting that
one term to its reference category for everyone (all of its
design-matrix dummy columns to 0) if it is categorical, or to the sample
mean for everyone if it is continuous, again leaving all other
covariates as observed; the population attributable fraction is
\eqn{\mathrm{PAF}_i = (P - P_i) / P}. Multi-level categorical predictors
contribute a single row (one counterfactual scenario: everyone at the
reference level), and only simple, untransformed, non-interaction
predictor terms are supported for \code{"PAF"} (a term must match a raw
column name in \code{data}; interactions or transformations such as
\code{poly()}/\code{log()} raise an error).

For all three g-computation cases, each bootstrap replicate resamples
whole rows (outcome and covariates together), refits the same lean IRLS
solver used for \code{"coef"}, and repeats the standardization step
within that replicate's own resampled covariate distribution and refit
coefficients. The BCa acceleration constant uses a closed-form (no-refit)
leave-one-out approximation built from the same one-step Newton update
used for \code{"coef"} (for \code{"PAF"}'s continuous-predictor recipes,
the reference value used in this jackknife step is fixed at the
full-data mean rather than recomputed per leave-one-out pseudo-replicate,
a small approximation that keeps the computation vectorized). The
p-value is obtained by CI inversion against the natural null value for
each scale (1 for the risk ratio, 0 for the risk difference or the
attributable fraction). \code{exposure} must resolve to exactly one
design-matrix column (a plain 0/1 numeric/logical variable, or a
two-level factor with the usual treatment contrasts) -- multi-level
factors, interactions involving \code{exposure}, or non-0/1-coded
columns are not supported and raise an error.
}
\value{
An object of class \code{c("boot.glm", "glm", "lm")}: identical to what
\code{\link[stats]{glm}} returns, with an additional \code{boot} element
(the bootstrap replicates of whichever quantity \code{effect} selected,
its BCa CI(s), CI-inversion p-value(s), \code{R}, and \code{effect}
itself). All standard \code{glm} methods keep working unchanged,
regardless of \code{effect}; use \code{summary()} for the bootstrap
inference table (columns \code{Estimate}, \code{CI lower}, \code{CI
upper}, and \code{Pr(>|z|)}; no separate bootstrap standard error
column).
}
\examples{
fit <- boot.glm(am ~ wt + hp, data = mtcars, family = binomial(), R = 500)
summary(fit)

## odds ratios (requires a logit link)
summary(boot.glm(am ~ wt + hp, data = mtcars, family = binomial(),
                  R = 500, effect = "OR"))

## conditional risk ratios for every predictor (requires a log link)
summary(boot.glm(am ~ wt + hp, data = mtcars, family = binomial(link = "log"),
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
