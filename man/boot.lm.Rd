\name{boot.lm}
\alias{boot.lm}
\alias{summary.boot.lm}
\alias{print.boot.lm}
\alias{print.summary.boot.lm}
\title{Linear Regression with BCa Bootstrap Confidence Intervals}
\description{
Fits a linear model exactly as \code{\link[stats]{lm}} does (same
formula interface, same fitted \code{"lm"} object), but additionally
performs a bootstrap and attaches bias-corrected and accelerated (BCa)
confidence intervals and CI-inversion p-values, retrievable via
\code{summary()}. By default the bootstrapped quantity is the slope
coefficients (as before); set \code{effect} to instead report, for
each predictor, its partial correlation, partial eta-squared, or
eta-squared.
}
\usage{
boot.lm(formula, data, subset, weights, na.action, offset,
        conf.level = 0.95, R = 10000,
        boot.method = c("wild", "case"),
        wild.dist = c("rademacher", "mammen", "normal"),
        ci.type = c("bca", "percentile"),
        effect = c("coef", "partial.cor", "partial.eta2", "eta2"),
        pred.r.squared = FALSE, seed = 123, \dots)

\method{summary}{boot.lm}(object, \dots)

\method{print}{boot.lm}(x, \dots)

\method{print}{summary.boot.lm}(x, digits = max(3L, getOption("digits") - 3L), \dots)
}
\arguments{
  \item{formula}{an object of class \code{"formula"}, as in \code{lm}.}
  \item{data}{a data frame containing the variables in the model.}
  \item{subset, weights, na.action, offset}{as in \code{\link[stats]{lm}}.}
  \item{conf.level}{confidence level for the bootstrap intervals.}
  \item{R}{number of bootstrap replicates (default \code{10000}).}
  \item{boot.method}{resampling method: \code{"wild"} (default) or
    \code{"case"}. \code{"wild"} keeps the design matrix fixed and
    perturbs each observation's own residual by an independent random
    multiplier (see \code{wild.dist}); because the design matrix never
    changes across replicates, it needs no explicit loop over replicates
    at all, making it both the default and the faster option in most
    cases. \code{"case"} is the classical resample-the-rows-with-
    replacement bootstrap used by every earlier version of this
    function; it remains available because it also resamples the
    covariate distribution itself, which the wild bootstrap does not.}
  \item{wild.dist}{distribution of the wild-bootstrap multipliers
    \eqn{v_i} (only used when \code{boot.method = "wild"}):
    \code{"rademacher"} (default; \eqn{v = \pm 1} with equal
    probability), \code{"mammen"} (a skewed two-point distribution
    matching the third moment), or \code{"normal"} (\eqn{v \sim N(0,1)}).}
  \item{ci.type}{how the confidence interval (and, by CI inversion, the
    p-value) is obtained from the bootstrap distribution: \code{"bca"}
    (default) for the bias-corrected and accelerated interval, or
    \code{"percentile"} for the plain empirical
    \eqn{[\alpha/2, 1-\alpha/2]} quantiles with no bias-correction or
    acceleration adjustment. Applies uniformly regardless of
    \code{boot.method} or \code{effect}.}
  \item{effect}{which quantity to bootstrap and report, per predictor:
    \code{"coef"} (default) for the ordinary slope coefficients;
    \code{"partial.cor"} for the partial correlation between the
    response and each predictor (controlling for the others);
    \code{"partial.eta2"} for partial eta-squared (the square of the
    partial correlation); or \code{"eta2"} for eta-squared (the square
    of the semi-partial/"part" correlation, i.e. the predictor's unique
    share of the total variance in the response). The intercept is
    never reported for the three effect-size options (only for
    \code{"coef"}).}
  \item{pred.r.squared}{logical; if \code{TRUE}, additionally compute the
    textbook "predicted R-squared" (\code{1 - PRESS/TSS}, leave-one-out
    cross-validated), reported alongside whichever \code{effect} was
    requested. Unlike ordinary or partial/semi-partial eta-squared, this
    is not bounded below by 0 -- a model that predicts worse than just
    the response mean legitimately has a negative predicted R-squared,
    and this is preserved rather than treated as bootstrap noise around
    a non-negative quantity. \code{FALSE} by default (not computed). See
    Details.}
  \item{seed}{integer seed used to make the bootstrap resampling
    reproducible; set via \code{\link[base]{set.seed}} at the start of
    the function. Defaults to \code{123}; pass \code{NULL} to use
    whatever random state is currently active (no reseeding), or any
    other integer for a different reproducible draw.}
  \item{object}{a \code{boot.lm} object, for \code{summary}.}
  \item{x}{a \code{boot.lm} object (for \code{print}, which delegates to
    \code{summary()} so the console default always reflects \code{effect}),
    or a \code{summary.boot.lm} object (for \code{print.summary.boot.lm}).}
  \item{digits}{number of significant digits for printing.}
  \item{\dots}{further arguments passed to \code{\link[stats]{lm}} (e.g.
    \code{contrasts}).}
}
\details{
\strong{\code{boot.method = "wild"}} (default): keeps the design matrix
\eqn{X} (and hence \eqn{X'X}, its one-time Cholesky factorization, and
the hat-matrix diagonal) fixed across every replicate, and only
perturbs the fitted model's own residuals,
\eqn{y^*_i = \hat y_i + v_i(y_i - \hat y_i)}, where \eqn{v_i} are
i.i.d. mean-0, variance-1 multipliers drawn from \code{wild.dist}.
Because \eqn{X} never changes, the entire \code{R}-replicate bootstrap
for every \code{effect} option (and for \code{pred.r.squared}) reduces
to a handful of full-matrix operations -- one \code{crossprod()}, one
\code{backsolve()} handling all \code{R} replicates' right-hand sides
at once, and elementwise arithmetic on \code{n x R}/length-\code{R}
objects -- with no explicit loop over replicates at all. The
leave-one-out predictions needed for \code{pred.r.squared} likewise
reuse the fixed hat values via the closed-form identity
\eqn{\hat y_{(-i)} = \hat y_i - h_i e_i/(1-h_i)}.

\strong{\code{boot.method = "case"}}: resamples rows (both the
response and every predictor) with replacement. Coefficients are
recomputed on each resample via a direct Cholesky solve of the normal
equations, one per replicate, since the design matrix itself changes
every time.

For \strong{both} \code{boot.method}s, the BCa acceleration constant
is estimated from the same closed-form leave-one-out jackknife on the
original, unperturbed sample -- this jackknife's role is to estimate
the curvature/skewness of the estimator's sampling distribution from
each observation's influence, which does not depend on which
resampling scheme generates the main bootstrap distribution, so it is
shared rather than duplicated.

\strong{\code{effect = "partial.cor"} / \code{"partial.eta2"} /
\code{"eta2"}}: for predictor \eqn{j}, the partial correlation is
\code{sign(t_j) * sqrt(t_j^2 / (t_j^2 + df))} where \code{t_j} is that
predictor's t-statistic and \code{df} the residual degrees of freedom
(both from the replicate's own refit); the semi-partial ("part")
correlation is instead \code{t_j * sqrt((1 - R^2) / df)} of that same
refit. This is equivalent to \code{sign(t_j) * sqrt(R^2_full -
R^2_reduced)} -- the (signed) square root of the unique drop in
R-squared from omitting predictor \eqn{j} -- via the general-linear-model
identity \eqn{F = t_j^2} for dropping a single predictor; it is
\emph{not} the partial correlation itself further rescaled by
\code{sqrt(1 - R^2)}, which only coincides with the correct value in the
trivial \code{t_j -> 0} limit and otherwise understates the effect.
Both come directly from quantities already produced by the same
per-replicate refit used for \code{"coef"} (no separate reduced-model
refits are needed), and their jackknife (for the BCa acceleration
constant) is likewise closed-form.

Because partial eta-squared and eta-squared are squares of a
statistic that can be negative, their confidence intervals are built by
first computing the interval on the signed (partial/semi-partial
correlation) scale, then converting to the squared scale by taking the
smallest absolute value attained within that interval as the lower
bound (0 whenever the signed interval straddles zero) and the largest
absolute value as the upper bound, and squaring both. The p-value (by
CI inversion) is computed once, on the signed scale, since squaring is
monotonic in absolute value and does not change which side of zero a
replicate falls on.

\strong{\code{pred.r.squared = TRUE}}: for each observation \eqn{j}, its
leave-one-out predicted value is \eqn{x_j'\hat\beta_{(-j)}}, where
\eqn{\hat\beta_{(-j)}} is the coefficient vector from the same model
refit excluding observation \eqn{j} -- obtained in closed form via the
hat-matrix identity (no refitting). From these, \eqn{\mathrm{PRESS} =
\sum_j (y_j - x_j'\hat\beta_{(-j)})^2} and predicted R-squared is
\eqn{1 - \mathrm{PRESS}/\mathrm{TSS}} -- the standard definition of
predicted R-squared. This quantity is already on its natural, signed
scale (it is not a square of anything), so it is bootstrapped directly
and given a confidence interval and CI-inversion p-value with no
separate sign-handling transform of any kind: predicted R-squared's
negative values are directly meaningful (the model predicts worse than
the mean) and are not an artifact to be corrected for.
}
\value{
An object of class \code{c("boot.lm", "lm")}: identical to what
\code{\link[stats]{lm}} returns, with an additional \code{boot} element
(the bootstrap replicates of whichever quantity \code{effect} selected,
its confidence intervals, CI-inversion p-values, \code{R},
\code{boot.method}, \code{ci.type}, and \code{effect} itself). All
standard \code{lm} methods (\code{predict}, \code{fitted},
\code{residuals}, \code{plot}, ...) keep working unchanged, regardless
of \code{effect} or \code{boot.method}; use \code{summary()} to see the
bootstrap inference table (columns \code{Estimate}, \code{CI lower},
\code{CI upper}, and \code{Pr(>|z|)}).
}
\examples{
fit <- boot.lm(mpg ~ wt + hp, data = mtcars, R = 500)
summary(fit)

## the classical case-resampling bootstrap, instead of the wild-
## bootstrap default
summary(boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, boot.method = "case"))

## partial correlations / effect sizes instead of slopes
summary(boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, effect = "partial.cor"))
summary(boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, effect = "eta2"))

## also report predicted (leave-one-out) R-squared
summary(boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, pred.r.squared = TRUE))
}
\seealso{
\code{\link[stats]{lm}}
}
