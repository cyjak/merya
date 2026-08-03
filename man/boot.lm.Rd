\name{boot.lm}
\alias{boot.lm}
\alias{summary.boot.lm}
\alias{print.boot.lm}
\alias{print.summary.boot.lm}
\title{Linear Regression with BCa Bootstrap Confidence Intervals}
\description{
Fits a linear model exactly as \code{\link[stats]{lm}} does (same
formula interface, same fitted \code{"lm"} object), but additionally
performs a case-resampling bootstrap and attaches bias-corrected and
accelerated (BCa) confidence intervals and CI-inversion p-values,
retrievable via \code{summary()}. By default the bootstrapped quantity is
the slope coefficients (as before); set \code{effect} to instead report,
for each predictor, its partial correlation, partial eta-squared, or
eta-squared; optionally also report a bootstrap predicted (leave-one-out)
R-squared.
}
\usage{
boot.lm(formula, data, subset, weights, na.action, offset,
        conf.level = 0.95, R = 10000,
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
\strong{\code{effect = "coef"}} (default): coefficients are recomputed
on each bootstrap resample via a direct Cholesky solve of the normal
equations, and the BCa acceleration constant comes from the closed-form
OLS leave-one-out formula (hat-matrix diagonal) -- no refitting needed
for the jackknife step.

\strong{\code{effect = "partial.cor"} / \code{"partial.eta2"} /
\code{"eta2"}}: for predictor \eqn{j}, the partial correlation is
\code{sign(t_j) * sqrt(t_j^2 / (t_j^2 + df))} where \code{t_j} is that
predictor's t-statistic and \code{df} the residual degrees of freedom
(both from the replicate's own refit); the semi-partial ("part")
correlation additionally rescales by \code{sqrt(1 - R^2)} of that same
refit. Because partial eta-squared and eta-squared are squares of a
statistic that can be negative, their confidence intervals are built by
first computing the BCa interval on the signed (partial/semi-partial
correlation) scale, then converting to the squared scale by taking the
smallest absolute value attained within that interval as the lower
bound (0 whenever the signed interval straddles zero) and the largest
absolute value as the upper bound, and squaring both.

\strong{\code{pred.r.squared = TRUE}}: for each observation \eqn{j}, its
leave-one-out predicted value is \eqn{x_j'\hat\beta_{(-j)}}, where
\eqn{\hat\beta_{(-j)}} is the coefficient vector from the same model
refit excluding observation \eqn{j} -- obtained in closed form via the
hat-matrix identity (no refitting), the same one used for \code{effect =
"coef"}'s BCa acceleration constant. From these, \eqn{\mathrm{PRESS} =
\sum_j (y_j - x_j'\hat\beta_{(-j)})^2} (the leave-one-out predicted error
sum of squares) and predicted R-squared is \eqn{1 -
\mathrm{PRESS}/\mathrm{TSS}}, where \eqn{\mathrm{TSS}} is the total sum
of squares -- the standard definition of predicted R-squared (as
reported by, e.g., Minitab and most regression-diagnostics references).
This quantity is already on its natural, signed scale (it is not a
square of anything), so it is bootstrapped directly (each replicate
resamples rows, refits, and repeats the same
leave-one-out-within-the-replicate PRESS/TSS calculation on its own
resampled data) and given a BCa confidence interval and CI-inversion
p-value with no separate sign-handling transform of any kind -- unlike
\code{"partial.eta2"}/\code{"eta2"}, which square a signed correlation
that is known to estimate a non-negative population quantity, predicted
R-squared's negative values are directly meaningful (the model predicts
worse than the mean) and are not an artifact to be corrected for. The BCa
acceleration constant is approximated with a fast closed-form
"delete-one" leave-one-out PRESS/TSS (excluding one observation's own
term from the full-sample PRESS sum, and using the same closed-form
leave-one-out TSS identity as \code{"eta2"}, rather than a full nested
double-jackknife), consistent with the closed-form/one-step
approximations used throughout this package.
}
\value{
An object of class \code{c("boot.lm", "lm")}: identical to what
\code{\link[stats]{lm}} returns, with an additional \code{boot} element
(the bootstrap replicates of whichever quantity \code{effect} selected,
its BCa CIs, CI-inversion p-values, \code{R}, \code{effect} itself, and,
if \code{pred.r.squared = TRUE}, a \code{pred.r.squared} sub-list). All
standard \code{lm} methods (\code{predict}, \code{fitted},
\code{residuals}, \code{plot}, ...) keep working unchanged, regardless of
\code{effect}; use \code{summary()} to see the bootstrap inference table
(columns \code{Estimate}, \code{CI lower}, \code{CI upper}, and
\code{Pr(>|z|)}; no separate bootstrap standard error column), plus the
predicted R-squared line when requested.
}
\examples{
fit <- boot.lm(mpg ~ wt + hp, data = mtcars, R = 500)
summary(fit)

## partial correlations / effect sizes instead of slopes
summary(boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, effect = "partial.cor"))
summary(boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, effect = "eta2"))

## also report predicted (leave-one-out) R-squared
summary(boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, pred.r.squared = TRUE))
}
\seealso{
\code{\link[stats]{lm}}
}
