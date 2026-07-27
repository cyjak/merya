\name{boot.lm}
\alias{boot.lm}
\alias{summary.boot.lm}
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
eta-squared.
}
\usage{
boot.lm(formula, data, subset, weights, na.action, offset,
        conf.level = 0.95, R = 10000,
        effect = c("coef", "partial.cor", "partial.eta2", "eta2"),
        seed = 123, \dots)

\method{summary}{boot.lm}(object, \dots)

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
  \item{seed}{integer seed used to make the bootstrap resampling
    reproducible; set via \code{\link[base]{set.seed}} at the start of
    the function. Defaults to \code{123}; pass \code{NULL} to use
    whatever random state is currently active (no reseeding), or any
    other integer for a different reproducible draw.}
  \item{object}{a \code{boot.lm} object, for \code{summary}.}
  \item{x}{a \code{summary.boot.lm} object, for \code{print}.}
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
refit. Both come directly from quantities already produced by the same
per-replicate Cholesky solve used for \code{"coef"} (no separate reduced-
model refits are needed), and their jackknife (for the BCa acceleration
constant) is likewise closed-form, extending the same leave-one-out
identities used for \code{"coef"} (leave-one-out RSS via the PRESS
identity, leave-one-out \eqn{(X'X)^{-1}} via a rank-one/Sherman-Morrison
update, and a closed-form leave-one-out total sum of squares) -- so,
exactly as for \code{"coef"}, no observation or replicate is ever refit
for the jackknife step.

Because partial eta-squared and eta-squared are \emph{squares} of a
statistic that can be negative, their confidence intervals are built in
two steps: first the BCa interval is computed on the signed (partial/
semi-partial correlation) scale, exactly as for \code{"partial.cor"};
then that signed interval is converted to the squared scale by taking
the smallest absolute value attained within it as the lower bound (0
whenever the signed interval straddles zero) and the largest absolute
value as the upper bound, and squaring both -- rather than naively
squaring the signed endpoints, which would misrepresent the interval
whenever it straddles zero. The p-value (by CI inversion) is computed
once, on the signed scale, since squaring is monotonic in absolute value
and does not change which side of zero a replicate falls on.
}
\value{
An object of class \code{c("boot.lm", "lm")}: identical to what
\code{\link[stats]{lm}} returns, with an additional \code{boot} element
(the bootstrap replicates of whichever quantity \code{effect} selected,
its BCa CIs, CI-inversion p-values, \code{R}, and \code{effect} itself).
All standard \code{lm} methods (\code{predict}, \code{fitted},
\code{residuals}, \code{plot}, ...) keep working unchanged, regardless of
\code{effect}; use \code{summary()} to see the bootstrap inference table
(columns \code{Estimate}, \code{CI lower}, \code{CI upper}, and
\code{Pr(>|z|)}).
}
\examples{
fit <- boot.lm(mpg ~ wt + hp, data = mtcars, R = 500)
summary(fit)

## partial correlations / effect sizes instead of slopes
summary(boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, effect = "partial.cor"))
summary(boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, effect = "eta2"))
}
\seealso{
\code{\link[stats]{lm}}
}
