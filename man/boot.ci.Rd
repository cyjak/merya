\name{boot.ci}
\alias{boot.ci}
\alias{boot.ci.default}
\alias{boot.ci.htest}
\alias{boot.ci.lm}
\alias{boot.ci.glm}
\title{BCa Bootstrap Confidence Intervals for an Already-Fitted Object}
\description{
A single generic entry point that takes an object already returned by
\code{\link[stats]{t.test}}, \code{\link[stats]{cor.test}},
\code{\link[stats]{lm}}, or \code{\link[stats]{glm}}, and returns the
same object with its confidence interval(s) replaced by a
bias-corrected and accelerated (BCa) bootstrap interval and its
p-value(s) replaced by the CI-inversion p-value -- i.e. exactly what
the corresponding
\code{\link{boot.t.test}}/\code{\link{boot.cor.test}}/\code{\link{boot.lm}}/
\code{\link{boot.glm}} function would have produced had it been used to
fit the model in the first place.
}
\usage{
boot.ci(object, \dots)

\method{boot.ci}{htest}(object, x, y, R = 10000L, conf.level = NULL,
         seed = 123, \dots)

\method{boot.ci}{lm}(object, R = 10000L, conf.level = 0.95, seed = 123, \dots)

\method{boot.ci}{glm}(object, R = 10000L, conf.level = 0.95, seed = 123, \dots)
}
\arguments{
  \item{object}{a fitted object of class \code{"htest"} (from
    \code{\link[stats]{t.test}} or \code{\link[stats]{cor.test}}),
    \code{"lm"}, or \code{"glm"}.}
  \item{x, y}{(\code{htest} method only) optionally supply the original
    data explicitly, overriding automatic recovery from
    \code{object$data.name}.}
  \item{R}{number of bootstrap replicates.}
  \item{conf.level}{confidence level. For the \code{htest} method,
    \code{NULL} (the default) keeps the original test's confidence level
    when recoverable, else falls back to \code{0.95}.}
  \item{seed}{integer seed used to make the bootstrap resampling
    reproducible; forwarded to the underlying \code{boot.*} call, which
    sets it via \code{\link[base]{set.seed}}. Defaults to \code{123};
    pass \code{NULL} to use whatever random state is currently active
    (no reseeding), or any other integer for a different reproducible
    draw.}
  \item{\dots}{further arguments forwarded to the underlying
    \code{boot.*} function -- e.g. \code{boot.method}, \code{wild.dist},
    \code{ci.type} (all three, see \code{\link{boot.lm}}), and
    \code{effect} for every method; \code{paired} for the \code{htest}
    method; \code{irls.maxit}/\code{irls.tol}/\code{exposure} for
    \code{glm}.}
}
\value{
An object of the \strong{same class, and in the same format}, as
\code{object}: a \code{"htest"} list for \code{t.test}/\code{cor.test}
input, or an \code{"lm"}/\code{"glm"} object for \code{lm}/\code{glm}
input.
}
\details{
For \code{"lm"} and \code{"glm"} objects, \code{boot.ci()} reads the
original fitting call off the object (\code{object$call}), substitutes
in the corresponding \code{boot.*} function, and re-evaluates it in the
environment that called \code{boot.ci()} -- the same mechanism
\code{\link[stats]{update}} uses to refit a model. This means
\code{boot.ci(lm(formula, data = d))} runs the \emph{exact same}
resampling/BCa/CI-inversion code as \code{boot.lm(formula, data = d)}
(down to using the same internal helper functions), so with the same
seed the two are equivalent; it also means the original variables (e.g.
\code{d} above) must still be reachable from wherever \code{boot.ci()}
is called, exactly as they would need to be for \code{update(fit)} to
work.

\code{"htest"} objects from \code{t.test()}/\code{cor.test()} are a
partial exception: base R's \code{t.test}/\code{cor.test} do not store a
\code{$call}, only a deparsed \code{data.name} label (e.g. \code{"x and
y"}). \code{boot.ci()} re-parses and re-evaluates that label in the
calling environment by default, which works whenever \code{x}/\code{y}
were plain variable names (the overwhelmingly common case); pass the
data explicitly via \code{boot.ci(object, x = ..., y = ...)} to
sidestep this whenever \code{x}/\code{y} were one-off expressions (e.g.
\code{t.test(rnorm(10))}).

No new resampling machinery is introduced anywhere in this function:
\code{boot.ci()} is a thin (and therefore fast) dispatcher, so its own
overhead is negligible next to the cost of the bootstrap it triggers.
}
\examples{
## lm: boot.ci(lm(...)) and boot.lm(...) use the same machinery, and
## agree exactly under the shared default seed (123)
fit <- lm(mpg ~ wt + hp, data = mtcars)
b1 <- boot.ci(fit, R = 300)
b2 <- boot.lm(mpg ~ wt + hp, data = mtcars, R = 300)
summary(b1)

## htest
x <- rnorm(30, mean = 1)
boot.ci(t.test(x))
}
\seealso{
\code{\link{boot.t.test}}, \code{\link{boot.cor.test}},
\code{\link{boot.lm}}, \code{\link{boot.glm}}
}
