#' Student's t-Test with BCa Bootstrap Confidence Interval
#'
#' Performs one- and two-sample t-tests, with the same calling convention
#' and output structure as \code{\link[stats]{t.test}}, but the confidence
#' interval is obtained from a bias-corrected and accelerated (BCa)
#' bootstrap of the relevant mean / mean-difference, and the p-value is
#' obtained by inverting that same BCa confidence interval at the
#' hypothesised value \code{mu} (confidence-interval-inversion testing),
#' rather than from the classical t-distribution.
#'
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param alternative a character string specifying the alternative
#'   hypothesis, must be one of \code{"two.sided"} (default),
#'   \code{"greater"} or \code{"less"}.
#' @param mu a number indicating the true value of the mean (or difference
#'   in means if you are performing a two sample test).
#' @param paired a logical indicating whether you want a paired test.
#' @param var.equal currently ignored (present for interface compatibility
#'   with \code{\link[stats]{t.test}}); the bootstrap procedure does not
#'   assume equal variances.
#' @param conf.level confidence level of the interval.
#' @param R number of bootstrap replicates. Default \code{10000}.
#' @param seed integer seed used to make the bootstrap resampling
#'   reproducible; set via \code{\link[base]{set.seed}} at the start of
#'   the function. Defaults to \code{123}; pass \code{NULL} to use
#'   whatever random state is currently active (no reseeding), or any
#'   other integer for a different reproducible draw.
#' @param ... further arguments (ignored; present for compatibility).
#'
#' @return A list with class \code{"htest"} containing the same components
#'   as \code{\link[stats]{t.test}} (\code{statistic}, \code{parameter},
#'   \code{p.value}, \code{conf.int}, \code{estimate}, \code{null.value},
#'   \code{stderr}, \code{alternative}, \code{method}, \code{data.name}),
#'   plus \code{R} (number of bootstrap replicates used).
#'
#' @details The bootstrap statistic is the sample mean (one-sample /
#'   paired case) or the difference in sample means (two-sample case).
#'   Resampling is fully vectorized (a single \code{n x R} index matrix is
#'   drawn and all replicate means are computed with \code{colMeans}), and
#'   the acceleration constant of the BCa interval is obtained from a
#'   closed-form jackknife on the sample mean, so no per-replicate loop is
#'   required. The p-value is obtained by analytically inverting the BCa
#'   transformation at \code{mu} (see package internals), which is O(1)
#'   given the bootstrap distribution -- no search / bisection is used.
#'
#' @examples
#' x <- rnorm(30, mean = 1)
#' boot.t.test(x)
#'
#' y <- rnorm(30, mean = 0)
#' boot.t.test(x, y)
#'
#' boot.t.test(x, y, paired = TRUE)
#'
#' @seealso \code{\link[stats]{t.test}}
#' @export
boot.t.test <- function(x, y = NULL,
                         alternative = c("two.sided", "less", "greater"),
                         mu = 0, paired = FALSE, var.equal = FALSE,
                         conf.level = 0.95, R = 10000, seed = 123, ...) {
  if (!is.null(seed)) set.seed(seed)
  alternative <- match.arg(alternative)
  dname <- deparse(substitute(x))

  if (!missing(y)) dname <- paste(dname, "and", deparse(substitute(y)))

  if (!is.numeric(x)) stop("'x' must be numeric")
  x <- x[!is.na(x)]
  nx <- length(x)
  if (nx < 2) stop("not enough 'x' observations")

  if (!is.null(y)) {
    if (!is.numeric(y)) stop("'y' must be numeric")
    if (paired) {
      if (length(x) != length(y)) stop("'x' and 'y' must have the same length")
      ok <- stats::complete.cases(x, y)
      x <- x[ok]; y <- y[ok]
      nx <- length(x)
      if (nx < 2) stop("not enough 'x' observations")
    } else {
      y <- y[!is.na(y)]
      if (length(y) < 2) stop("not enough 'y' observations")
    }
  }

  if (!missing(conf.level) &&
      (length(conf.level) != 1 || !is.finite(conf.level) ||
       conf.level < 0 || conf.level > 1)) {
    stop("'conf.level' must be a single number between 0 and 1")
  }

  if (is.null(y) || paired) {
    if (paired) {
      z <- x - y
      method <- "Paired bootstrap t-test (BCa CI, CI-inversion p-value)"
      estimate_name <- "mean of the differences"
    } else {
      z <- x
      method <- "One Sample bootstrap t-test (BCa CI, CI-inversion p-value)"
      estimate_name <- "mean of x"
    }
    n <- length(z)
    theta_hat <- mean(z)

    theta_boot <- .boot_mean(z, R)
    theta_loo  <- .jack_mean(z)

    ci <- .bca_ci(theta_boot, theta_hat, theta_loo, conf.level)
    a  <- attr(ci, "a")
    p.value <- .bca_pvalue(mu, theta_boot, theta_hat, a, alternative)

    se <- stats::sd(theta_boot)
    tstat <- (theta_hat - mu) / (stats::sd(z) / sqrt(n))
    df <- n - 1

    if (alternative == "less") {
      ci <- c(-Inf, ci[2])
    } else if (alternative == "greater") {
      ci <- c(ci[1], Inf)
    }

    estimate <- theta_hat
    names(estimate) <- estimate_name

  } else {
    ny <- length(y)
    method <- "Two independent samples bootstrap t-test (BCa CI, CI-inversion p-value)"

    theta_hat <- mean(x) - mean(y)
    theta_boot <- .boot_mean_diff(x, y, R)

    ## acceleration via a combined jackknife: drop one observation at a time
    ## from the pooled sample, recomputing the appropriate group mean only
    jx <- .jack_mean(x); jy <- .jack_mean(y)
    theta_loo <- c(jx - mean(y), mean(x) - jy)

    ci <- .bca_ci(theta_boot, theta_hat, theta_loo, conf.level)
    a  <- attr(ci, "a")
    p.value <- .bca_pvalue(mu, theta_boot, theta_hat, a, alternative)

    se <- stats::sd(theta_boot)
    se_classic <- sqrt(stats::var(x) / length(x) + stats::var(y) / ny)
    tstat <- (theta_hat - mu) / se_classic
    df <- se_classic^4 / ((stats::var(x) / length(x))^2 / (length(x) - 1) +
                            (stats::var(y) / ny)^2 / (ny - 1))

    if (alternative == "less") {
      ci <- c(-Inf, ci[2])
    } else if (alternative == "greater") {
      ci <- c(ci[1], Inf)
    }

    estimate <- c(mean(x), mean(y))
    names(estimate) <- c("mean of x", "mean of y")
  }

  names(tstat) <- "t"
  names(df) <- "df"
  names(mu) <- if (is.null(y) || paired) "mean" else "difference in means"
  ci <- as.numeric(ci)
  attr(ci, "conf.level") <- conf.level

  out <- list(
    statistic   = tstat,
    parameter   = df,
    p.value     = p.value,
    conf.int    = ci,
    estimate    = estimate,
    null.value  = mu,
    stderr      = se,
    alternative = alternative,
    method      = method,
    data.name   = dname,
    R           = R
  )
  class(out) <- "htest"
  out
}
