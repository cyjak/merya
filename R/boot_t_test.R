#' Student's t-Test with BCa Bootstrap Confidence Intervals
#'
#' Performs one-sample, paired, and two independent samples t-tests,
#' using the same arguments as \code{\link[stats]{t.test}} and returning
#' an object of the same \code{"htest"} class, but with a bootstrap
#' confidence interval and CI-inversion p-value in place of the classical
#' asymptotic ones.
#'
#' Internally, every case is expressed as a linear regression and handed
#' to the exact same engines \code{\link{boot.lm}} uses (an intercept-only
#' regression of \code{x} for the one-sample case, of the paired
#' differences \code{x - y} for the paired case, or of the combined
#' sample on a 0/1 group indicator for the two-sample case, testing that
#' indicator's coefficient -- which is algebraically identical to the
#' difference in means). This means \code{boot.t.test(x, y)} and
#' \code{boot.lm(c(x, y) ~ group)}'s group coefficient, for the same
#' \code{boot.method}/\code{wild.dist}/\code{ci.type}/\code{seed}, are
#' not just conceptually comparable but numerically identical -- the
#' same bootstrap replicates, the same jackknife, and the same interval.
#'
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param alternative character string, one of \code{"two.sided"} (default),
#'   \code{"greater"} or \code{"less"}.
#' @param mu a number indicating the true value of the mean (or difference
#'   in means if performing a two-sample test).
#' @param paired logical; if \code{TRUE}, a paired test is performed.
#' @param var.equal currently ignored; kept for interface compatibility
#'   with \code{\link[stats]{t.test}}. Neither bootstrap method assumes
#'   equal variances in the classical (pooled-variance) sense; see
#'   Details.
#' @param conf.level confidence level of the returned interval.
#' @param R number of bootstrap replicates. Default \code{10000}.
#' @param boot.method resampling method: \code{"wild"} (default) or
#'   \code{"case"}; see \code{\link{boot.lm}}, whose engines this
#'   function reuses directly (see Details).
#' @param wild.dist distribution of the wild-bootstrap multipliers
#'   (only used when \code{boot.method = "wild"}); see
#'   \code{\link{boot.lm}}.
#' @param ci.type how the confidence interval (and CI-inversion p-value)
#'   is obtained from the bootstrap distribution: \code{"bca"} (default)
#'   or \code{"percentile"}; see \code{\link{boot.lm}}.
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
#'   plus \code{R}, \code{boot.method}, and \code{ci.type}.
#'
#' @details The bootstrap statistic is the sample mean (one-sample /
#'   paired case) or the difference in sample means (two-sample case),
#'   obtained as the relevant coefficient of an equivalent OLS regression
#'   (see above) via \code{\link{boot.lm}}'s own internal engines --
#'   \code{.boot_lm_coef}/\code{.wild_lm_core} for the bootstrap
#'   distribution, and \code{.jack_lm_coef} for the BCa acceleration
#'   constant's leave-one-out jackknife, all closed-form/vectorized with
#'   no explicit loop over observations or (for \code{"wild"})
#'   replicates. Because this reduces the two-sample case to a
#'   *combined*, unstratified regression on a group indicator, both
#'   \code{boot.method}s here resample/perturb the pooled sample rather
#'   than resampling each group separately at its own fixed size; this is
#'   what guarantees consistency with \code{\link{boot.lm}} for the same
#'   data, and is standard practice for the wild bootstrap in particular
#'   (which has no natural "resample within group" analogue, since it
#'   does not resample rows at all). The p-value is obtained by
#'   analytically inverting the BCa (or, for \code{ci.type =
#'   "percentile"}, empirical-quantile) transformation at \code{mu}, in
#'   the same way as every other function in this package.
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
#' ## the classical case-resampling bootstrap, instead of the wild-
#' ## bootstrap default
#' boot.t.test(x, y, boot.method = "case")
#'
#' ## consistency with boot.lm(): the "group" coefficient below matches
#' ## boot.t.test(x, y)'s estimate/CI/p-value exactly (same seed, method)
#' d <- data.frame(z = c(x, y), group = c(rep(1, length(x)), rep(0, length(y))))
#' boot.lm(z ~ group, data = d)
#'
#' @seealso \code{\link[stats]{t.test}}, \code{\link{boot.lm}}
#' @export
boot.t.test <- function(x, y = NULL,
                         alternative = c("two.sided", "less", "greater"),
                         mu = 0, paired = FALSE, var.equal = FALSE,
                         conf.level = 0.95, R = 10000,
                         boot.method = c("wild", "case"),
                         wild.dist = c("rademacher", "mammen", "normal"),
                         ci.type = c("bca", "percentile"),
                         seed = 123, ...) {
  if (!is.null(seed)) set.seed(seed)
  boot.method <- match.arg(boot.method)
  wild.dist <- match.arg(wild.dist)
  ci.type <- match.arg(ci.type)
  alternative <- match.arg(alternative)
  ci_fn   <- if (ci.type == "bca") .bca_ci else .percentile_ci
  pval_fn <- if (ci.type == "bca") .bca_pvalue else .percentile_pvalue

  dname <- deparse(substitute(x))
  has_y <- !is.null(y)
  if (paired && !has_y) {
    stop("'y' is missing for paired test")
  }
  if (has_y) dname <- paste(dname, "and", deparse(substitute(y)))

  x <- x[is.finite(x)]
  if (length(x) < 2) stop("not enough (finite) 'x' observations")
  if (has_y) {
    y <- y[is.finite(y)]
    if (paired && length(x) != length(y)) {
      stop("'x' and 'y' must have the same length for a paired test")
    }
    if (!paired && length(y) < 2) stop("not enough (finite) 'y' observations")
  }
  if (!missing(mu) && (length(mu) != 1 || !is.finite(mu))) {
    stop("'mu' must be a single finite number")
  }

  ## ---- build the equivalent regression, refit via boot.lm()'s own
  ## engines, and extract the coefficient of interest -- exactly what
  ## boot.lm() itself would do for the same data, guaranteeing identical
  ## results under a shared seed ----
  ci_label <- if (ci.type == "percentile") "percentile" else "BCa"
  if (!has_y || paired) {
    z <- if (paired) x - y else x
    n <- length(z)
    if (paired) {
      method <- sprintf("Paired bootstrap t-test (%s CI, CI-inversion p-value)", ci_label)
      estimate_name <- "mean of the differences"
    } else {
      method <- sprintf("One Sample bootstrap t-test (%s CI, CI-inversion p-value)", ci_label)
      estimate_name <- "mean of x"
    }
    Xd <- matrix(1, nrow = n, ncol = 1L, dimnames = list(NULL, "(Intercept)"))
    beta_hat <- mean(z)
    theta_hat <- beta_hat

    ch <- tryCatch(chol(crossprod(Xd)), error = function(e) NULL)
    XtX_inv <- if (!is.null(ch)) chol2inv(ch) else MASS_ginv_fallback(crossprod(Xd))

    theta_boot <- if (boot.method == "wild") {
      .wild_lm_core(Xd, z, R, XtX_inv, beta_hat, ch, wild.dist,
                    need_pred_loo = FALSE)$boot_coef[, 1L]
    } else {
      .boot_lm_coef(Xd, z, R)[, 1L]
    }
    theta_loo <- .jack_lm_coef(Xd, z, beta_hat, XtX_inv)[, 1L]

    se_classic <- stats::sd(z) / sqrt(n)
    df <- n - 1
    estimate <- theta_hat
    names(estimate) <- estimate_name

  } else {
    nx <- length(x); ny <- length(y)
    n <- nx + ny
    method <- sprintf("Two independent samples bootstrap t-test (%s CI, CI-inversion p-value)",
                       ci_label)

    z <- c(x, y)
    group <- c(rep(1, nx), rep(0, ny))   # group = 1 for x, 0 for y
    Xd <- cbind(`(Intercept)` = 1, group = group)

    ch <- tryCatch(chol(crossprod(Xd)), error = function(e) NULL)
    XtX_inv <- if (!is.null(ch)) chol2inv(ch) else MASS_ginv_fallback(crossprod(Xd))
    Xtz <- crossprod(Xd, z)
    beta_hat <- as.vector(if (!is.null(ch)) {
      backsolve(ch, backsolve(ch, Xtz, transpose = TRUE))
    } else {
      XtX_inv %*% Xtz
    })
    names(beta_hat) <- c("(Intercept)", "group")
    theta_hat <- unname(beta_hat["group"])   # = mean(x) - mean(y)

    theta_boot_mat <- if (boot.method == "wild") {
      .wild_lm_core(Xd, z, R, XtX_inv, beta_hat, ch, wild.dist,
                    need_pred_loo = FALSE)$boot_coef
    } else {
      .boot_lm_coef(Xd, z, R)
    }
    theta_boot <- theta_boot_mat[, 2L]
    theta_loo <- .jack_lm_coef(Xd, z, beta_hat, XtX_inv)[, 2L]

    se_classic <- sqrt(stats::var(x) / nx + stats::var(y) / ny)
    df <- se_classic^4 / ((stats::var(x) / nx)^2 / (nx - 1) +
                           (stats::var(y) / ny)^2 / (ny - 1))
    estimate <- c(mean(x), mean(y))
    names(estimate) <- c("mean of x", "mean of y")
  }

  ok <- stats::complete.cases(theta_boot)
  tb <- theta_boot[ok]
  if (length(tb) < 10) {
    stop("boot.t.test(): fewer than 10 valid bootstrap replicates; ",
         "cannot build a confidence interval.")
  }
  ci_out <- ci_fn(tb, theta_hat, theta_loo, conf.level)
  p.value <- pval_fn(mu, tb, theta_hat, attr(ci_out, "a"), alternative)
  ci <- as.numeric(ci_out)
  if (alternative == "less") {
    ci <- c(-Inf, ci[2L])
  } else if (alternative == "greater") {
    ci <- c(ci[1L], Inf)
  }
  attr(ci, "conf.level") <- conf.level

  tstat <- (theta_hat - mu) / se_classic
  names(tstat) <- "t"
  names(df) <- "df"
  names(mu) <- if (!has_y || paired) "mean" else "difference in means"

  out <- list(
    statistic   = tstat,
    parameter   = df,
    p.value     = p.value,
    conf.int    = ci,
    estimate    = estimate,
    null.value  = mu,
    stderr      = stats::sd(tb),
    alternative = alternative,
    method      = method,
    data.name   = dname,
    R           = R,
    boot.method = boot.method,
    ci.type     = ci.type
  )
  if (boot.method == "wild") out$wild.dist <- wild.dist
  class(out) <- "htest"
  out
}
