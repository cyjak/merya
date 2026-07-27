#' Correlation Test with BCa Bootstrap Confidence Interval
#'
#' Test for association between paired samples, using the same calling
#' convention and output structure as \code{\link[stats]{cor.test}}, but
#' with a bias-corrected and accelerated (BCa) bootstrap confidence
#' interval for the correlation coefficient, and a p-value obtained by
#' inverting that interval at zero correlation (confidence-interval
#' inversion testing).
#'
#' @param x,y numeric vectors of the same length (at least 3 complete
#'   observations required).
#' @param method a character string indicating which correlation
#'   coefficient to use: \code{"pearson"} (default), \code{"kendall"}, or
#'   \code{"spearman"}.
#' @param alternative a character string specifying the alternative
#'   hypothesis: \code{"two.sided"} (default), \code{"greater"} or
#'   \code{"less"}.
#' @param conf.level confidence level for the returned confidence interval.
#' @param R number of bootstrap replicates. Default \code{10000}.
#' @param seed integer seed used to make the bootstrap resampling
#'   reproducible; set via \code{\link[base]{set.seed}} at the start of
#'   the function. Defaults to \code{123}; pass \code{NULL} to use
#'   whatever random state is currently active (no reseeding), or any
#'   other integer for a different reproducible draw.
#' @param ... further arguments (ignored; present for compatibility).
#'
#' @return A list with class \code{"htest"} with the same components as
#'   \code{\link[stats]{cor.test}} (\code{statistic}, \code{parameter},
#'   \code{p.value}, \code{estimate}, \code{null.value}, \code{alternative},
#'   \code{method}, \code{data.name}, \code{conf.int}), plus \code{R}.
#'
#' @details For \code{method = "pearson"} and \code{"spearman"}, the
#'   bootstrap replicates are computed in a fully vectorized way: a single
#'   \code{n x R} resampling-index matrix is drawn, and all R replicate
#'   correlations are obtained from columnwise sufficient statistics
#'   (sums, sums of squares/products) rather than a per-replicate loop,
#'   which is what makes the procedure fast even for large \code{R}.
#'   \code{method = "kendall"} does not admit the same closed-form
#'   vectorization and falls back to one compiled \code{cor()} call per
#'   replicate. The acceleration constant of the BCa interval is obtained
#'   from a closed-form leave-one-out (jackknife) formula (Pearson /
#'   Spearman) evaluated in O(n).
#'
#' @examples
#' x <- rnorm(40); y <- x * 0.5 + rnorm(40)
#' boot.cor.test(x, y)
#' boot.cor.test(x, y, method = "spearman")
#'
#' @seealso \code{\link[stats]{cor.test}}
#' @export
boot.cor.test <- function(x, y,
                           method = c("pearson", "kendall", "spearman"),
                           alternative = c("two.sided", "less", "greater"),
                           conf.level = 0.95, R = 10000, seed = 123, ...) {
  if (!is.null(seed)) set.seed(seed)
  method <- match.arg(method)
  alternative <- match.arg(alternative)

  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  if (length(x) != length(y)) stop("'x' and 'y' must have the same length")
  ok <- stats::complete.cases(x, y)
  x <- as.numeric(x[ok]); y <- as.numeric(y[ok])
  n <- length(x)
  if (n < 3) stop("not enough finite observations")

  if (method == "pearson") {
    theta_hat  <- stats::cor(x, y)
    theta_boot <- .boot_cor_pearson(x, y, R)
    theta_loo  <- .jack_cor_pearson(x, y)
    method_name <- "Pearson"
  } else if (method == "spearman") {
    theta_hat  <- stats::cor(x, y, method = "spearman")
    theta_boot <- .boot_cor_spearman(x, y, R)
    rx <- rank(x); ry <- rank(y)
    theta_loo  <- .jack_cor_pearson(rx, ry)
    method_name <- "Spearman"
  } else {
    theta_hat  <- stats::cor(x, y, method = "kendall")
    theta_boot <- .boot_cor_kendall(x, y, R)
    ## jackknife via leave-one-out kendall tau (loop; kendall itself is
    ## already the slow path and n is typically << R)
    theta_loo <- vapply(seq_len(n), function(i) {
      stats::cor(x[-i], y[-i], method = "kendall")
    }, numeric(1))
    method_name <- "Kendall's rank correlation tau"
  }

  ci <- .bca_ci(theta_boot, theta_hat, theta_loo, conf.level)
  a  <- attr(ci, "a")
  p.value <- .bca_pvalue(0, theta_boot, theta_hat, a, alternative)

  if (alternative == "less") {
    ci <- c(-1, ci[2])
  } else if (alternative == "greater") {
    ci <- c(ci[1], 1)
  }
  ci <- pmin(pmax(as.numeric(ci), -1), 1)

  ## report a t-like statistic on the same scale as stats::cor.test for
  ## Pearson (informational only; the p-value itself is the BCa inversion)
  if (method == "pearson") {
    tstat <- theta_hat * sqrt((n - 2) / (1 - theta_hat^2))
    names(tstat) <- "t"
    df <- n - 2
    names(df) <- "df"
    statistic <- tstat
    parameter <- df
  } else if (method == "spearman") {
    statistic <- c(S = (1 - theta_hat) * (n^3 - n) / 6)
    parameter <- NULL
  } else {
    statistic <- c(T = if (n <= 1000) {
      sum(outer(x, x, "<") & outer(y, y, "<"))
    } else {
      NA_real_  # avoid an n x n memory blow-up for large samples; the
                # BCa p-value above does not depend on this display statistic
    })
    parameter <- NULL
  }

  estimate <- theta_hat
  names(estimate) <- if (method == "kendall") "tau" else if (method == "spearman") "rho" else "cor"

  null.value <- 0
  names(null.value) <- if (method == "kendall") "tau" else if (method == "spearman") "rho" else "correlation"

  out <- list(
    statistic   = statistic,
    parameter   = parameter,
    p.value     = p.value,
    estimate    = estimate,
    null.value  = null.value,
    alternative = alternative,
    method      = paste0(method_name, "'s product-moment correlation bootstrap test",
                          "\n(BCa CI, CI-inversion p-value)"),
    data.name   = dname,
    conf.int    = structure(ci, conf.level = conf.level),
    R           = R
  )
  class(out) <- "htest"
  out
}
