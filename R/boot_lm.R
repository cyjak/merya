#' Linear Regression with BCa Bootstrap Confidence Intervals
#'
#' Fits a linear model exactly as \code{\link[stats]{lm}} does (same
#' formula interface, same fitted \code{"lm"} object), but additionally
#' performs a case-resampling bootstrap and attaches bias-corrected and
#' accelerated (BCa) confidence intervals and CI-inversion p-values,
#' retrievable via \code{summary()}. By default the bootstrapped quantity
#' is the slope coefficients (as before); set \code{effect} to instead
#' report, for each predictor, its partial correlation, partial
#' eta-squared, or eta-squared.
#'
#' @param formula an object of class \code{"formula"}, as in \code{lm}.
#' @param data a data frame containing the variables in the model.
#' @param subset,weights,na.action,offset as in \code{\link[stats]{lm}}.
#' @param conf.level confidence level for the bootstrap intervals.
#' @param R number of bootstrap replicates. Default \code{10000}.
#' @param seed integer seed used to make the bootstrap resampling
#'   reproducible; set via \code{\link[base]{set.seed}} at the start of
#'   the function. Defaults to \code{123}; pass \code{NULL} to use
#'   whatever random state is currently active (no reseeding), or any
#'   other integer for a different reproducible draw.
#' @param effect which quantity to bootstrap and report, per predictor:
#'   \code{"coef"} (default) for the ordinary slope coefficients;
#'   \code{"partial.cor"} for the partial correlation between the
#'   response and each predictor (controlling for the others);
#'   \code{"partial.eta2"} for partial eta-squared (the square of the
#'   partial correlation); or \code{"eta2"} for eta-squared (the square
#'   of the semi-partial/"part" correlation, i.e. the predictor's unique
#'   share of the total variance in the response). The intercept is
#'   never reported for the three effect-size options (only for
#'   \code{"coef"}).
#' @param ... further arguments passed to \code{\link[stats]{lm}} (e.g.
#'   \code{contrasts}).
#'
#' @return An object of class \code{c("boot.lm", "lm")}: identical to what
#'   \code{\link[stats]{lm}} returns, with an additional \code{boot}
#'   element (the bootstrap replicates of whichever quantity \code{effect}
#'   selected, its BCa CIs, CI-inversion p-values, \code{R}, and
#'   \code{effect} itself). All standard \code{lm} methods
#'   (\code{predict}, \code{fitted}, \code{residuals}, \code{plot}, ...)
#'   keep working unchanged, regardless of \code{effect}; use
#'   \code{summary()} to see the bootstrap inference table.
#'
#' @details \strong{\code{effect = "coef"}} (default): coefficients are
#'   recomputed on each bootstrap resample via a direct Cholesky solve of
#'   the normal equations (\code{crossprod}/\code{backsolve}), and the
#'   BCa acceleration constant comes from the closed-form OLS
#'   leave-one-out formula (hat-matrix diagonal) -- no refitting needed
#'   for the jackknife step.
#'
#'   \strong{\code{effect = "partial.cor"} / \code{"partial.eta2"} /
#'   \code{"eta2"}}: for predictor \eqn{j}, the partial correlation is
#'   \code{sign(t_j) * sqrt(t_j^2 / (t_j^2 + df))} where \code{t_j} is
#'   that predictor's t-statistic and \code{df} the residual degrees of
#'   freedom (both from the replicate's own refit); the semi-partial
#'   ("part") correlation additionally rescales by
#'   \code{sqrt(1 - R^2)} of that same refit. Both come directly from
#'   quantities already produced by the same per-replicate Cholesky
#'   solve used for \code{"coef"} (no separate reduced-model refits are
#'   needed), and their jackknife (for the BCa acceleration constant) is
#'   likewise closed-form, extending the same leave-one-out identities
#'   used for \code{"coef"} (leave-one-out RSS via the PRESS identity,
#'   leave-one-out \eqn{(X'X)^{-1}} via a rank-one/Sherman-Morrison
#'   update, and a closed-form leave-one-out total sum of squares) --
#'   so, exactly as for \code{"coef"}, no observation or replicate is
#'   ever refit for the jackknife step.
#'
#'   Because partial eta-squared and eta-squared are *squares* of a
#'   statistic that can be negative, their confidence intervals are
#'   built in two steps: first the BCa interval is computed on the
#'   signed (partial/semi-partial correlation) scale, exactly as for
#'   \code{"partial.cor"}; then that signed interval is converted to the
#'   squared scale by taking the smallest absolute value attained within
#'   it as the lower bound (0 whenever the signed interval straddles
#'   zero) and the largest absolute value as the upper bound, and
#'   squaring both -- rather than naively squaring the signed endpoints,
#'   which would misrepresent the interval whenever it straddles zero.
#'   The p-value (by CI inversion) is computed once, on the signed
#'   scale, since squaring is monotonic in absolute value and does not
#'   change which side of zero a replicate falls on.
#'
#' @examples
#' fit <- boot.lm(mpg ~ wt + hp, data = mtcars, R = 500)
#' summary(fit)
#'
#' ## partial correlations / effect sizes instead of slopes
#' summary(boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, effect = "partial.cor"))
#' summary(boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, effect = "eta2"))
#'
#' @seealso \code{\link[stats]{lm}}
#' @export
boot.lm <- function(formula, data, subset, weights, na.action, offset,
                     conf.level = 0.95, R = 10000,
                     effect = c("coef", "partial.cor", "partial.eta2", "eta2"),
                     seed = 123, ...) {
  if (!is.null(seed)) set.seed(seed)
  effect <- match.arg(effect)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  y <- stats::model.response(mf, "numeric")
  X <- stats::model.matrix(mt, mf)
  w <- as.vector(stats::model.weights(mf))
  off <- as.vector(stats::model.offset(mf))

  ## fit via the *original* call, verbatim, so subset/weights/na.action/
  ## offset/contrasts are all honoured exactly as stats::lm() would
  lm_call <- cl
  lm_call[[1L]] <- quote(stats::lm)
  lm_call$conf.level <- NULL
  lm_call$R <- NULL
  lm_call$effect <- NULL
  lm_call$seed <- NULL
  fit <- eval(lm_call, parent.frame())
  fit$call <- cl

  beta_hat <- stats::coef(fit)
  p <- length(beta_hat)
  n <- nrow(X)

  if (!is.null(w)) {
    sw <- sqrt(w)
    Xs <- X * sw
    ys <- y * sw
  } else {
    Xs <- X
    ys <- y
  }
  if (!is.null(off)) ys <- ys - off * if (!is.null(w)) sqrt(w) else 1

  XtX_inv <- tryCatch(chol2inv(chol(crossprod(Xs))),
                       error = function(e) MASS_ginv_fallback(crossprod(Xs)))

  if (effect == "coef") {
    boot_coef <- .boot_lm_coef(Xs, ys, R)
    loo_coef  <- .jack_lm_coef(Xs, ys, beta_hat, XtX_inv)

    ci <- matrix(NA_real_, nrow = p, ncol = 2,
                 dimnames = list(names(beta_hat), c("lower", "upper")))
    pval <- setNames(numeric(p), names(beta_hat))
    for (j in seq_len(p)) {
      ok <- stats::complete.cases(boot_coef[, j])
      tb <- boot_coef[ok, j]
      if (length(tb) < 10) {
        ci[j, ] <- c(NA, NA); pval[j] <- NA
        next
      }
      out <- .bca_ci(tb, beta_hat[j], loo_coef[, j], conf.level)
      ci[j, ] <- as.numeric(out)
      pval[j] <- .bca_pvalue(0, tb, beta_hat[j], attr(out, "a"), "two.sided")
    }

    fit$boot <- list(
      effect       = "coef",
      coefficients = boot_coef,
      conf.int     = ci,
      p.value      = pval,
      conf.level   = conf.level,
      R            = R
    )
  } else {
    nm_all <- colnames(Xs)
    has_icpt <- "(Intercept)" %in% nm_all
    eff_names <- if (has_icpt) setdiff(nm_all, "(Intercept)") else nm_all
    if (length(eff_names) == 0L) {
      stop("boot.lm(): effect = \"", effect, "\" requires at least one ",
           "non-intercept predictor.")
    }
    eff_idx <- match(eff_names, nm_all)

    ## point estimates from the ordinary lm summary (t-values, residual df,
    ## R-squared) -- 'fit' still has class "lm" only at this point, so this
    ## calls stats::summary.lm(), not summary.boot.lm()
    sfit <- summary(fit)
    tvals <- sfit$coefficients[, "t value"]
    df_resid <- fit$df.residual
    pr_hat_all <- sign(tvals) * sqrt(tvals^2 / (tvals^2 + df_resid))
    r2_hat <- sfit$r.squared
    sr_hat_all <- pr_hat_all * sqrt(max(1 - r2_hat, 0))

    stat_kind <- if (effect == "eta2") "eta2" else "pcor"
    theta_hat_all <- if (stat_kind == "eta2") sr_hat_all else pr_hat_all

    boot_stat <- .boot_lm_effect(Xs, ys, R, stat_kind)
    loo_stat  <- .jack_lm_effect(Xs, ys, beta_hat, XtX_inv, stat_kind)

    ci  <- matrix(NA_real_, nrow = length(eff_names), ncol = 2,
                  dimnames = list(eff_names, c("lower", "upper")))
    pval    <- setNames(numeric(length(eff_names)), eff_names)
    est     <- setNames(numeric(length(eff_names)), eff_names)

    for (k in seq_along(eff_idx)) {
      j <- eff_idx[k]
      ok <- stats::complete.cases(boot_stat[, j])
      tb <- boot_stat[ok, j]
      if (length(tb) < 10) {
        ci[k, ] <- c(NA, NA); pval[k] <- NA
        est[k] <- theta_hat_all[j]
        next
      }
      out <- .bca_ci(tb, theta_hat_all[j], loo_stat[, j], conf.level)
      pval[k] <- .bca_pvalue(0, tb, theta_hat_all[j], attr(out, "a"), "two.sided")

      if (effect == "partial.cor") {
        ci[k, ] <- as.numeric(out)
        est[k] <- theta_hat_all[j]
      } else {
        ## partial.eta2 / eta2: square the signed CI endpoints (see
        ## .signed_ci_to_squared) and report the squared point estimate
        ci[k, ] <- .signed_ci_to_squared(as.numeric(out))
        est[k] <- theta_hat_all[j]^2
      }
    }

    fit$boot <- list(
      effect       = effect,
      coefficients = boot_stat[, eff_idx, drop = FALSE],
      estimate     = est,
      conf.int     = ci,
      p.value      = pval,
      conf.level   = conf.level,
      R            = R
    )
  }

  class(fit) <- c("boot.lm", class(fit))
  fit
}

#' @export
#' @method summary boot.lm
summary.boot.lm <- function(object, ...) {
  s <- NextMethod()
  b <- object$boot
  effect <- if (is.null(b$effect)) "coef" else b$effect

  if (effect == "coef") {
    cf <- s$coefficients
    est <- stats::coef(object)[rownames(cf)]
    new_cf <- cbind(
      Estimate   = est,
      `CI lower` = b$conf.int[rownames(cf), "lower"],
      `CI upper` = b$conf.int[rownames(cf), "upper"],
      `Pr(>|z|)` = b$p.value[rownames(cf)]
    )
  } else {
    nm <- names(b$estimate)
    new_cf <- cbind(
      Estimate   = b$estimate[nm],
      `CI lower` = b$conf.int[nm, "lower"],
      `CI upper` = b$conf.int[nm, "upper"],
      `Pr(>|z|)` = b$p.value[nm]
    )
  }
  s$coefficients <- new_cf
  s$conf.level <- b$conf.level
  s$R <- b$R
  s$effect <- effect
  class(s) <- c("summary.boot.lm", class(s))
  s
}

#' @export
#' @method print summary.boot.lm
print.summary.boot.lm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n")
  print(x$call)
  eff_label <- switch(x$effect,
    coef         = "Coefficients",
    partial.cor  = "Partial correlations",
    partial.eta2 = "Partial eta-squared",
    eta2         = "Eta-squared",
    "Coefficients"
  )
  cat(sprintf(
    "\n%s (BCa bootstrap, R = %d, %.0f%% CI; p-values via CI inversion):\n",
    eff_label, x$R, 100 * x$conf.level))
  stats::printCoefmat(x$coefficients, digits = digits, has.Pvalue = TRUE,
                       cs.ind = 1L, tst.ind = integer(0), P.values = TRUE,
                       signif.stars = getOption("show.signif.stars", TRUE))
  cat(sprintf(
    "\nResidual standard error: %s on %d degrees of freedom\n",
    format(signif(x$sigma, digits)), x$df[2L]))
  cat(sprintf("Multiple R-squared: %s,  Adjusted R-squared: %s\n",
              format(signif(x$r.squared, digits)),
              format(signif(x$adj.r.squared, digits))))
  invisible(x)
}
