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
#' @param pred.r.squared logical; if \code{TRUE}, additionally compute the
#'   textbook "predicted R-squared" (\code{1 - PRESS/TSS}, leave-one-out
#'   cross-validated), reported alongside whichever \code{effect} was
#'   requested. Unlike ordinary or partial/semi-partial eta-squared, this
#'   is not bounded below by 0 -- a model that predicts worse than just
#'   the response mean legitimately has a negative predicted R-squared,
#'   and this is preserved rather than treated as bootstrap noise around
#'   a non-negative quantity. \code{FALSE} by default (not computed). See
#'   Details.
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
#'   ("part") correlation is instead \code{t_j * sqrt((1 - R^2) / df)}
#'   of that same refit. This is equivalent to \code{sign(t_j) *
#'   sqrt(R^2_full - R^2_reduced)} -- the (signed) square root of the
#'   unique drop in R-squared from omitting predictor \eqn{j} -- via the
#'   general-linear-model identity \eqn{F = t_j^2} for dropping a single
#'   predictor; it is \emph{not} the partial correlation itself further
#'   rescaled by \code{sqrt(1 - R^2)}, which only coincides with the
#'   correct value in the trivial \code{t_j -> 0} limit and otherwise
#'   understates the effect. Both come directly from quantities already
#'   produced by the same per-replicate Cholesky solve used for
#'   \code{"coef"} (no separate reduced-model refits are needed), and
#'   their jackknife (for the BCa acceleration constant) is likewise
#'   closed-form, extending the same leave-one-out identities
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
#'   \strong{\code{pred.r.squared = TRUE}}: for each observation \eqn{j},
#'   its leave-one-out predicted value is \eqn{x_j'\hat\beta_{(-j)}},
#'   where \eqn{\hat\beta_{(-j)}} is the coefficient vector from the same
#'   model refit excluding observation \eqn{j} -- obtained in closed form
#'   via the hat-matrix identity (no refitting), the same one used for
#'   \code{effect = "coef"}'s BCa acceleration constant. From these,
#'   \eqn{\mathrm{PRESS} = \sum_j (y_j - x_j'\hat\beta_{(-j)})^2} (the
#'   leave-one-out predicted error sum of squares) and predicted
#'   R-squared is \eqn{1 - \mathrm{PRESS}/\mathrm{TSS}}, where
#'   \eqn{\mathrm{TSS}} is the total sum of squares -- the standard
#'   definition of predicted R-squared (as reported by, e.g., Minitab and
#'   most regression-diagnostics references). This quantity is already on
#'   its natural, signed scale (it is not a square of anything), so it is
#'   bootstrapped directly (each replicate resamples rows, refits, and
#'   repeats the same leave-one-out-within-the-replicate PRESS/TSS
#'   calculation on its own resampled data) and given a BCa confidence
#'   interval and CI-inversion p-value with no separate sign-handling
#'   transform of any kind -- unlike \code{"partial.eta2"}/\code{"eta2"},
#'   which square a signed correlation that is known to estimate a
#'   non-negative population quantity, predicted R-squared's negative
#'   values are directly meaningful (the model predicts worse than the
#'   mean) and are not an artifact to be corrected for. The BCa
#'   acceleration constant is approximated with a fast closed-form
#'   "delete-one" leave-one-out PRESS/TSS (excluding one observation's
#'   own term from the full-sample PRESS sum, and using the same
#'   closed-form leave-one-out TSS identity as \code{"eta2"}, rather than
#'   a full nested double-jackknife), consistent with the closed-form/
#'   one-step approximations used throughout this package.
#'
#' @examples
#' fit <- boot.lm(mpg ~ wt + hp, data = mtcars, R = 500)
#' summary(fit)
#'
#' ## partial correlations / effect sizes instead of slopes
#' summary(boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, effect = "partial.cor"))
#' summary(boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, effect = "eta2"))
#'
#' ## also report predicted (leave-one-out) R-squared
#' summary(boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, pred.r.squared = TRUE))
#'
#' @seealso \code{\link[stats]{lm}}
#' @export
boot.lm <- function(formula, data, subset, weights, na.action, offset,
                     conf.level = 0.95, R = 10000,
                     effect = c("coef", "partial.cor", "partial.eta2", "eta2"),
                     pred.r.squared = FALSE, seed = 123, ...) {
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
  lm_call$pred.r.squared <- NULL
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

  ## shared by effect = "coef" and by pred.r.squared (leave-one-out OLS
  ## coefficients, closed-form via the hat-matrix identity -- no refits)
  need_loo_coef <- (effect == "coef") || isTRUE(pred.r.squared)
  loo_coef <- if (need_loo_coef) .jack_lm_coef(Xs, ys, beta_hat, XtX_inv) else NULL

  if (effect == "coef") {
    boot_coef <- .boot_lm_coef(Xs, ys, R)
    tbl <- .coef_bca_table(boot_coef, loo_coef, beta_hat, conf.level)

    fit$boot <- list(
      effect       = "coef",
      coefficients = boot_coef,
      conf.int     = tbl$conf.int,
      p.value      = tbl$p.value,
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
    # sr_i = t_i * sqrt((1-R^2)/df_resid) -- NOT pr_i * sqrt(1-R^2); see
    # .boot_lm_effect() for the derivation (F = t_i^2 for a single-df drop).
    sr_hat_all <- tvals * sqrt(pmax(1 - r2_hat, 0) / df_resid)

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

  if (isTRUE(pred.r.squared)) {
    ## Textbook "predicted R-squared": 1 - PRESS/TSS, where PRESS is the
    ## leave-one-out predicted error sum of squares. This is already on
    ## its natural (signed, unbounded-below) scale -- a model that
    ## predicts worse than the mean legitimately gives a negative value
    ## here -- so, unlike partial.eta2/eta2, no separate sign-preserving
    ## square transform is applied anywhere: the BCa interval is built
    ## directly on this statistic.
    pred_loo_full <- rowSums(Xs * loo_coef)
    has_icpt_predR <- "(Intercept)" %in% colnames(Xs)
    tss_hat <- sum((ys - if (has_icpt_predR) mean(ys) else 0)^2)
    press_hat <- sum((ys - pred_loo_full)^2)
    theta_hat_predR <- if (tss_hat > 0) 1 - press_hat / tss_hat else NA_real_

    boot_predR <- .boot_lm_predR2(Xs, ys, R)
    loo_predR  <- .jack_lm_predR2(ys, pred_loo_full, has_icpt_predR)

    ok <- stats::complete.cases(boot_predR)
    tb <- boot_predR[ok]
    if (length(tb) < 10 || is.na(theta_hat_predR)) {
      ci_predR <- c(NA_real_, NA_real_)
      pval_predR <- NA_real_
    } else {
      out <- .bca_ci(tb, theta_hat_predR, loo_predR, conf.level)
      ci_predR <- as.numeric(out)
      pval_predR <- .bca_pvalue(0, tb, theta_hat_predR, attr(out, "a"), "two.sided")
    }

    fit$boot$pred.r.squared <- list(
      estimate   = theta_hat_predR,
      conf.int   = matrix(ci_predR, nrow = 1L, dimnames = list("pred.r.squared", c("lower", "upper"))),
      p.value    = pval_predR,
      conf.level = conf.level,
      R          = R
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
  s$pred.r.squared <- b$pred.r.squared
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
  if (!is.null(x$pred.r.squared)) {
    pr <- x$pred.r.squared
    cat(sprintf(
      "Predicted R-squared (leave-one-out, BCa bootstrap, R = %d, %.0f%% CI; ",
      pr$R, 100 * pr$conf.level))
    cat(sprintf(
      "p-value via CI inversion):\n  %s  [%s, %s]  p = %s\n",
      format(signif(pr$estimate, digits)),
      format(signif(pr$conf.int[1L, "lower"], digits)),
      format(signif(pr$conf.int[1L, "upper"], digits)),
      format(signif(pr$p.value, digits))))
  }
  invisible(x)
}

#' @export
#' @method print boot.lm
print.boot.lm <- function(x, ...) {
  ## Without this method, print(x) / auto-print at the console would fall
  ## through to the inherited print.lm(), which always shows the RAW,
  ## un-transformed model coefficients regardless of 'effect' -- e.g. still
  ## negative log-odds-scale-like slopes even when effect = "partial.eta2"
  ## was requested. Delegating to summary() ensures the console default
  ## always matches whichever quantity 'effect' selected.
  print(summary(x, ...))
  invisible(x)
}
