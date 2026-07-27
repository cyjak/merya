#' Generalized Linear Regression with BCa Bootstrap Confidence Intervals
#'
#' Fits a generalized linear model exactly as \code{\link[stats]{glm}}
#' does (same formula interface, same fitted \code{"glm"} object), but
#' additionally performs a case-resampling bootstrap and attaches
#' bias-corrected and accelerated (BCa) confidence intervals and
#' CI-inversion p-values, retrievable via \code{summary()}. By default the
#' bootstrapped quantity is the model coefficients (as before); for a
#' binomial model, set \code{effect} to instead report a marginal
#' (population-averaged) risk ratio or risk difference for a given binary
#' predictor, obtained by g-computation.
#'
#' @param formula an object of class \code{"formula"}, as in \code{glm}.
#' @param family as in \code{\link[stats]{glm}} (e.g. \code{binomial()},
#'   \code{poisson()}, \code{gaussian()}).
#' @param data a data frame containing the variables in the model.
#' @param subset,weights,na.action,offset as in \code{\link[stats]{glm}}.
#' @param conf.level confidence level for the bootstrap intervals.
#' @param R number of bootstrap replicates. Default \code{10000}.
#' @param seed integer seed used to make the bootstrap resampling
#'   reproducible; set via \code{\link[base]{set.seed}} at the start of
#'   the function. Defaults to \code{123}; pass \code{NULL} to use
#'   whatever random state is currently active (no reseeding), or any
#'   other integer for a different reproducible draw.
#' @param irls.maxit maximum number of IRLS (Newton) iterations per
#'   bootstrap replicate. Default \code{25}; in practice, because each
#'   replicate is warm-started from the full-data MLE, convergence is
#'   typically reached in a handful of iterations.
#' @param irls.tol relative convergence tolerance on the coefficient
#'   update between IRLS iterations. Default \code{1e-8}.
#' @param effect which quantity to bootstrap and report: \code{"coef"}
#'   (default) for the ordinary model coefficients; \code{"rr"} for the
#'   marginal (population-averaged) risk ratio associated with
#'   \code{exposure}, obtained by g-computation; or \code{"rd"} for the
#'   corresponding marginal risk difference. \code{"rr"}/\code{"rd"} are
#'   only available when \code{family} is binomial.
#' @param exposure the name (character string) of a single binary
#'   predictor in \code{formula} for which the marginal risk ratio/
#'   difference is computed; required (and only used) when
#'   \code{effect} is \code{"rr"} or \code{"rd"}. Must correspond to a
#'   two-level factor or a 0/1-coded main-effect term (no interactions);
#'   see Details.
#' @param ... further arguments passed to \code{\link[stats]{glm}} for
#'   the single full-data fit (e.g. \code{contrasts}); does not affect
#'   the bootstrap replicates.
#'
#' @return An object of class \code{c("boot.glm", "glm", "lm")}: identical
#'   to what \code{\link[stats]{glm}} returns, with an additional
#'   \code{boot} element (the bootstrap replicates of whichever quantity
#'   \code{effect} selected, its BCa CI, CI-inversion p-value, \code{R},
#'   and \code{effect} itself). All standard \code{glm} methods keep
#'   working unchanged, regardless of \code{effect}; use \code{summary()}
#'   for the bootstrap inference table.
#'
#' @details \strong{\code{effect = "coef"}} (default): each bootstrap
#'   replicate is refit with a lean, hand-rolled IRLS solver (not
#'   \code{stats::glm.fit()}), warm-started from the full-data MLE
#'   coefficients. \code{stats::glm.fit()} is comparatively costly to call
#'   \code{R} times because every call revalidates inputs, tracks
#'   deviance with step-halving checks at every iteration, solves via a
#'   full rank-revealing QR decomposition, and builds a complete result
#'   object -- all overhead that is unnecessary when only the coefficient
#'   vector of each replicate is needed. The internal solver instead:
#'   warm-starts from the MLE (so convergence typically takes 2-4
#'   iterations rather than glm.fit's ~6-10 from a cold start), solves
#'   the weighted normal equations with a Cholesky factorization and two
#'   backsolves, and checks convergence via the size of the coefficient
#'   update rather than recomputing the deviance. The special case
#'   \code{family = gaussian(link = "identity")} skips iteration entirely
#'   and solves the (weighted) normal equations exactly in one step. The
#'   acceleration constant for the BCa interval is approximated from a
#'   one-step Newton (infinitesimal jackknife) update using the converged
#'   IRLS working weights, avoiding n full leave-one-out refits.
#'
#'   \strong{\code{effect = "rr"} / \code{"rd"}} (binomial only):
#'   g-computation ("standardization") estimates the marginal effect of
#'   \code{exposure} by predicting, for every row of the data, the
#'   response probability that model would give if that row's exposure
#'   were set to 1, and again as if it were set to 0 -- leaving every
#'   other covariate at its observed value -- and averaging each set of
#'   predictions over the whole sample: \eqn{R_1 = \mathrm{mean}(\hat p_i
#'   \mid \mathrm{exposure}_i = 1)}, \eqn{R_0 = \mathrm{mean}(\hat p_i
#'   \mid \mathrm{exposure}_i = 0)}. The marginal risk ratio is
#'   \eqn{R_1 / R_0} and the marginal risk difference is \eqn{R_1 - R_0}.
#'   Each bootstrap replicate resamples whole rows (outcome and
#'   covariates together), refits the same lean IRLS solver used for
#'   \code{"coef"}, and then repeats the standardization step within
#'   that replicate's own resampled covariate distribution -- the
#'   standard nonparametric bootstrap for a g-computed effect. The BCa
#'   acceleration constant uses a closed-form (no-refit) leave-one-out
#'   approximation built from the same one-step Newton update used for
#'   \code{"coef"}, applied to the g-computation formula above. The
#'   p-value is obtained by CI inversion against the natural null value
#'   for each scale (1 for the risk ratio, 0 for the risk difference).
#'   \code{exposure} must resolve to exactly one design-matrix column
#'   (a plain 0/1 numeric/logical variable, or a two-level factor with
#'   the usual treatment contrasts) -- multi-level factors, interactions
#'   involving \code{exposure}, or non-0/1-coded columns are not
#'   supported and raise an error.
#'
#' @examples
#' fit <- boot.glm(am ~ wt + hp, data = mtcars, family = binomial(), R = 500)
#' summary(fit)
#'
#' ## marginal risk ratio / risk difference for a binary predictor,
#' ## via g-computation
#' d <- mtcars; d$vs <- factor(d$vs)
#' rr_fit <- boot.glm(am ~ vs + wt, data = d, family = binomial(),
#'                     R = 500, effect = "rr", exposure = "vs")
#' summary(rr_fit)
#'
#' @seealso \code{\link[stats]{glm}}
#' @export
boot.glm <- function(formula, family = stats::gaussian(), data, subset,
                      weights, na.action, offset,
                      conf.level = 0.95, R = 10000,
                      irls.maxit = 25L, irls.tol = 1e-8,
                      effect = c("coef", "rr", "rd"), exposure = NULL,
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
  y <- stats::model.response(mf, "any")
  if (is.matrix(y)) {
    stop("boot.glm() does not support matrix-style (cbind(success, failure)) ",
         "responses; please supply a 0/1 or proportion response with 'weights' ",
         "for binomial models instead.")
  }
  X <- stats::model.matrix(mt, mf)
  w <- stats::model.weights(mf)
  off <- stats::model.offset(mf)
  if (is.character(family)) family <- get(family, mode = "function")
  if (is.function(family)) family <- family()

  if (effect != "coef") {
    if (!identical(family$family, "binomial")) {
      stop("boot.glm(): effect = \"", effect, "\" (marginal risk ratio/",
           "difference by g-computation) is only available for family = ",
           "binomial(); got family = \"", family$family, "\".")
    }
    if (is.null(exposure) || !is.character(exposure) || length(exposure) != 1L) {
      stop("boot.glm(): effect = \"", effect, "\" requires 'exposure' to ",
           "be the (single) name of a binary predictor in 'formula'.")
    }
  }

  ## fit via the *original* call, verbatim, so subset/weights/na.action/
  ## offset are all honoured exactly as stats::glm() would
  glm_call <- cl
  glm_call[[1L]] <- quote(stats::glm)
  glm_call$conf.level <- NULL
  glm_call$R <- NULL
  glm_call$irls.maxit <- NULL
  glm_call$irls.tol <- NULL
  glm_call$effect <- NULL
  glm_call$exposure <- NULL
  glm_call$seed <- NULL
  glm_call$family <- family
  fit <- eval(glm_call, parent.frame())
  fit$call <- cl

  beta_hat <- stats::coef(fit)
  p <- length(beta_hat)
  n <- nrow(X)

  yy <- if (is.factor(y)) as.numeric(y) - 1 else as.numeric(y)
  ww <- if (is.null(w)) rep(1, n) else w
  oo <- if (is.null(off)) rep(0, n) else off

  if (effect == "coef") {
    boot_coef <- .boot_glm_coef(X, yy, R, family, weights = ww, offset = oo,
                                 start = beta_hat, irls.maxit = irls.maxit,
                                 irls.tol = irls.tol)
    loo_coef <- .jack_glm_coef(X, yy, fit)

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
    nm_all <- colnames(X)
    term.labels <- attr(mt, "term.labels")
    term_idx <- match(exposure, term.labels)
    if (is.na(term_idx)) {
      stop("boot.glm(): exposure = \"", exposure, "\" was not found as a ",
           "main-effect term in 'formula'. Supported term labels are: ",
           paste(term.labels, collapse = ", "), ".")
    }
    assign_vec <- attr(X, "assign")
    expo_cols <- which(assign_vec == term_idx)
    if (length(expo_cols) != 1L) {
      stop("boot.glm(): exposure = \"", exposure, "\" must correspond to ",
           "exactly one design-matrix column (a 0/1 numeric/logical ",
           "variable, or a two-level factor with treatment contrasts); ",
           "it currently maps to ", length(expo_cols), " column(s). ",
           "Multi-level factors and interactions involving 'exposure' ",
           "are not supported.")
    }
    expo_col <- expo_cols[1L]
    expo_vals <- X[, expo_col]
    if (!all(expo_vals %in% c(0, 1))) {
      stop("boot.glm(): the design-matrix column for exposure = \"", exposure,
           "\" is not coded as 0/1; g-computation requires a binary ",
           "(0/1-coded) exposure.")
    }

    ## point estimate: g-computation on the observed data with the MLE
    X1 <- X; X1[, expo_col] <- 1
    X0 <- X; X0[, expo_col] <- 0
    p1_hat <- family$linkinv(as.vector(X1 %*% beta_hat) + oo)
    p0_hat <- family$linkinv(as.vector(X0 %*% beta_hat) + oo)
    R1_hat <- mean(p1_hat); R0_hat <- mean(p0_hat)
    theta_hat <- if (effect == "rr") R1_hat / R0_hat else R1_hat - R0_hat

    boot_stat <- .boot_glm_gcomp(X, yy, R, family, weights = ww, offset = oo,
                                  start = beta_hat, irls.maxit = irls.maxit,
                                  irls.tol = irls.tol, expo_col = expo_col,
                                  effect = effect)
    loo_beta <- .jack_glm_coef(X, yy, fit)
    loo_stat <- .jack_glm_gcomp(X, loo_beta, expo_col, family$linkinv, effect,
                                 offset = oo)

    ok <- stats::complete.cases(boot_stat)
    tb <- boot_stat[ok]
    if (length(tb) < 10) {
      ci <- c(NA_real_, NA_real_)
      pval_val <- NA_real_
    } else {
      out <- .bca_ci(tb, theta_hat, loo_stat, conf.level)
      ci <- as.numeric(out)
      null_val <- if (effect == "rr") 1 else 0
      pval_val <- .bca_pvalue(null_val, tb, theta_hat, attr(out, "a"), "two.sided")
    }
    nm <- exposure
    ci <- matrix(ci, nrow = 1L, ncol = 2, dimnames = list(nm, c("lower", "upper")))
    pval <- setNames(pval_val, nm)
    est <- setNames(theta_hat, nm)

    fit$boot <- list(
      effect       = effect,
      exposure     = exposure,
      coefficients = matrix(boot_stat, ncol = 1L, dimnames = list(NULL, nm)),
      estimate     = est,
      conf.int     = ci,
      p.value      = pval,
      conf.level   = conf.level,
      R            = R
    )
  }

  class(fit) <- c("boot.glm", class(fit))
  fit
}

#' @export
#' @method summary boot.glm
summary.boot.glm <- function(object, ...) {
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
  s$exposure <- b$exposure
  class(s) <- c("summary.boot.glm", class(s))
  s
}

#' @export
#' @method print summary.boot.glm
print.summary.boot.glm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n")
  print(x$call)
  eff_label <- switch(x$effect,
    coef = "Coefficients",
    rr   = sprintf("Marginal risk ratio (exposure: '%s', g-computation)", x$exposure),
    rd   = sprintf("Marginal risk difference (exposure: '%s', g-computation)", x$exposure),
    "Coefficients"
  )
  cat(sprintf(
    "\n%s (BCa bootstrap, R = %d, %.0f%% CI; p-values via CI inversion):\n",
    eff_label, x$R, 100 * x$conf.level))
  stats::printCoefmat(x$coefficients, digits = digits, has.Pvalue = TRUE,
                       cs.ind = 1L, tst.ind = integer(0), P.values = TRUE,
                       signif.stars = getOption("show.signif.stars", TRUE))
  cat(sprintf("\n(Dispersion parameter taken to be %s)\n",
              format(signif(x$dispersion, digits))))
  cat(sprintf("Null deviance: %s on %d degrees of freedom\n",
              format(signif(x$null.deviance, digits)), x$df.null))
  cat(sprintf("Residual deviance: %s on %d degrees of freedom\n",
              format(signif(x$deviance, digits)), x$df.residual))
  cat(sprintf("AIC: %s\n", format(signif(x$aic, digits))))
  invisible(x)
}
