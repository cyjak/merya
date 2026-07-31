#' Generalized Linear Regression with BCa Bootstrap Confidence Intervals
#'
#' Fits a generalized linear model exactly as \code{\link[stats]{glm}}
#' does (same formula interface, same fitted \code{"glm"} object), but
#' additionally performs a case-resampling bootstrap and attaches
#' bias-corrected and accelerated (BCa) confidence intervals and
#' CI-inversion p-values, retrievable via \code{summary()}. By default the
#' bootstrapped quantity is the model coefficients (as before); \code{effect}
#' can instead request odds ratios, risk ratios/differences, or the
#' population attributable fraction -- see Details.
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
#' @param effect which quantity to bootstrap and report:
#'   \describe{
#'     \item{\code{"coef"}}{(default) the ordinary model coefficients.}
#'     \item{\code{"OR"}}{odds ratios: the exact same bootstrap/jackknife
#'       coefficient replicates as \code{"coef"}, simply exponentiated
#'       (estimate and CI); the p-value is untouched. Requires
#'       \code{family = binomial(link = "logit")}.}
#'     \item{\code{"RR"}}{risk ratio. If \code{exposure} is \code{NULL},
#'       reports the *conditional* risk ratio for every predictor as the
#'       exponentiated coefficient (again the same replicates as
#'       \code{"coef"}, just exponentiated, p-value untouched); this
#'       requires a log link (any family, e.g. \code{binomial(link =
#'       "log")} or \code{poisson(link = "log")}). If \code{exposure} is
#'       given, instead reports the *marginal* risk ratio for that one
#'       binary predictor via g-computation (as for \code{"RD"} below);
#'       this requires \code{family = binomial()} (any link).}
#'     \item{\code{"RD"}}{the marginal risk difference for
#'       \code{exposure} via g-computation; requires \code{family =
#'       binomial()} and a non-\code{NULL} \code{exposure}.}
#'     \item{\code{"PAF"}}{the population attributable fraction for every
#'       eligible predictor term via g-computation; requires \code{family
#'       = binomial()}. \code{exposure} is not used.}
#'   }
#' @param exposure the name (character string) of a single binary
#'   predictor in \code{formula}; required for \code{effect = "RD"},
#'   optional for \code{effect = "RR"} (see above), and unused otherwise.
#'   Must correspond to a two-level factor or a 0/1-coded main-effect
#'   term (no interactions); see Details.
#' @param ... further arguments passed to \code{\link[stats]{glm}} for
#'   the single full-data fit (e.g. \code{contrasts}); does not affect
#'   the bootstrap replicates.
#'
#' @return An object of class \code{c("boot.glm", "glm", "lm")}: identical
#'   to what \code{\link[stats]{glm}} returns, with an additional
#'   \code{boot} element (the bootstrap replicates of whichever quantity
#'   \code{effect} selected, its BCa CI(s), CI-inversion p-value(s),
#'   \code{R}, and \code{effect} itself). All standard \code{glm} methods
#'   keep working unchanged, regardless of \code{effect}; use
#'   \code{summary()} for the bootstrap inference table (columns
#'   \code{Estimate}, \code{CI lower}, \code{CI upper}, and
#'   \code{Pr(>|z|)}; no separate bootstrap standard error column).
#'
#' @details \strong{\code{effect = "coef"}} (default): each bootstrap
#'   replicate is refit with a lean, hand-rolled IRLS solver (not
#'   \code{stats::glm.fit()}), warm-started from the full-data MLE
#'   coefficients, converging in a handful of iterations. The
#'   acceleration constant for the BCa interval is approximated from a
#'   one-step Newton (infinitesimal jackknife) update using the converged
#'   IRLS working weights, avoiding n full leave-one-out refits.
#'
#'   \strong{\code{effect = "OR"}} and the conditional-risk-ratio form of
#'   \strong{\code{effect = "RR"}} (\code{exposure = NULL}): these reuse
#'   the identical bootstrap/jackknife coefficient replicates as
#'   \code{"coef"} -- no separate resampling is done -- and simply report
#'   \code{exp(estimate)}/\code{exp(CI)} instead of the raw coefficient
#'   and its CI. Because the BCa interval's endpoints are specific
#'   (bias-and-acceleration-corrected) quantiles of the bootstrap
#'   distribution, and \code{exp()} is monotonic, exponentiating those
#'   endpoints gives the exact corresponding BCa interval on the
#'   OR/RR scale -- there is no need to (and this does not) rebuild the
#'   BCa calculation from an exponentiated bootstrap distribution. The
#'   p-value is left as computed on the coefficient scale (testing
#'   coefficient = 0, equivalently OR/RR = 1), since exponentiation does
#'   not change which side of the null a replicate falls on.
#'
#'   \strong{The g-computation form of \code{effect = "RR"}}
#'   (\code{exposure} given) \strong{, \code{effect = "RD"}, and
#'   \code{effect = "PAF"}}: all three use g-computation
#'   ("standardization"). For \code{"RR"}/\code{"RD"}, every row's
#'   \code{exposure} is set to 1 for everyone, then to 0 for everyone
#'   (leaving every other covariate at its observed value), giving
#'   population-averaged risks \eqn{R_1 = \mathrm{mean}(\hat p_i \mid
#'   \mathrm{exposure}_i = 1)} and \eqn{R_0 = \mathrm{mean}(\hat p_i \mid
#'   \mathrm{exposure}_i = 0)}; the marginal risk ratio is
#'   \eqn{R_1 / R_0} and the marginal risk difference is
#'   \eqn{R_1 - R_0}. For \code{"PAF"}, the true (observed) prevalence
#'   \eqn{P} is compared, for every eligible predictor term \eqn{i} in
#'   turn, against the counterfactual prevalence \eqn{P_i} obtained by
#'   setting that one term to its reference category for everyone (all
#'   of its design-matrix dummy columns to 0) if it is categorical, or to
#'   the sample mean for everyone if it is continuous, again leaving all
#'   other covariates as observed; the population attributable fraction
#'   is \eqn{\mathrm{PAF}_i = (P - P_i) / P}. Multi-level categorical
#'   predictors contribute a single row (one counterfactual scenario:
#'   everyone at the reference level), and only simple, untransformed,
#'   non-interaction predictor terms are supported for \code{"PAF"} (a
#'   term must match a raw column name in \code{data}; interactions or
#'   transformations such as \code{poly()}/\code{log()} raise an error).
#'
#'   For all three g-computation cases, each bootstrap replicate
#'   resamples whole rows (outcome and covariates together), refits the
#'   same lean IRLS solver used for \code{"coef"}, and repeats the
#'   standardization step within that replicate's own resampled
#'   covariate distribution and refit coefficients -- the standard
#'   nonparametric bootstrap for a g-computed effect. The BCa
#'   acceleration constant uses a closed-form (no-refit) leave-one-out
#'   approximation built from the same one-step Newton update used for
#'   \code{"coef"} (for \code{"PAF"}'s continuous-predictor recipes, the
#'   reference value used in this jackknife step is fixed at the
#'   full-data mean rather than recomputed per leave-one-out pseudo-
#'   replicate -- a small approximation that keeps the computation fully
#'   vectorized). The p-value is obtained by CI inversion against the
#'   natural null value for each scale (1 for the risk ratio, 0 for the
#'   risk difference or the attributable fraction). \code{exposure} must
#'   resolve to exactly one design-matrix column (a plain 0/1
#'   numeric/logical variable, or a two-level factor with the usual
#'   treatment contrasts) -- multi-level factors, interactions involving
#'   \code{exposure}, or non-0/1-coded columns are not supported and
#'   raise an error.
#'
#' @examples
#' fit <- boot.glm(am ~ wt + hp, data = mtcars, family = binomial(), R = 500)
#' summary(fit)
#'
#' ## odds ratios (requires a logit link)
#' summary(boot.glm(am ~ wt + hp, data = mtcars, family = binomial(),
#'                   R = 500, effect = "OR"))
#'
#' ## conditional risk ratios for every predictor (requires a log link,
#' ## but works for any family -- e.g. Poisson, not just binomial)
#' summary(boot.glm(carb ~ wt, data = mtcars, family = poisson(link = "log"),
#'                   R = 500, effect = "RR"))
#'
#' ## marginal risk ratio / risk difference for a binary predictor,
#' ## via g-computation
#' d <- mtcars; d$vs <- factor(d$vs)
#' rr_fit <- boot.glm(am ~ vs + wt, data = d, family = binomial(),
#'                     R = 500, effect = "RR", exposure = "vs")
#' summary(rr_fit)
#'
#' ## population attributable fraction for every eligible predictor
#' summary(boot.glm(am ~ vs + wt, data = d, family = binomial(),
#'                   R = 500, effect = "PAF"))
#'
#' @seealso \code{\link[stats]{glm}}
#' @export
boot.glm <- function(formula, family = stats::gaussian(), data, subset,
                      weights, na.action, offset,
                      conf.level = 0.95, R = 10000,
                      irls.maxit = 25L, irls.tol = 1e-8,
                      effect = c("coef", "OR", "RR", "RD", "PAF"),
                      exposure = NULL, seed = 123, ...) {
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

  ## ---- validate the effect/family/link/exposure combination up front ----
  use_gcomp_rr_rd <- (effect == "RD") || (effect == "RR" && !is.null(exposure))
  use_direct_rr    <- (effect == "RR" && is.null(exposure))

  if (effect == "OR") {
    if (!identical(family$family, "binomial") || !identical(family$link, "logit")) {
      stop("boot.glm(): effect = \"OR\" (odds ratio) requires family = ",
           "binomial(link = \"logit\"); got family = \"", family$family,
           "\", link = \"", family$link, "\".")
    }
  }
  if (use_direct_rr) {
    if (!identical(family$link, "log")) {
      stop("boot.glm(): effect = \"RR\" without 'exposure' (conditional ",
           "risk ratio, one per predictor) requires a log link, e.g. ",
           "family = binomial(link = \"log\") or poisson(link = \"log\"); ",
           "got link = \"", family$link, "\". Specify 'exposure' instead ",
           "for the marginal risk ratio via g-computation, which works ",
           "with any link (family = binomial()).")
    }
  }
  if (use_gcomp_rr_rd) {
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
  if (effect == "PAF" && !identical(family$family, "binomial")) {
    stop("boot.glm(): effect = \"PAF\" (population attributable fraction) ",
         "is only available for family = binomial(); got family = \"",
         family$family, "\".")
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

  if (effect %in% c("coef", "OR") || use_direct_rr) {
    ## "coef", "OR", and the direct (no-exposure) form of "RR" all bootstrap
    ## the exact same coefficient replicates; only the final reporting scale
    ## (raw vs. exponentiated) differs, and the p-value is always computed
    ## on the coefficient scale.
    boot_coef <- .boot_glm_coef(X, yy, R, family, weights = ww, offset = oo,
                                 start = beta_hat, irls.maxit = irls.maxit,
                                 irls.tol = irls.tol)
    loo_coef <- .jack_glm_coef(X, yy, fit)
    tbl <- .coef_bca_table(boot_coef, loo_coef, beta_hat, conf.level)

    if (effect == "coef") {
      fit$boot <- list(
        effect       = "coef",
        coefficients = boot_coef,
        conf.int     = tbl$conf.int,
        p.value      = tbl$p.value,
        conf.level   = conf.level,
        R            = R
      )
    } else {
      fit$boot <- list(
        effect       = effect,
        method       = "direct",
        coefficients = boot_coef,
        estimate     = exp(beta_hat),
        conf.int     = exp(tbl$conf.int),
        p.value      = tbl$p.value,
        conf.level   = conf.level,
        R            = R
      )
    }
  } else if (use_gcomp_rr_rd) {
    gc_effect <- if (effect == "RR") "rr" else "rd"
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
    theta_hat <- if (gc_effect == "rr") R1_hat / R0_hat else R1_hat - R0_hat

    boot_stat <- .boot_glm_gcomp(X, yy, R, family, weights = ww, offset = oo,
                                  start = beta_hat, irls.maxit = irls.maxit,
                                  irls.tol = irls.tol, expo_col = expo_col,
                                  effect = gc_effect)
    loo_beta <- .jack_glm_coef(X, yy, fit)
    loo_stat <- .jack_glm_gcomp(X, loo_beta, expo_col, family$linkinv, gc_effect,
                                 offset = oo)

    ok <- stats::complete.cases(boot_stat)
    tb <- boot_stat[ok]
    if (length(tb) < 10) {
      ci <- c(NA_real_, NA_real_)
      pval_val <- NA_real_
    } else {
      out <- .bca_ci(tb, theta_hat, loo_stat, conf.level)
      ci <- as.numeric(out)
      null_val <- if (gc_effect == "rr") 1 else 0
      pval_val <- .bca_pvalue(null_val, tb, theta_hat, attr(out, "a"), "two.sided")
    }
    nm <- exposure
    ci <- matrix(ci, nrow = 1L, ncol = 2, dimnames = list(nm, c("lower", "upper")))
    pval <- setNames(pval_val, nm)
    est <- setNames(theta_hat, nm)

    fit$boot <- list(
      effect       = effect,
      method       = "gcomputation",
      exposure     = exposure,
      coefficients = matrix(boot_stat, ncol = 1L, dimnames = list(NULL, nm)),
      estimate     = est,
      conf.int     = ci,
      p.value      = pval,
      conf.level   = conf.level,
      R            = R
    )
  } else {
    ## effect == "PAF"
    recipes <- .paf_recipes(mt, data, X)
    nm <- vapply(recipes, function(r) r$name, character(1))
    k <- length(recipes)

    p_true_hat <- mean(yy)
    if (!(p_true_hat > 0)) {
      stop("boot.glm(): effect = \"PAF\" requires at least one observed ",
           "event (mean(y) > 0).")
    }
    est <- setNames(numeric(k), nm)
    for (kk in seq_len(k)) {
      rec <- recipes[[kk]]
      Xcf <- X
      if (rec$type == "categorical") Xcf[, rec$cols] <- 0
      else Xcf[, rec$cols] <- mean(X[, rec$cols])
      p_cf_hat <- mean(family$linkinv(as.vector(Xcf %*% beta_hat) + oo))
      est[kk] <- (p_true_hat - p_cf_hat) / p_true_hat
    }

    boot_stat <- .boot_glm_paf(X, yy, R, family, weights = ww, offset = oo,
                                start = beta_hat, irls.maxit = irls.maxit,
                                irls.tol = irls.tol, recipes = recipes)

    loo_beta <- .jack_glm_coef(X, yy, fit)
    loo_cf <- .jack_glm_cf_prevalence(X, loo_beta, family$linkinv, recipes,
                                       offset = oo)
    Sy <- sum(yy)
    p_true_loo <- (Sy - yy) / (n - 1)          # closed-form LOO prevalence
    loo_stat <- (p_true_loo - loo_cf) / p_true_loo   # n x k, row-recycled

    ci <- matrix(NA_real_, nrow = k, ncol = 2, dimnames = list(nm, c("lower", "upper")))
    pval <- setNames(numeric(k), nm)
    for (kk in seq_len(k)) {
      ok <- stats::complete.cases(boot_stat[, kk])
      tb <- boot_stat[ok, kk]
      if (length(tb) < 10) { ci[kk, ] <- c(NA, NA); pval[kk] <- NA; next }
      out <- .bca_ci(tb, est[kk], loo_stat[, kk], conf.level)
      ci[kk, ] <- as.numeric(out)
      pval[kk] <- .bca_pvalue(0, tb, est[kk], attr(out, "a"), "two.sided")
    }

    fit$boot <- list(
      effect       = "PAF",
      method       = "gcomputation",
      coefficients = boot_stat,
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
  s$method <- b$method
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
    OR   = "Odds ratios (exponentiated coefficients)",
    RR   = if (identical(x$method, "gcomputation"))
             sprintf("Marginal risk ratio (exposure: '%s', g-computation)", x$exposure)
           else
             "Conditional risk ratios (exponentiated coefficients)",
    RD   = sprintf("Marginal risk difference (exposure: '%s', g-computation)", x$exposure),
    PAF  = "Population attributable fraction (g-computation)",
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

#' @export
#' @method print boot.glm
print.boot.glm <- function(x, ...) {
  ## Without this method, print(x) / auto-print at the console would fall
  ## through to the inherited print.glm(), which always shows the RAW
  ## model coefficients on the linear-predictor scale (e.g. still negative
  ## log-odds) regardless of 'effect' -- so simply typing a fitted
  ## effect = "OR"/"RR" object at the console would look like it returned
  ## a negative odds/risk ratio, when in fact only summary() had been
  ## updated to show the exponentiated values. Delegating to summary()
  ## ensures the console default always matches whichever quantity
  ## 'effect' selected.
  print(summary(x, ...))
  invisible(x)
}
