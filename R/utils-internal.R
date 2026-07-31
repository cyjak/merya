## Internal engine for BCa bootstrap confidence intervals and
## confidence-interval-inversion p-values. Not exported.
##
## Design goals:
##  - base R only (stats, which ships with every R installation)
##  - vectorized resampling (single index matrix, matrix algebra) instead of
##    per-replicate loops wherever the statistic admits a closed form
##  - closed-form BCa inversion for p-values (no bisection / root-finding)

#' @keywords internal
#' @noRd
.boot_idx <- function(n, R) {
  ## one (n x R) matrix of bootstrap resample indices, drawn in a single call
  matrix(sample.int(n, size = n * R, replace = TRUE), nrow = n, ncol = R)
}

#' Normal quantile helper (avoids repeated qnorm() dispatch overhead)
#' @keywords internal
#' @noRd
.qnorm <- function(p) stats::qnorm(p)
.pnorm <- function(q) stats::pnorm(q)

#' Compute BCa acceleration constant from jackknife (leave-one-out) values
#'
#' @param theta_loo numeric vector of leave-one-out statistic values (length n)
#' @keywords internal
#' @noRd
.bca_accel <- function(theta_loo) {
  theta_dot <- mean(theta_loo)
  d <- theta_dot - theta_loo
  num <- sum(d^3)
  den <- 6 * (sum(d^2))^1.5
  if (den <= .Machine$double.eps) return(0)
  num / den
}

#' Bias-correction constant z0 from bootstrap replicates and the observed stat
#' @keywords internal
#' @noRd
.bca_z0 <- function(theta_boot, theta_hat) {
  p0 <- mean(theta_boot < theta_hat) + 0.5 * mean(theta_boot == theta_hat)
  p0 <- min(max(p0, 1 / (length(theta_boot) + 1)), length(theta_boot) / (length(theta_boot) + 1))
  .qnorm(p0)
}

#' BCa confidence interval from bootstrap replicates
#'
#' @param theta_boot numeric vector of bootstrap replicate statistics
#' @param theta_hat observed statistic
#' @param theta_loo leave-one-out (jackknife) statistic values used for
#'   the acceleration constant. May also be a precomputed acceleration value
#'   if `accel_only = TRUE`.
#' @param conf.level confidence level
#' @return named numeric vector c(lower, upper), plus attributes z0, a
#' @keywords internal
#' @noRd
.bca_ci <- function(theta_boot, theta_hat, theta_loo, conf.level = 0.95,
                     accel_only = FALSE) {
  z0 <- .bca_z0(theta_boot, theta_hat)
  a  <- if (accel_only) theta_loo else .bca_accel(theta_loo)

  alpha <- (1 - conf.level) / 2
  z_alpha <- .qnorm(c(alpha, 1 - alpha))

  adj <- z0 + (z0 + z_alpha) / (1 - a * (z0 + z_alpha))
  probs <- .pnorm(adj)
  probs <- pmin(pmax(probs, 0), 1)

  ci <- stats::quantile(theta_boot, probs = probs, type = 7, names = FALSE)
  out <- sort(ci)
  attr(out, "z0") <- z0
  attr(out, "a") <- a
  out
}

#' CI-inversion p-value from the BCa bootstrap distribution (closed form)
#'
#' Uses the analytic inverse of the BCa transformation: rather than
#' bisecting on alpha to find where the BCa CI bound crosses the null value
#' (theta0), we invert the transform directly, which is O(1) given z0 and a.
#'
#' @param theta0 the null-hypothesis value of the statistic
#' @param theta_boot bootstrap replicate statistics
#' @param theta_hat observed statistic
#' @param a precomputed acceleration constant
#' @param alternative "two.sided", "less", or "greater"
#' @keywords internal
#' @noRd
.bca_pvalue <- function(theta0, theta_boot, theta_hat, a,
                         alternative = "two.sided") {
  n_boot <- length(theta_boot)
  z0 <- .bca_z0(theta_boot, theta_hat)

  ## empirical proportion of bootstrap replicates at/below theta0
  p0 <- mean(theta_boot <= theta0)
  p0 <- min(max(p0, 1 / (n_boot + 1)), n_boot / (n_boot + 1))
  u <- .qnorm(p0)

  v <- u - z0
  denom <- 1 + v * a
  if (abs(denom) < .Machine$double.eps) denom <- sign(denom) * .Machine$double.eps
  z_alpha <- v / denom - z0

  alpha_one_sided <- .pnorm(z_alpha)   # P(theta_true <= theta0)-type tail prob

  p <- switch(alternative,
              two.sided = 2 * min(alpha_one_sided, 1 - alpha_one_sided),
              less      = 1 - alpha_one_sided,
              greater   = alpha_one_sided,
              stop("invalid 'alternative'"))
  min(max(p, 0), 1)
}

#' Fast vectorized bootstrap mean (single sample) for use in t-tests
#' @keywords internal
#' @noRd
.boot_mean <- function(x, R) {
  n <- length(x)
  idx <- .boot_idx(n, R)
  colMeans(matrix(x[idx], nrow = n, ncol = R))
}

#' Fast vectorized bootstrap mean-difference (two independent samples)
#' @keywords internal
#' @noRd
.boot_mean_diff <- function(x, y, R) {
  .boot_mean(x, R) - .boot_mean(y, R)
}

#' Fast vectorized jackknife means (leave-one-out), used for acceleration
#' @keywords internal
#' @noRd
.jack_mean <- function(x) {
  n <- length(x)
  s <- sum(x)
  (s - x) / (n - 1)
}

#' Fast vectorized bootstrap Pearson correlation using sufficient statistics
#' (sums, sums of squares/products) computed columnwise over the resample
#' index matrix -- avoids any per-replicate loop.
#' @keywords internal
#' @noRd
.boot_cor_pearson <- function(x, y, R) {
  n <- length(x)
  idx <- .boot_idx(n, R)
  Xb <- matrix(x[idx], nrow = n, ncol = R)
  Yb <- matrix(y[idx], nrow = n, ncol = R)
  sx  <- colSums(Xb);  sy  <- colSums(Yb)
  sxx <- colSums(Xb * Xb); syy <- colSums(Yb * Yb); sxy <- colSums(Xb * Yb)
  num <- n * sxy - sx * sy
  den <- sqrt((n * sxx - sx^2) * (n * syy - sy^2))
  ifelse(den <= .Machine$double.eps, 0, num / den)
}

#' Jackknife Pearson correlations (leave-one-out), closed form, O(n)
#' @keywords internal
#' @noRd
.jack_cor_pearson <- function(x, y) {
  n <- length(x)
  sx <- sum(x); sy <- sum(y)
  sxx <- sum(x * x); syy <- sum(y * y); sxy <- sum(x * y)
  ## leave-one-out sums
  sx_i  <- sx  - x
  sy_i  <- sy  - y
  sxx_i <- sxx - x * x
  syy_i <- syy - y * y
  sxy_i <- sxy - x * y
  m <- n - 1
  num <- m * sxy_i - sx_i * sy_i
  den <- sqrt((m * sxx_i - sx_i^2) * (m * syy_i - sy_i^2))
  ifelse(den <= .Machine$double.eps, 0, num / den)
}

#' Bootstrap Spearman correlation (rank-based). Ranks are recomputed on each
#' resample (ties can change under resampling), then the vectorized Pearson
#' routine above is reused on the ranks.
#' @keywords internal
#' @noRd
.boot_cor_spearman <- function(x, y, R) {
  n <- length(x)
  idx <- .boot_idx(n, R)
  Xb <- matrix(x[idx], nrow = n, ncol = R)
  Yb <- matrix(y[idx], nrow = n, ncol = R)
  Rx <- apply(Xb, 2, rank)
  Ry <- apply(Yb, 2, rank)
  sx  <- colSums(Rx);  sy  <- colSums(Ry)
  sxx <- colSums(Rx * Rx); syy <- colSums(Ry * Ry); sxy <- colSums(Rx * Ry)
  num <- n * sxy - sx * sy
  den <- sqrt((n * sxx - sx^2) * (n * syy - sy^2))
  ifelse(den <= .Machine$double.eps, 0, num / den)
}

#' Bootstrap Kendall's tau (no closed-form vectorization; looped but uses the
#' compiled base implementation `cor(method = "kendall")` per replicate).
#' @keywords internal
#' @noRd
.boot_cor_kendall <- function(x, y, R) {
  n <- length(x)
  idx <- .boot_idx(n, R)
  out <- numeric(R)
  for (b in seq_len(R)) {
    out[b] <- stats::cor(x[idx[, b]], y[idx[, b]], method = "kendall")
  }
  out
}

#' Fast case-resampling bootstrap for (weighted) least squares coefficients.
#'
#' Uses a pre-factorized-per-replicate crossprod/solve on the (typically
#' small) p x p system, which is the practical speed bottleneck-free
#' approach in base R (BLAS/LAPACK compiled code does the heavy lifting;
#' looping only over R replicates of small linear solves is fast).
#'
#' @param X model matrix (n x p)
#' @param y response vector (n)
#' @param R number of bootstrap replicates
#' @return (R x p) matrix of bootstrap coefficient estimates
#' @keywords internal
#' @noRd
.boot_lm_coef <- function(X, y, R) {
  n <- nrow(X); p <- ncol(X)
  idx <- .boot_idx(n, R)
  out <- matrix(NA_real_, nrow = R, ncol = p)
  for (b in seq_len(R)) {
    rows <- idx[, b]
    Xb <- X[rows, , drop = FALSE]
    yb <- y[rows]
    XtX <- crossprod(Xb)
    Xty <- crossprod(Xb, yb)
    ch <- tryCatch(chol(XtX), error = function(e) NULL)
    if (is.null(ch)) {
      out[b, ] <- tryCatch(qr.solve(Xb, yb), error = function(e) rep(NA_real_, p))
    } else {
      out[b, ] <- backsolve(ch, backsolve(ch, Xty, transpose = TRUE))
    }
  }
  out
}

#' Leverage-based (infinitesimal-jackknife) approximate leave-one-out
#' coefficients for OLS, used only to compute the BCa acceleration constant
#' quickly (O(n) given the hat-matrix diagonal, no refitting).
#'
#' beta_(i) approx beta_hat - H_ii/(1-H_ii) * (X_i' ...) -- we use the exact
#' closed-form leave-one-out update available for OLS via the hat matrix:
#'   beta_(i) = beta_hat - (1 / (1 - h_ii)) * (XtX_inv %*% x_i) * e_i
#' where e_i is the ordinary residual and h_ii the leverage.
#'
#' @keywords internal
#' @noRd
.jack_lm_coef <- function(X, y, beta_hat, XtX_inv) {
  n <- nrow(X)
  fitted <- as.vector(X %*% beta_hat)
  e <- y - fitted
  h <- rowSums((X %*% XtX_inv) * X)
  h <- pmin(h, 1 - 1e-10)
  ## (n x p) matrix of leave-one-out coefficient deviations
  XtXinvXt <- X %*% XtX_inv          # n x p  (this is (X'X)^{-1} x_i as rows)
  delta <- (e / (1 - h)) * XtXinvXt  # n x p, row i = correction for obs i
  ## beta_(i) = beta_hat - delta_i  ->  return as (n x p) matrix of LOO betas
  sweep(-delta, 2, beta_hat, "+")
}

#' Assemble BCa confidence intervals and CI-inversion p-values for a whole
#' matrix of bootstrap / leave-one-out coefficient replicates in one
#' pass, on the natural (unt transformed) coefficient scale. Shared by
#' boot.lm()'s and boot.glm()'s "coef" branch, and by boot.glm()'s "OR"
#' and direct ("no exposure", log-link) "RR" branches, which reuse the
#' exact same bootstrap/jackknife coefficient replicates and simply
#' exponentiate the resulting estimate/CI afterwards (the p-value is left
#' untouched, since exponentiation is monotonic and does not change which
#' side of the null a replicate falls on).
#'
#' @param boot_mat (R x p) matrix of bootstrap replicate coefficients
#'   (failed replicates should be NA rows)
#' @param loo_mat (n x p) matrix of leave-one-out coefficients
#' @param beta_hat named numeric vector, the full-data estimate
#' @keywords internal
#' @noRd
.coef_bca_table <- function(boot_mat, loo_mat, beta_hat, conf.level = 0.95) {
  p  <- length(beta_hat)
  nm <- names(beta_hat)
  ci   <- matrix(NA_real_, nrow = p, ncol = 2, dimnames = list(nm, c("lower", "upper")))
  pval <- stats::setNames(rep(NA_real_, p), nm)
  for (j in seq_len(p)) {
    ok <- stats::complete.cases(boot_mat[, j])
    tb <- boot_mat[ok, j]
    if (length(tb) < 10) next
    out <- .bca_ci(tb, beta_hat[j], loo_mat[, j], conf.level)
    ci[j, ] <- as.numeric(out)
    pval[j] <- .bca_pvalue(0, tb, beta_hat[j], attr(out, "a"), "two.sided")
  }
  list(conf.int = ci, p.value = pval)
}

#' Fast case-resampling bootstrap of the "predicted R" statistic used by
#' boot.lm()'s \code{pred.r.squared} option: for each replicate,
#' resample rows, refit OLS via the same Cholesky solve used elsewhere
#' (\code{.boot_lm_coef}), then -- reusing the exact leave-one-out
#' identity already used for the BCa acceleration constant
#' (\code{.jack_lm_coef}), but now applied *within* this one replicate's
#' own resampled data -- get every resampled observation's leave-one-out
#' predicted value, and return the (signed) Pearson correlation between
#' the resampled response and those leave-one-out predictions. This is
#' the quantity whose square is reported as "predicted R-squared".
#'
#' @return length-R numeric vector (NA for replicates whose refit failed
#'   or were otherwise degenerate)
#' @keywords internal
#' @noRd
.boot_lm_predR <- function(X, y, R) {
  n <- nrow(X); p <- ncol(X)
  df_resid <- n - p
  idx <- .boot_idx(n, R)
  out <- rep(NA_real_, R)
  if (df_resid <= 1L) return(out)

  for (b in seq_len(R)) {
    rows <- idx[, b]
    Xb <- X[rows, , drop = FALSE]
    yb <- y[rows]
    XtX <- crossprod(Xb)
    ch <- tryCatch(chol(XtX), error = function(e) NULL)
    if (is.null(ch)) next
    Xty <- crossprod(Xb, yb)
    beta_b <- backsolve(ch, backsolve(ch, Xty, transpose = TRUE))
    XtX_inv_b <- tryCatch(chol2inv(ch), error = function(e) NULL)
    if (is.null(XtX_inv_b)) next

    fitted_b <- as.vector(Xb %*% beta_b)
    e_b <- yb - fitted_b
    h_b <- rowSums((Xb %*% XtX_inv_b) * Xb)
    h_b <- pmin(h_b, 1 - 1e-10)
    V_b <- Xb %*% XtX_inv_b
    delta_b <- (e_b / (1 - h_b)) * V_b
    beta_loo_b <- sweep(-delta_b, 2, beta_b, "+")     # n x p, within-replicate LOO coefs
    pred_loo_b <- rowSums(Xb * beta_loo_b)            # n, within-replicate LOO predictions

    if (stats::sd(pred_loo_b) <= 0 || stats::sd(yb) <= 0) next
    out[b] <- stats::cor(yb, pred_loo_b)
  }
  out
}

#' Closed-form "delete-one-pair" leave-one-out correlation between two
#' vectors, used only for the BCa acceleration constant of
#' \code{pred.r.squared}: rather than a full nested double-jackknife (for
#' each held-out observation, refitting *and* recomputing every other
#' observation's own leave-one-out prediction on the reduced (n-1)-point
#' sample -- an O(n^2)-ish computation), this reuses the already-computed
#' full-sample leave-one-out predictions and simply excludes one
#' (response, prediction) pair at a time from the correlation formula,
#' via the standard incremental sum-of-squares/sum-of-products identity.
#' This is an approximation of the true nested jackknife (it does not
#' re-derive each remaining point's leave-one-out prediction under the
#' (n-1)-point sample), consistent with the fast, one-step/closed-form
#' jackknife approximations already used elsewhere in this package, and
#' keeps the whole computation O(n) instead of O(n^2).
#'
#' @param u,v numeric vectors of the same length (response, leave-one-out
#'   predictions)
#' @return length-n numeric vector of leave-one-out correlations
#' @keywords internal
#' @noRd
.loo_cor_exclude_pairs <- function(u, v) {
  n <- length(u)
  m <- n - 1
  Su <- sum(u); Sv <- sum(v)
  Suu <- sum(u^2); Svv <- sum(v^2); Suv <- sum(u * v)
  Su_i <- Su - u;   Sv_i <- Sv - v
  Suu_i <- Suu - u^2; Svv_i <- Svv - v^2; Suv_i <- Suv - u * v
  cov_i  <- Suv_i - Su_i * Sv_i / m
  varu_i <- Suu_i - Su_i^2 / m
  varv_i <- Svv_i - Sv_i^2 / m
  denom <- sqrt(varu_i * varv_i)
  ifelse(denom > 0, cov_i / denom, NA_real_)
}

#' Given a BCa confidence interval computed on a *signed* statistic (one
#' that can be negative -- a partial or semi-partial correlation), derive
#' the confidence interval for its *square* (partial eta-squared /
#' eta-squared) without re-running the bootstrap: naively squaring the
#' endpoints (lower^2, upper^2) is wrong whenever the signed interval
#' spans zero, since the squared quantity can then get arbitrarily close
#' to 0 even though neither endpoint is 0. The correct endpoints are the
#' smallest and largest *absolute* value attainable within the signed
#' interval, each then squared: 0 whenever the interval contains 0 (for
#' the lower bound), and max(|lower|, |upper|) otherwise/always (for the
#' upper bound).
#'
#' @param ci length-2 numeric vector c(lower, upper) on the signed scale
#' @return length-2 numeric vector c(lower, upper) on the squared scale
#' @keywords internal
#' @noRd
.signed_ci_to_squared <- function(ci) {
  lower <- ci[1]; upper <- ci[2]
  if (anyNA(c(lower, upper))) return(c(NA_real_, NA_real_))
  min_abs <- if (lower <= 0 && upper >= 0) 0 else min(abs(lower), abs(upper))
  max_abs <- max(abs(lower), abs(upper))
  c(min_abs^2, max_abs^2)
}

#' Fast case-resampling bootstrap of a *signed* per-predictor effect-size
#' statistic for OLS: the partial correlation (used as-is for
#' \code{effect = "partial.cor"}, and as the pre-squaring basis for
#' \code{effect = "partial.eta2"} since partial-eta2_j = partial-cor_j^2),
#' or the semi-partial ("part") correlation (the pre-squaring basis for
#' \code{effect = "eta2"}, since eta2_j = semi-partial-cor_j^2).
#'
#' For predictor j with t-statistic t_j (from the *replicate's own*
#' refit) and residual df, the partial correlation is
#' \code{sign(t_j) * sqrt(t_j^2 / (t_j^2 + df))}; the semi-partial
#' correlation additionally rescales by \code{sqrt(1 - R^2)} of the
#' replicate's own fit. Both follow directly from the refit coefficients,
#' so (as with \code{.boot_lm_coef}) each replicate only needs one
#' Cholesky solve -- no separate model comparisons or refits with
#' predictors dropped are needed.
#'
#' @param X,y (weighted, offset-adjusted) design matrix / response, as
#'   passed to \code{.boot_lm_coef}
#' @param R number of bootstrap replicates
#' @param stat_kind \code{"pcor"} (partial correlation basis, for
#'   \code{partial.cor}/\code{partial.eta2}) or \code{"eta2"}
#'   (semi-partial correlation basis, for plain \code{eta2})
#' @return (R x p) matrix; column j holds the signed statistic for
#'   predictor j on each replicate (NA row if that replicate's fit failed
#'   or was degenerate)
#' @keywords internal
#' @noRd
.boot_lm_effect <- function(X, y, R, stat_kind = c("pcor", "eta2")) {
  stat_kind <- match.arg(stat_kind)
  n <- nrow(X); p <- ncol(X)
  has_icpt <- identical(colnames(X)[1L], "(Intercept)")
  df_resid <- n - p
  idx <- .boot_idx(n, R)
  out <- matrix(NA_real_, nrow = R, ncol = p)
  if (df_resid <= 0) return(out)

  for (b in seq_len(R)) {
    rows <- idx[, b]
    Xb <- X[rows, , drop = FALSE]
    yb <- y[rows]
    XtX <- crossprod(Xb)
    Xty <- crossprod(Xb, yb)
    ch <- tryCatch(chol(XtX), error = function(e) NULL)
    if (is.null(ch)) next
    beta_b <- backsolve(ch, backsolve(ch, Xty, transpose = TRUE))
    fitted_b <- as.vector(Xb %*% beta_b)
    resid_b <- yb - fitted_b
    rss_b <- sum(resid_b^2)
    if (rss_b <= 0) next
    XtX_inv_b <- tryCatch(chol2inv(ch), error = function(e) NULL)
    if (is.null(XtX_inv_b)) next
    sigma2_b <- rss_b / df_resid
    se_b <- sqrt(sigma2_b * diag(XtX_inv_b))
    t_b <- as.vector(beta_b) / se_b
    pr_b <- sign(t_b) * sqrt(t_b^2 / (t_b^2 + df_resid))
    if (stat_kind == "eta2") {
      tss_b <- if (has_icpt) sum((yb - mean(yb))^2) else sum(yb^2)
      r2_b <- if (tss_b > 0) 1 - rss_b / tss_b else NA_real_
      out[b, ] <- pr_b * sqrt(pmax(1 - r2_b, 0))
    } else {
      out[b, ] <- pr_b
    }
  }
  out
}

#' Closed-form (no-refit) leave-one-out values of the same signed
#' effect-size statistic as \code{.boot_lm_effect}, used only for the BCa
#' acceleration constant. Reuses the exact leave-one-out beta/RSS/TSS
#' identities already used by \code{.jack_lm_coef} (hat-matrix diagonal,
#' PRESS-type leave-one-out RSS, and a closed-form leave-one-out TSS), so
#' no observation ever needs to be refit.
#'
#' @inheritParams .boot_lm_effect
#' @param beta_hat,XtX_inv full-data OLS coefficients and (X'X)^{-1}
#' @return (n x p) matrix of leave-one-out statistic values
#' @keywords internal
#' @noRd
.jack_lm_effect <- function(X, y, beta_hat, XtX_inv, stat_kind = c("pcor", "eta2")) {
  stat_kind <- match.arg(stat_kind)
  n <- nrow(X); p <- ncol(X)
  has_icpt <- identical(colnames(X)[1L], "(Intercept)")
  df_full <- n - p
  df_loo <- (n - 1L) - p
  if (df_loo <= 0) return(matrix(NA_real_, nrow = n, ncol = p))

  fitted <- as.vector(X %*% beta_hat)
  e <- y - fitted
  h <- rowSums((X %*% XtX_inv) * X)
  h <- pmin(h, 1 - 1e-10)
  V <- X %*% XtX_inv                         # n x p, row i = (XtX_inv %*% x_i)'
  delta <- (e / (1 - h)) * V                 # n x p
  beta_loo <- sweep(-delta, 2, beta_hat, "+") # n x p leave-one-out betas

  rss_full <- sum(e^2)
  rss_loo  <- rss_full - e^2 / (1 - h)        # length n (PRESS-type identity)
  sigma2_loo <- rss_loo / df_loo              # length n

  diagXtXinv <- diag(XtX_inv)                                  # length p
  diagLOO <- sweep((V^2) / (1 - h), 2, diagXtXinv, "+")         # n x p
  se_loo <- sqrt(sigma2_loo * diagLOO)                          # n x p (row-recycled)
  t_loo  <- beta_loo / se_loo                                   # n x p
  pr_loo <- sign(t_loo) * sqrt(t_loo^2 / (t_loo^2 + df_loo))    # n x p

  if (stat_kind == "eta2") {
    Syy <- sum(y^2); Sy <- sum(y)
    tss_loo <- if (has_icpt) {
      (Syy - y^2) - (Sy - y)^2 / (n - 1L)      # closed-form leave-one-out TSS
    } else {
      Syy - y^2
    }
    r2_loo <- pmin(pmax(1 - rss_loo / tss_loo, 0), 1)  # length n
    pr_loo * sqrt(1 - r2_loo)                          # n x p (row-recycled)
  } else {
    pr_loo
  }
}

#' Fast case-resampling bootstrap of a marginal risk ratio or risk
#' difference for a binary exposure, via g-computation ("standardization"):
#' each replicate resamples whole rows (outcome *and* covariates
#' together), refits the binomial GLM on that replicate with a lean IRLS
#' solve (mirroring \code{.boot_glm_coef}), and then, still within that
#' same replicate's resampled covariate distribution, predicts the
#' response probability with the exposure column set to 1 for every
#' row and again set to 0 for every row. The marginal (population-
#' averaged) risks are the means of those two predicted-probability
#' vectors, and the reported effect is their ratio (\code{"rr"}) or
#' difference (\code{"rd"}). Resampling the standardization population
#' together with the outcome model on every replicate (rather than
#' fixing it at the observed data) is what makes this the standard
#' nonparametric bootstrap for a g-computed effect.
#'
#' @param X,y,weights,offset,start,irls.maxit,irls.tol as in
#'   \code{.boot_glm_coef} (binomial family only)
#' @param expo_col integer column index of the 0/1 exposure indicator in X
#' @param effect \code{"rr"} or \code{"rd"}
#' @return length-R numeric vector of replicate risk ratios/differences
#'   (\code{NA} for replicates whose refit failed or whose control-arm
#'   risk was 0, for \code{"rr"})
#' @keywords internal
#' @noRd
.boot_glm_gcomp <- function(X, y, R, family, weights, offset, start,
                             irls.maxit, irls.tol, expo_col,
                             effect = c("rr", "rd")) {
  effect <- match.arg(effect)
  n <- nrow(X); p <- ncol(X)
  idx <- .boot_idx(n, R)
  out <- rep(NA_real_, R)

  w_all <- weights
  o_all <- offset
  beta0 <- as.numeric(start)
  variance <- family$variance
  linkinv  <- family$linkinv
  mu.eta   <- family$mu.eta

  for (b in seq_len(R)) {
    rows <- idx[, b]
    Xb <- X[rows, , drop = FALSE]
    yb <- y[rows]
    wb <- w_all[rows]
    ob <- o_all[rows]

    beta_b <- tryCatch({
      beta <- beta0
      eta <- as.vector(Xb %*% beta) + ob
      for (it in seq_len(irls.maxit)) {
        mu <- linkinv(eta)
        mu.eta.val <- mu.eta(eta)
        mu.eta.val[mu.eta.val == 0] <- .Machine$double.eps
        z <- (eta - ob) + (yb - mu) / mu.eta.val
        W <- (mu.eta.val^2 / variance(mu)) * wb
        XW <- Xb * W
        ch <- chol(crossprod(Xb, XW))
        beta_new <- backsolve(ch, backsolve(ch, crossprod(XW, z), transpose = TRUE))
        if (max(abs(beta_new - beta)) < irls.tol * (max(abs(beta)) + irls.tol)) {
          beta <- beta_new
          break
        }
        beta <- beta_new
        eta <- as.vector(Xb %*% beta) + ob
        if (!all(is.finite(eta))) stop("non-finite linear predictor")
      }
      beta
    }, error = function(e) NULL)

    if (is.null(beta_b)) next

    X1 <- Xb; X1[, expo_col] <- 1
    X0 <- Xb; X0[, expo_col] <- 0
    R1 <- mean(linkinv(as.vector(X1 %*% beta_b) + ob))
    R0 <- mean(linkinv(as.vector(X0 %*% beta_b) + ob))
    out[b] <- if (effect == "rr") {
      if (R0 <= 0) NA_real_ else R1 / R0
    } else {
      R1 - R0
    }
  }
  out
}

#' Closed-form (no-refit) leave-one-out values of the g-computed risk
#' ratio/difference, used only for the BCa acceleration constant. Reuses
#' the one-step Newton (infinitesimal-jackknife) leave-one-out betas
#' already computed by \code{.jack_glm_coef} and applies the same
#' g-computation formula as \code{.boot_glm_gcomp}, standardizing over
#' the full (observed) covariate distribution for every leave-one-out
#' replicate -- consistent with the fact that this jackknife is itself
#' already a fast one-step approximation rather than an exact refit.
#'
#' Vectorized as a single (chunked) matrix product rather than one
#' n x p matrix-vector multiply per observation: for a block of k
#' leave-one-out beta vectors, predictions for all n standardization
#' rows under exposure = 1 (or 0) are the n x k product of the
#' counterfactual design matrix and the transposed beta block, so the
#' whole leave-one-out risk vector is obtained from one matrix multiply
#' per chunk instead of a length-n R loop. Chunking bounds peak memory
#' for large n (each chunk holds two n x chunk matrices at a time).
#'
#' @param X full design matrix (n x p) with the exposure column as
#'   observed
#' @param beta_loo (n x p) leave-one-out coefficient matrix from
#'   \code{.jack_glm_coef}
#' @param expo_col integer column index of the exposure indicator
#' @param linkinv the family's linkinv function
#' @param effect \code{"rr"} or \code{"rd"}
#' @return length-n numeric vector of leave-one-out risk ratios/differences
#' @keywords internal
#' @noRd
.jack_glm_gcomp <- function(X, beta_loo, expo_col, linkinv,
                             effect = c("rr", "rd"), offset = NULL) {
  effect <- match.arg(effect)
  n <- nrow(beta_loo)
  off <- if (is.null(offset)) rep(0, nrow(X)) else offset
  X1 <- X; X1[, expo_col] <- 1
  X0 <- X; X0[, expo_col] <- 0
  R1 <- numeric(n); R0 <- numeric(n)
  chunk <- if (n > 2000L) 500L else n
  starts <- seq(1L, n, by = chunk)
  for (s in starts) {
    e <- min(s + chunk - 1L, n)
    Bt   <- t(beta_loo[s:e, , drop = FALSE])  # p x k
    Eta1 <- X1 %*% Bt + off                    # n x k (off recycled by row)
    Eta0 <- X0 %*% Bt + off
    R1[s:e] <- colMeans(linkinv(Eta1))
    R0[s:e] <- colMeans(linkinv(Eta0))
  }
  if (effect == "rr") ifelse(R0 <= 0, NA_real_, R1 / R0) else R1 - R0
}

#' Identify the simple, untransformed, non-interaction predictor terms
#' eligible for population attributable fraction (PAF) reporting, and
#' classify each as "categorical" (factor/character/logical -- the
#' reference/counterfactual value is "everyone at the reference level",
#' i.e. all of that term's design-matrix dummy columns set to 0 under the
#' usual treatment contrasts) or "continuous" (numeric -- the
#' counterfactual value is "everyone at the sample mean").
#'
#' Interaction terms and terms that don't match a raw column name in
#' \code{data} (e.g. \code{poly(x, 2)}, \code{log(x)}) are not supported
#' -- the "set to reference/mean for everyone" counterfactual has no
#' unambiguous meaning for a transformed or interaction term -- so this
#' errors out (rather than silently skipping them) if any are found,
#' listing the offending term(s).
#'
#' @param mt the model \code{terms} object
#' @param data the original (raw, untransformed) data actually used by the
#'   fit
#' @param X the fitted model matrix (for its \code{"assign"} attribute)
#' @return a list; each element is \code{list(name, type, cols)} where
#'   \code{cols} are the design-matrix column indices for that term
#' @keywords internal
#' @noRd
.paf_recipes <- function(mt, data, X) {
  term.labels <- attr(mt, "term.labels")
  if (length(term.labels) == 0L) {
    stop("boot.glm(): effect = \"PAF\" requires at least one predictor term.",
         call. = FALSE)
  }
  assign_vec <- attr(X, "assign")
  bad <- character(0)
  recipes <- list()
  for (i in seq_along(term.labels)) {
    term <- term.labels[i]
    if (grepl(":", term, fixed = TRUE) || !term %in% names(data)) {
      bad <- c(bad, term)
      next
    }
    cols <- which(assign_vec == i)
    if (length(cols) == 0L) {
      bad <- c(bad, term)
      next
    }
    v <- data[[term]]
    type <- if (is.factor(v) || is.character(v) || is.logical(v)) "categorical" else "continuous"
    recipes[[length(recipes) + 1L]] <- list(name = term, type = type, cols = cols)
  }
  if (length(bad) > 0L) {
    stop("boot.glm(): effect = \"PAF\" only supports simple, untransformed, ",
         "non-interaction predictor terms (each term must match a raw ",
         "column name in 'data', with no interactions or transformations ",
         "such as poly()/log()); please refit without: ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  recipes
}

#' Fast case-resampling bootstrap of the population attributable fraction
#' (PAF) for every eligible predictor term at once, via g-computation:
#' each replicate resamples whole rows (outcome and covariates together),
#' refits the binomial GLM with the same lean IRLS solver used elsewhere
#' (\code{\link{.boot_glm_coef}}), computes that replicate's own observed
#' ("true") prevalence \code{mean(yb)}, and then -- still within that
#' same replicate's resampled covariate distribution and refit
#' coefficients -- predicts the counterfactual prevalence for every
#' recipe (categorical terms: all dummy columns set to 0; continuous
#' terms: set to that replicate's own resampled column mean). The PAF for
#' each term is \code{(true - counterfactual) / true}. Resampling the
#' standardization population together with the outcome model on every
#' replicate is the same "cases bootstrap" already used for
#' \code{.boot_glm_gcomp}'s risk ratio/difference.
#'
#' @param recipes list from \code{.paf_recipes}
#' @return (R x k) matrix, k = length(recipes); NA row for replicates
#'   whose refit failed or whose resampled prevalence was 0
#' @keywords internal
#' @noRd
.boot_glm_paf <- function(X, y, R, family, weights, offset, start,
                           irls.maxit, irls.tol, recipes) {
  n <- nrow(X); k <- length(recipes)
  idx <- .boot_idx(n, R)
  out <- matrix(NA_real_, nrow = R, ncol = k)

  w_all <- weights
  o_all <- offset
  beta0 <- as.numeric(start)
  variance <- family$variance
  linkinv  <- family$linkinv
  mu.eta   <- family$mu.eta

  for (b in seq_len(R)) {
    rows <- idx[, b]
    Xb <- X[rows, , drop = FALSE]
    yb <- y[rows]
    wb <- w_all[rows]
    ob <- o_all[rows]

    beta_b <- tryCatch({
      beta <- beta0
      eta <- as.vector(Xb %*% beta) + ob
      for (it in seq_len(irls.maxit)) {
        mu <- linkinv(eta)
        mu.eta.val <- mu.eta(eta)
        mu.eta.val[mu.eta.val == 0] <- .Machine$double.eps
        z <- (eta - ob) + (yb - mu) / mu.eta.val
        W <- (mu.eta.val^2 / variance(mu)) * wb
        XW <- Xb * W
        ch <- chol(crossprod(Xb, XW))
        beta_new <- backsolve(ch, backsolve(ch, crossprod(XW, z), transpose = TRUE))
        if (max(abs(beta_new - beta)) < irls.tol * (max(abs(beta)) + irls.tol)) {
          beta <- beta_new
          break
        }
        beta <- beta_new
        eta <- as.vector(Xb %*% beta) + ob
        if (!all(is.finite(eta))) stop("non-finite linear predictor")
      }
      beta
    }, error = function(e) NULL)

    if (is.null(beta_b)) next

    p_true_b <- mean(yb)
    if (!(p_true_b > 0)) next

    for (kk in seq_len(k)) {
      rec <- recipes[[kk]]
      Xcf <- Xb
      if (rec$type == "categorical") {
        Xcf[, rec$cols] <- 0
      } else {
        Xcf[, rec$cols] <- mean(Xb[, rec$cols])
      }
      p_cf <- mean(linkinv(as.vector(Xcf %*% beta_b) + ob))
      out[b, kk] <- (p_true_b - p_cf) / p_true_b
    }
  }
  out
}

#' Closed-form (no-refit) leave-one-out *counterfactual prevalence* for
#' every PAF recipe, used only for the BCa acceleration constant. Reuses
#' the one-step Newton leave-one-out betas already computed by
#' \code{.jack_glm_coef} (see \code{.jack_glm_gcomp} for the same
#' technique applied to a single exposure). For continuous recipes the
#' reference value is fixed at the *full-data* column mean rather than a
#' recomputed leave-one-out mean -- a small approximation, consistent
#' with this jackknife already being a fast one-step approximation rather
#' than an exact refit, that keeps the computation fully vectorized
#' (chunked matrix products) instead of needing a different reference
#' value for every one of the n pseudo-replicates.
#'
#' Combining this with the (separately computed, closed-form) leave-one-
#' out *observed* prevalence gives the leave-one-out PAF used for the BCa
#' acceleration constant; see boot_glm.R.
#'
#' @param recipes list from \code{.paf_recipes}
#' @return (n x k) matrix of leave-one-out counterfactual prevalences
#' @keywords internal
#' @noRd
.jack_glm_cf_prevalence <- function(X, beta_loo, linkinv, recipes, offset = NULL) {
  n <- nrow(X); k <- length(recipes)
  off <- if (is.null(offset)) rep(0, n) else offset
  out <- matrix(NA_real_, nrow = n, ncol = k)
  chunk <- if (n > 2000L) 500L else n
  starts <- seq(1L, n, by = chunk)
  for (kk in seq_len(k)) {
    rec <- recipes[[kk]]
    Xcf <- X
    if (rec$type == "categorical") {
      Xcf[, rec$cols] <- 0
    } else {
      Xcf[, rec$cols] <- mean(X[, rec$cols])
    }
    p_cf <- numeric(n)
    for (s in starts) {
      e <- min(s + chunk - 1L, n)
      Bt  <- t(beta_loo[s:e, , drop = FALSE])
      Eta <- Xcf %*% Bt + off
      p_cf[s:e] <- colMeans(linkinv(Eta))
    }
    out[, kk] <- p_cf
  }
  out
}

#' Fast case-resampling bootstrap for GLM coefficients using a lean,
#' warm-started, coefficients-only IRLS solver -- NOT stats::glm.fit().
#'
#' stats::glm.fit() is comparatively slow to call R times because, on every
#' single call, it: re-validates inputs and family functions, tracks
#' deviance and performs step-halving checks every iteration, computes a
#' full (rank-revealing) QR decomposition rather than a simple Cholesky
#' solve, and builds a complete result list (weights, R, qr, effects,
#' rank, ...) even though only the coefficient vector is needed here. All
#' of that is skipped below:
#'  - each replicate is warm-started from the full-data MLE (`start`), so
#'    convergence typically takes only 2-4 Newton steps instead of the
#'    ~6-10 glm.fit needs from a cold start;
#'  - the p x p weighted normal equations are solved with a Cholesky +
#'    two backsolves (cheaper than a QR solve for the small, well
#'    conditioned systems typical of a design matrix bootstrap);
#'  - convergence is checked via the size of the coefficient update
#'    (cheap) rather than recomputing the deviance (which requires an
#'    extra pass over dev.resids) every iteration;
#'  - a single tryCatch wraps the whole per-replicate iteration (not one
#'    per iteration), and failed/non-converged replicates are simply left
#'    as NA rows, which downstream BCa code already filters out.
#'  - the common special case family = gaussian(link = "identity") is
#'    just weighted OLS, solved exactly in one step (no IRLS iterations
#'    at all).
#'
#' @keywords internal
#' @noRd
.boot_glm_coef <- function(X, y, R, family, weights = NULL, offset = NULL,
                            start = NULL, irls.maxit = 25L, irls.tol = 1e-8) {
  n <- nrow(X); p <- ncol(X)
  idx <- .boot_idx(n, R)
  out <- matrix(NA_real_, nrow = R, ncol = p)

  w_all <- if (is.null(weights)) rep(1, n) else weights
  o_all <- if (is.null(offset)) rep(0, n) else offset
  beta0 <- if (is.null(start)) rep(0, p) else as.numeric(start)

  if (identical(family$family, "gaussian") && identical(family$link, "identity")) {
    ## exact weighted-least-squares solve, no iteration needed
    for (b in seq_len(R)) {
      rows <- idx[, b]
      Xb <- X[rows, , drop = FALSE]
      yb <- y[rows] - o_all[rows]
      sw <- sqrt(w_all[rows])
      Xw <- Xb * sw
      yw <- yb * sw
      out[b, ] <- tryCatch({
        ch <- chol(crossprod(Xw))
        backsolve(ch, backsolve(ch, crossprod(Xw, yw), transpose = TRUE))
      }, error = function(e) rep(NA_real_, p))
    }
    return(out)
  }

  variance <- family$variance
  linkinv  <- family$linkinv
  mu.eta   <- family$mu.eta

  for (b in seq_len(R)) {
    rows <- idx[, b]
    Xb <- X[rows, , drop = FALSE]
    yb <- y[rows]
    wb <- w_all[rows]
    ob <- o_all[rows]

    out[b, ] <- tryCatch({
      beta <- beta0
      eta <- as.vector(Xb %*% beta) + ob
      for (it in seq_len(irls.maxit)) {
        mu <- linkinv(eta)
        mu.eta.val <- mu.eta(eta)
        mu.eta.val[mu.eta.val == 0] <- .Machine$double.eps
        z <- (eta - ob) + (yb - mu) / mu.eta.val
        W <- (mu.eta.val^2 / variance(mu)) * wb
        XW <- Xb * W                       # row i scaled by W[i]
        ch <- chol(crossprod(Xb, XW))      # = t(X) %*% (W * X), p x p
        beta_new <- backsolve(ch, backsolve(ch, crossprod(XW, z), transpose = TRUE))
        if (max(abs(beta_new - beta)) < irls.tol * (max(abs(beta)) + irls.tol)) {
          beta <- beta_new
          break
        }
        beta <- beta_new
        eta <- as.vector(Xb %*% beta) + ob
        if (!all(is.finite(eta))) stop("non-finite linear predictor")
      }
      beta
    }, error = function(e) rep(NA_real_, p))
  }
  out
}

#' Approximate leave-one-out GLM coefficients via one-step Newton update
#' (uses working weights from the fitted model; avoids n full refits).
#' @keywords internal
#' @noRd
.jack_glm_coef <- function(X, y, fit) {
  n <- nrow(X)
  beta_hat <- stats::coef(fit)
  w <- fit$weights                       # IRLS working weights
  XtWX_inv <- tryCatch(chol2inv(chol(crossprod(X, X * w))),
                        error = function(e) MASS_ginv_fallback(crossprod(X, X * w)))
  eta <- fit$linear.predictors
  mu <- fit$fitted.values
  fam <- fit$family
  mu.eta <- fam$mu.eta(eta)
  ## working residual on the linear-predictor scale
  z_resid <- (y - mu) / ifelse(mu.eta == 0, .Machine$double.eps, mu.eta)
  h <- rowSums((X %*% XtWX_inv) * X) * w
  h <- pmin(h, 1 - 1e-10)
  XtXinvXt <- X %*% XtWX_inv
  delta <- (w * z_resid / (1 - h)) * XtXinvXt
  sweep(-delta, 2, beta_hat, "+")
}

#' Fallback generalized inverse (base R only) used only if a crossproduct
#' is numerically singular for the Cholesky route.
#' @keywords internal
#' @noRd
MASS_ginv_fallback <- function(m) {
  s <- svd(m)
  tol <- max(dim(m)) * max(s$d) * .Machine$double.eps
  pos <- s$d > tol
  s$v[, pos, drop = FALSE] %*% (1 / s$d[pos] * t(s$u[, pos, drop = FALSE]))
}
