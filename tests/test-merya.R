## Simple base-R regression tests (no external test framework). The
## package has zero dependencies beyond stats/utils, so every test here
## runs unconditionally (no Suggests-gated skips). Run automatically by
## R CMD check as part of tests/.

set.seed(1)

## ---- boot.t.test ----------------------------------------------------
x <- rnorm(40, mean = 2)
y <- rnorm(40, mean = 0)

r1 <- merya::boot.t.test(x, R = 300)
stopifnot(inherits(r1, "htest"))
stopifnot(is.finite(r1$p.value), r1$p.value >= 0, r1$p.value <= 1)
stopifnot(length(r1$conf.int) == 2, r1$conf.int[1] < r1$conf.int[2])
## true mean is clearly away from 0 -> should reject at alpha = .05
stopifnot(r1$p.value < 0.05)

r2 <- merya::boot.t.test(x, y, R = 300)
stopifnot(inherits(r2, "htest"))
stopifnot(r2$p.value < 0.05)
## "Welch" was renamed to "Two independent samples"
stopifnot(grepl("Two independent samples", r2$method))
stopifnot(!grepl("Welch", r2$method))

r3 <- merya::boot.t.test(x, x + rnorm(40, sd = 0.01), paired = TRUE, R = 300)
stopifnot(inherits(r3, "htest"))
stopifnot(grepl("Paired", r3$method))

## one-sided alternatives produce a one-infinite-bound interval
r4 <- merya::boot.t.test(x, alternative = "greater", R = 300)
stopifnot(is.infinite(r4$conf.int[2]))
r5 <- merya::boot.t.test(x, alternative = "less", R = 300)
stopifnot(is.infinite(r5$conf.int[1]))

## ---- seed: default of 123 makes calls reproducible; NULL / explicit
## values change the outcome -----------------------------------------
r_a <- merya::boot.t.test(x, R = 300)          # default seed = 123
r_b <- merya::boot.t.test(x, R = 300)          # default seed = 123 again
stopifnot(isTRUE(all.equal(r_a$conf.int, r_b$conf.int)))
stopifnot(isTRUE(all.equal(r_a$p.value, r_b$p.value)))

r_c <- merya::boot.t.test(x, R = 300, seed = 999)
## a different seed should (with overwhelming probability) give a
## different bootstrap CI even though the data are identical
stopifnot(!isTRUE(all.equal(r_a$conf.int, r_c$conf.int)))

set.seed(777)
r_d1 <- merya::boot.t.test(x, R = 300, seed = NULL)  # uses current RNG stream
set.seed(777)
r_d2 <- merya::boot.t.test(x, R = 300, seed = NULL)  # same stream reproduced manually
stopifnot(isTRUE(all.equal(r_d1$conf.int, r_d2$conf.int)))

cat("boot.t.test / seed tests passed.\n")

## ---- boot.cor.test ---------------------------------------------------
n <- 60
a <- rnorm(n)
b <- 0.6 * a + rnorm(n, sd = 0.5)

rc1 <- merya::boot.cor.test(a, b, R = 300)
stopifnot(inherits(rc1, "htest"))
stopifnot(rc1$p.value < 0.05)
stopifnot(rc1$estimate > 0)

rc2 <- merya::boot.cor.test(a, b, method = "spearman", R = 300)
stopifnot(inherits(rc2, "htest"))

rc3 <- merya::boot.cor.test(a, b, method = "kendall", R = 100)
stopifnot(inherits(rc3, "htest"))

## independent variables -> should typically not reject
set.seed(2)
u <- rnorm(200); v <- rnorm(200)
rc4 <- merya::boot.cor.test(u, v, R = 300)
stopifnot(rc4$p.value > 0.01)

## default seed reproducibility
rc5a <- merya::boot.cor.test(a, b, R = 300)
rc5b <- merya::boot.cor.test(a, b, R = 300)
stopifnot(isTRUE(all.equal(rc5a$conf.int, rc5b$conf.int)))

## ---- boot.lm ----------------------------------------------------------
fit <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 300)
stopifnot(inherits(fit, "lm"))
stopifnot(inherits(fit, "boot.lm"))
s <- summary(fit)
stopifnot(inherits(s, "summary.boot.lm"))
## "Boot SE" must no longer appear in the summary table
stopifnot(all(c("Estimate", "CI lower", "CI upper", "Pr(>|z|)") %in%
              colnames(s$coefficients)))
stopifnot(!("Boot SE" %in% colnames(s$coefficients)))
## predict() and other lm methods must keep working
p <- predict(fit, newdata = mtcars[1:3, ])
stopifnot(length(p) == 3)
stopifnot(all(fit$boot$conf.int[, "lower"] <= fit$boot$conf.int[, "upper"]))

## default seed reproducibility for boot.lm
fit_r1 <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 300)
fit_r2 <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 300)
stopifnot(isTRUE(all.equal(fit_r1$boot$conf.int, fit_r2$boot$conf.int)))
fit_r3 <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 300, seed = 999)
stopifnot(!isTRUE(all.equal(fit_r1$boot$conf.int, fit_r3$boot$conf.int)))

## ---- boot.lm: effect = "partial.cor" / "partial.eta2" / "eta2" ----------
fit_pcor <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, effect = "partial.cor")
stopifnot(inherits(fit_pcor, "boot.lm"))
stopifnot(identical(fit_pcor$boot$effect, "partial.cor"))
stopifnot(!("(Intercept)" %in% rownames(fit_pcor$boot$conf.int)))
stopifnot(all(c("wt", "hp") %in% rownames(fit_pcor$boot$conf.int)))
## partial correlations are signed and bounded in [-1, 1]
stopifnot(all(abs(fit_pcor$boot$estimate) <= 1))
sp <- summary(fit_pcor)
stopifnot(identical(sp$effect, "partial.cor"))
stopifnot(!("Boot SE" %in% colnames(sp$coefficients)))

fit_peta2 <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, effect = "partial.eta2")
stopifnot(identical(fit_peta2$boot$effect, "partial.eta2"))
## eta2-type quantities are non-negative, as are their CI bounds
stopifnot(all(fit_peta2$boot$estimate >= 0))
stopifnot(all(fit_peta2$boot$conf.int >= 0, na.rm = TRUE))
## partial eta2 must equal (partial correlation)^2 for the same predictors
## (both use the default seed = 123, so their underlying resamples match)
fit_pcor2 <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, effect = "partial.cor")
fit_peta2_2 <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, effect = "partial.eta2")
stopifnot(isTRUE(all.equal(unname(fit_pcor2$boot$estimate^2),
                            unname(fit_peta2_2$boot$estimate))))

fit_eta2 <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 500, effect = "eta2")
stopifnot(identical(fit_eta2$boot$effect, "eta2"))
stopifnot(all(fit_eta2$boot$estimate >= 0))
stopifnot(all(fit_eta2$boot$conf.int >= 0, na.rm = TRUE))
## eta2 (non-partial) should not exceed partial eta2 for the same predictor
stopifnot(all(fit_eta2$boot$estimate <= fit_peta2_2$boot$estimate + 1e-8))
se <- summary(fit_eta2)
stopifnot(identical(se$effect, "eta2"))
stopifnot(!("Boot SE" %in% colnames(se$coefficients)))

## effect = "coef" (explicit) must still match the pre-existing behaviour
fit_default_a <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 300)
fit_default_b <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 300, effect = "coef")
stopifnot(isTRUE(all.equal(fit_default_a$boot$conf.int, fit_default_b$boot$conf.int)))

cat("boot.lm effect-size tests passed.\n")

## ---- boot.lm: pred.r.squared ---------------------------------------------
fit_predR_off <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 300)
stopifnot(is.null(fit_predR_off$boot$pred.r.squared))  # off by default

fit_predR <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 500,
                             pred.r.squared = TRUE)
pr <- fit_predR$boot$pred.r.squared
stopifnot(!is.null(pr))
stopifnot(pr$estimate >= 0)  # it's a square, must be non-negative
stopifnot(all(pr$conf.int >= 0, na.rm = TRUE))
stopifnot(pr$conf.int[1, "lower"] <= pr$conf.int[1, "upper"])
stopifnot(is.finite(pr$p.value), pr$p.value >= 0, pr$p.value <= 1)
## should coexist with any effect= choice
fit_predR_eta2 <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 500,
                                  effect = "eta2", pred.r.squared = TRUE)
stopifnot(!is.null(fit_predR_eta2$boot$pred.r.squared))
stopifnot(identical(fit_predR_eta2$boot$effect, "eta2"))
## predicted R-squared should not exceed ordinary (in-sample) R-squared
stopifnot(pr$estimate <= summary.lm(fit_predR)$r.squared + 1e-8)
## summary() should print without erroring and expose it
s_predR <- summary(fit_predR)
stopifnot(!is.null(s_predR$pred.r.squared))

## default seed reproducibility for pred.r.squared
fit_predR_r1 <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 300,
                                pred.r.squared = TRUE)
fit_predR_r2 <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 300,
                                pred.r.squared = TRUE)
stopifnot(isTRUE(all.equal(fit_predR_r1$boot$pred.r.squared$conf.int,
                            fit_predR_r2$boot$pred.r.squared$conf.int)))

cat("boot.lm pred.r.squared tests passed.\n")

## ---- boot.glm ----------------------------------------------------------
gfit <- merya::boot.glm(am ~ wt + hp, data = mtcars, family = binomial(), R = 300)
stopifnot(inherits(gfit, "glm"))
stopifnot(inherits(gfit, "boot.glm"))
gs <- summary(gfit)
stopifnot(inherits(gs, "summary.boot.glm"))
stopifnot(all(c("Estimate", "CI lower", "CI upper", "Pr(>|z|)") %in%
              colnames(gs$coefficients)))
stopifnot(!("Boot SE" %in% colnames(gs$coefficients)))
gp <- predict(gfit, newdata = mtcars[1:3, ], type = "response")
stopifnot(length(gp) == 3)
stopifnot(all(gfit$boot$conf.int[, "lower"] <= gfit$boot$conf.int[, "upper"]))

## default seed reproducibility for boot.glm
gfit_r1 <- merya::boot.glm(am ~ wt + hp, data = mtcars, family = binomial(), R = 300)
gfit_r2 <- merya::boot.glm(am ~ wt + hp, data = mtcars, family = binomial(), R = 300)
stopifnot(isTRUE(all.equal(gfit_r1$boot$conf.int, gfit_r2$boot$conf.int)))

## ---- boot.glm: effect = "OR" ---------------------------------------------
or_fit <- merya::boot.glm(am ~ wt + hp, data = mtcars,
                           family = binomial(link = "logit"), R = 500, effect = "OR")
stopifnot(identical(or_fit$boot$effect, "OR"))
stopifnot(identical(or_fit$boot$method, "direct"))
stopifnot(all(or_fit$boot$conf.int > 0, na.rm = TRUE))  # ORs are non-negative
## OR must be exp(coef) and its CI exp(coef CI); p-value must be identical
## to plain "coef" mode (same underlying bootstrap/jackknife replicates)
coef_fit <- merya::boot.glm(am ~ wt + hp, data = mtcars,
                             family = binomial(link = "logit"), R = 500)
stopifnot(isTRUE(all.equal(unname(exp(stats::coef(coef_fit))),
                            unname(or_fit$boot$estimate[names(stats::coef(coef_fit))]))))
stopifnot(isTRUE(all.equal(exp(coef_fit$boot$conf.int), or_fit$boot$conf.int,
                            check.attributes = FALSE)))
stopifnot(isTRUE(all.equal(coef_fit$boot$p.value, or_fit$boot$p.value,
                            check.attributes = FALSE)))
ors <- summary(or_fit)
stopifnot(identical(ors$effect, "OR"))
stopifnot(!("Boot SE" %in% colnames(ors$coefficients)))

## OR requires family = binomial(link = "logit")
stopifnot(inherits(
  tryCatch(merya::boot.glm(am ~ wt, data = mtcars, family = binomial(link = "probit"),
                            effect = "OR", R = 50),
           error = function(e) e),
  "error"))
stopifnot(inherits(
  tryCatch(merya::boot.glm(mpg ~ wt, data = mtcars, family = gaussian(),
                            effect = "OR", R = 50),
           error = function(e) e),
  "error"))

cat("boot.glm OR tests passed.\n")

## ---- boot.glm: effect = "RR" (conditional, direct exp(coef)) ------------
rr_direct <- merya::boot.glm(am ~ wt + hp, data = mtcars,
                              family = binomial(link = "log"), R = 500, effect = "RR")
stopifnot(identical(rr_direct$boot$effect, "RR"))
stopifnot(identical(rr_direct$boot$method, "direct"))
stopifnot(all(rr_direct$boot$conf.int > 0, na.rm = TRUE))
## must match exp(coef) from the same log-link model fit with effect="coef"
coef_log_fit <- merya::boot.glm(am ~ wt + hp, data = mtcars,
                                 family = binomial(link = "log"), R = 500)
stopifnot(isTRUE(all.equal(exp(coef_log_fit$boot$conf.int), rr_direct$boot$conf.int,
                            check.attributes = FALSE)))
stopifnot(isTRUE(all.equal(coef_log_fit$boot$p.value, rr_direct$boot$p.value,
                            check.attributes = FALSE)))

## conditional RR (no exposure) requires a log link, but works for ANY family
rr_poisson <- merya::boot.glm(carb ~ wt, data = mtcars,
                               family = poisson(link = "log"), R = 300, effect = "RR")
stopifnot(identical(rr_poisson$boot$effect, "RR"))
stopifnot(identical(rr_poisson$boot$method, "direct"))

## conditional RR without 'exposure' must be rejected for a non-log link
stopifnot(inherits(
  tryCatch(merya::boot.glm(am ~ wt, data = mtcars, family = binomial(link = "logit"),
                            effect = "RR", R = 50),
           error = function(e) e),
  "error"))

cat("boot.glm RR (direct) tests passed.\n")

## ---- boot.glm: effect = "RR" / "RD" (marginal, g-computation) -----------
d_rr <- mtcars
d_rr$vs <- factor(d_rr$vs)
rr_fit <- merya::boot.glm(am ~ vs + wt, data = d_rr, family = binomial(),
                           R = 500, effect = "RR", exposure = "vs")
stopifnot(inherits(rr_fit, "boot.glm"))
stopifnot(identical(rr_fit$boot$effect, "RR"))
stopifnot(identical(rr_fit$boot$method, "gcomputation"))
stopifnot(identical(rr_fit$boot$exposure, "vs"))
stopifnot(rr_fit$boot$estimate > 0)  # risk ratio is non-negative
stopifnot(rr_fit$boot$conf.int[1, "lower"] <= rr_fit$boot$conf.int[1, "upper"])
rrs <- summary(rr_fit)
stopifnot(identical(rrs$effect, "RR"))
stopifnot("vs" %in% rownames(rrs$coefficients))
stopifnot(!("Boot SE" %in% colnames(rrs$coefficients)))

## RR with 'exposure' given must run g-computation even under a log link
## (per spec: exposure given -> always marginal RR via g-computation)
rr_fit_log <- merya::boot.glm(am ~ vs + wt, data = d_rr, family = binomial(link = "log"),
                               R = 300, effect = "RR", exposure = "vs")
stopifnot(identical(rr_fit_log$boot$method, "gcomputation"))

rd_fit <- merya::boot.glm(am ~ vs + wt, data = d_rr, family = binomial(),
                           R = 500, effect = "RD", exposure = "vs")
stopifnot(identical(rd_fit$boot$effect, "RD"))
stopifnot(identical(rd_fit$boot$method, "gcomputation"))
stopifnot(rd_fit$boot$estimate >= -1 && rd_fit$boot$estimate <= 1)
stopifnot(rd_fit$boot$conf.int[1, "lower"] <= rd_fit$boot$conf.int[1, "upper"])

## marginal RR/RD (exposure given) must be rejected for non-binomial families
stopifnot(inherits(
  tryCatch(merya::boot.glm(mpg ~ wt, data = mtcars, family = gaussian(),
                            effect = "RR", exposure = "wt", R = 50),
           error = function(e) e),
  "error"))
stopifnot(inherits(
  tryCatch(merya::boot.glm(am ~ vs + wt, data = d_rr, family = binomial(),
                            effect = "RD", R = 50),
           error = function(e) e),
  "error"))

cat("boot.glm RR/RD (g-computation) tests passed.\n")

## ---- boot.glm: effect = "PAF" --------------------------------------------
paf_fit <- merya::boot.glm(am ~ vs + wt, data = d_rr, family = binomial(),
                            R = 500, effect = "PAF")
stopifnot(identical(paf_fit$boot$effect, "PAF"))
stopifnot(identical(paf_fit$boot$method, "gcomputation"))
## one row per eligible term: "vs" (categorical) and "wt" (continuous)
stopifnot(setequal(rownames(paf_fit$boot$conf.int), c("vs", "wt")))
pafs <- summary(paf_fit)
stopifnot(identical(pafs$effect, "PAF"))
stopifnot(!("Boot SE" %in% colnames(pafs$coefficients)))

## PAF requires family = binomial()
stopifnot(inherits(
  tryCatch(merya::boot.glm(mpg ~ wt, data = mtcars, family = gaussian(),
                            effect = "PAF", R = 50),
           error = function(e) e),
  "error"))

## PAF must reject interaction / transformed terms
stopifnot(inherits(
  tryCatch(merya::boot.glm(am ~ vs * wt, data = d_rr, family = binomial(),
                            effect = "PAF", R = 50),
           error = function(e) e),
  "error"))
stopifnot(inherits(
  tryCatch(merya::boot.glm(am ~ vs + log(wt), data = d_rr, family = binomial(),
                            effect = "PAF", R = 50),
           error = function(e) e),
  "error"))

cat("boot.glm PAF tests passed.\n")

cat("All merya smoke tests passed.\n")

## ---- boot.ci: generic dispatch + equivalence with boot.xxx() ------------

## htest: t.test (one-sample)
xt <- rnorm(30, mean = 1.5)
r_direct <- merya::boot.t.test(xt, R = 300)
tt <- t.test(xt)
r_via_ci <- merya::boot.ci(tt, R = 300)
stopifnot(inherits(r_via_ci, "htest"))
stopifnot(is.finite(r_via_ci$p.value), r_via_ci$p.value >= 0, r_via_ci$p.value <= 1)
stopifnot(r_via_ci$conf.int[1] < r_via_ci$conf.int[2])
stopifnot(r_via_ci$p.value < 0.05)  # true mean is clearly away from 0
## both used the default seed = 123, so results agree exactly
stopifnot(isTRUE(all.equal(r_direct$conf.int, r_via_ci$conf.int)))

## htest: t.test (two-sample), with x/y passed explicitly to sidestep
## data.name recovery entirely
xt2 <- rnorm(25, mean = 2); yt2 <- rnorm(25, mean = 0)
tt2 <- t.test(xt2, yt2)
r_via_ci2 <- merya::boot.ci(tt2, x = xt2, y = yt2, R = 300)
stopifnot(inherits(r_via_ci2, "htest"))
stopifnot(r_via_ci2$p.value < 0.05)

## htest: t.test with data.name auto-recovery (plain variable names)
xa <- rnorm(30, mean = 1.2); yb <- rnorm(30, mean = 0)
tt3 <- t.test(xa, yb)
r_via_ci3 <- merya::boot.ci(tt3, R = 300)
stopifnot(inherits(r_via_ci3, "htest"))
stopifnot(r_via_ci3$p.value < 0.05)

## htest: cor.test
ac <- rnorm(50); bc <- 0.6 * ac + rnorm(50, sd = 0.5)
ct <- cor.test(ac, bc)
r_via_ci4 <- merya::boot.ci(ct, R = 300)
stopifnot(inherits(r_via_ci4, "htest"))
stopifnot(r_via_ci4$p.value < 0.05)
stopifnot(r_via_ci4$estimate > 0)

ct_sp <- cor.test(ac, bc, method = "spearman")
r_via_ci5 <- merya::boot.ci(ct_sp, R = 200)
stopifnot(inherits(r_via_ci5, "htest"))

## lm: boot.ci(lm(...)) should exactly match boot.lm(...) under the same
## default seed (both default to seed = 123)
fit_lm <- lm(mpg ~ wt + hp, data = mtcars)
lm_direct <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 400)
lm_via_ci <- merya::boot.ci(fit_lm, R = 400)
stopifnot(inherits(lm_via_ci, "lm"))
stopifnot(inherits(lm_via_ci, "boot.lm"))
stopifnot(isTRUE(all.equal(lm_direct$boot$conf.int, lm_via_ci$boot$conf.int)))
stopifnot(isTRUE(all.equal(lm_direct$boot$p.value, lm_via_ci$boot$p.value)))
stopifnot(isTRUE(all.equal(unname(stats::coef(lm_direct)), unname(stats::coef(lm_via_ci)))))
## an explicit, matching, non-default seed must also agree
lm_direct_s <- merya::boot.lm(mpg ~ wt + hp, data = mtcars, R = 400, seed = 55)
lm_via_ci_s <- merya::boot.ci(fit_lm, R = 400, seed = 55)
stopifnot(isTRUE(all.equal(lm_direct_s$boot$conf.int, lm_via_ci_s$boot$conf.int)))

## glm: same equivalence check
fit_glm <- glm(am ~ wt + hp, data = mtcars, family = binomial())
glm_direct <- merya::boot.glm(am ~ wt + hp, data = mtcars, family = binomial(), R = 400)
glm_via_ci <- merya::boot.ci(fit_glm, R = 400)
stopifnot(inherits(glm_via_ci, "glm"))
stopifnot(inherits(glm_via_ci, "boot.glm"))
stopifnot(isTRUE(all.equal(glm_direct$boot$conf.int, glm_via_ci$boot$conf.int)))
stopifnot(isTRUE(all.equal(glm_direct$boot$p.value, glm_via_ci$boot$p.value)))

## boot.ci() on an unsupported class should error informatively
stopifnot(inherits(tryCatch(merya::boot.ci(list(a = 1)), error = function(e) e),
                    "error"))

cat("boot.ci (htest/lm/glm) smoke tests passed.\n")

cat("All boot.ci smoke tests passed.\n")
