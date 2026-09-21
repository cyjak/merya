# merya

<!-- badges: start -->
<!-- badges: end -->

**merya** provides drop-in replacements for `t.test()`, `cor.test()`,
`lm()`, and `glm()` that report **bias-corrected and accelerated (BCa)
bootstrap confidence intervals** and **p-values obtained by inverting
that confidence interval**, instead of relying on classical asymptotic
theory.

- Same arguments as the functions you already know.
- Same output classes (`"htest"`, `"lm"`, `"glm"`) — existing code,
  `print()`, `summary()`, `predict()`, etc. keep working.
- **Base R only.** No dependencies beyond `stats`, which ships with
  every R installation.
- **Wild bootstrap by default.** `boot.lm()`, `boot.glm()`, and
  `boot.t.test()` default to `boot.method = "wild"`: instead of
  resampling rows, the design matrix is kept fixed and each
  observation's own residual is perturbed by an independent random
  multiplier. This is the standard choice for heteroskedasticity-robust
  inference, and because the design (or, for `boot.glm()`, the weighted
  design) never changes across replicates, it needs **no explicit loop
  over replicates at all** — the whole bootstrap reduces to a handful
  of full-matrix operations, which is what makes it both the default
  and usually the faster option. The classical
  resample-the-rows-with-replacement bootstrap remains available via
  `boot.method = "case"`. The multiplier distribution (`wild.dist`) can
  be `"rademacher"` (default), `"mammen"`, or `"normal"`.
- **BCa or percentile intervals.** Every bootstrap-based function
  accepts `ci.type = "bca"` (default) or `"percentile"`, applied
  uniformly regardless of `boot.method`.
- **`boot.t.test()` matches `boot.lm()` exactly.** Internally,
  `boot.t.test()` is reformulated as a regression (intercept-only for
  one-sample/paired tests, on a 0/1 group indicator for two-sample
  tests) and reuses `boot.lm()`'s own engines directly — so
  `boot.t.test(x, y)` and `boot.lm(c(x, y) ~ group)`'s group
  coefficient agree exactly, not just approximately, for the same
  `boot.method`/`wild.dist`/`ci.type`/`seed`.
- **`boot.lm()` effect sizes.** Besides slope coefficients (the
  default), `boot.lm(..., effect = "partial.cor" | "partial.eta2" |
  "eta2")` reports, per predictor, its partial correlation, partial
  eta-squared, or eta-squared, each with a confidence interval and
  CI-inversion p-value. `boot.lm(..., pred.r.squared = TRUE)`
  additionally reports a bootstrapped, leave-one-out predicted
  R-squared for the model as a whole.
- **`boot.glm()` odds ratios, risk ratios/differences, and PAF.**
  `boot.glm(..., effect = "OR")` reports odds ratios (requires a logit
  link); `effect = "RR"` reports either the *conditional* risk ratio
  for every predictor (any family with a log link, if `exposure` is
  omitted) or the *marginal* risk ratio for one binary predictor via
  g-computation (binomial family, if `exposure` is given);
  `effect = "RD"` reports the marginal risk difference the same way;
  and `effect = "PAF"` reports the population attributable fraction for
  every eligible predictor (binomial family).
- **`boot.ci()`** takes an object you already fit with `t.test()`,
  `cor.test()`, `lm()`, or `glm()` and substitutes in bootstrap
  inference, using the exact same internal code as the matching
  `boot.*()` function above -- so `boot.ci(lm(...))` and `boot.lm(...)`
  agree exactly under the shared default settings.

## Installation

```r
# from GitHub
# install.packages("remotes")
remotes::install_github("cyjak/merya")

# from CRAN (once published)
install.packages("merya")
```

## Usage

```r
library(merya)

## t-test
x <- rnorm(30, mean = 1)
y <- rnorm(30, mean = 0)
boot.t.test(x, y)

## correlation
boot.cor.test(mtcars$wt, mtcars$mpg)

## linear regression (wild bootstrap by default)
fit <- boot.lm(mpg ~ wt + hp, data = mtcars)
summary(fit)

## ... or the classical case-resampling bootstrap
summary(boot.lm(mpg ~ wt + hp, data = mtcars, boot.method = "case"))

## ... or a plain percentile confidence interval instead of BCa
summary(boot.lm(mpg ~ wt + hp, data = mtcars, ci.type = "percentile"))

## ... and/or, for each predictor, an effect size instead of its slope
summary(boot.lm(mpg ~ wt + hp, data = mtcars, effect = "partial.cor"))
summary(boot.lm(mpg ~ wt + hp, data = mtcars, effect = "eta2"))

## ... and/or a leave-one-out predicted R-squared for the whole model
summary(boot.lm(mpg ~ wt + hp, data = mtcars, pred.r.squared = TRUE))

## logistic regression
gfit <- boot.glm(am ~ wt + hp, data = mtcars, family = binomial())
summary(gfit)

## odds ratios (requires a logit link)
summary(boot.glm(am ~ wt + hp, data = mtcars, family = binomial(),
                  effect = "OR"))

## conditional risk ratios for every predictor (any family, log link)
summary(boot.glm(carb ~ wt, data = mtcars, family = poisson(link = "log"),
                  effect = "RR"))

## marginal risk ratio / risk difference for a binary predictor,
## via g-computation (family = binomial(), any link)
d <- mtcars; d$vs <- factor(d$vs)
rr_fit <- boot.glm(am ~ vs + wt, data = d, family = binomial(),
                    effect = "RR", exposure = "vs")
summary(rr_fit)

## population attributable fraction for every eligible predictor
summary(boot.glm(am ~ vs + wt, data = d, family = binomial(), effect = "PAF"))

## consistency: boot.t.test(x, y) matches boot.lm() on the same data
## exactly (same seed, boot.method, ci.type)
d2 <- data.frame(z = c(x, y), group = c(rep(1, length(x)), rep(0, length(y))))
boot.t.test(x, y)
boot.lm(z ~ group, data = d2)

## reproducibility: same seed (default 123) -> same bootstrap result;
## a different seed, or seed = NULL, changes it
boot.lm(mpg ~ wt + hp, data = mtcars)                # seed = 123
boot.lm(mpg ~ wt + hp, data = mtcars)                # identical result
boot.lm(mpg ~ wt + hp, data = mtcars, seed = 456)    # different draw
boot.lm(mpg ~ wt + hp, data = mtcars, seed = NULL)   # uses current RNG state

## boot.ci(): apply bootstrap inference to an object you already fit
fit <- lm(mpg ~ wt + hp, data = mtcars)
summary(boot.ci(fit))          # same result boot.lm() would give
boot.ci(t.test(rnorm(30, 1)))  # works for t.test()/cor.test() output too
```

## How it works

`boot.t.test()`, `boot.cor.test()`, `boot.lm()`, and `boot.glm()`:

1. Draw bootstrap resamples (`boot.method = "case"`) or wild-bootstrap
   multipliers (`boot.method = "wild"`, the default) via a single
   vectorized call — for the wild bootstrap, this means a single
   `n x R` matrix of multipliers rather than an R-length loop, and
   because the (weighted) design matrix never changes across
   replicates, the entire bootstrap for every quantity this package
   supports reduces to a handful of full-matrix operations (one
   `crossprod()`, one `backsolve()` handling all `R` replicates at
   once) with no explicit loop over replicates at all.
2. Compute the statistic of interest for all replicates at once, using
   closed-form vectorized formulas where possible (means, mean
   differences, Pearson/Spearman correlation, OLS/GLM coefficients,
   partial/semi-partial correlations and leave-one-out predictions for
   `boot.lm()`'s effect sizes and `pred.r.squared`, exponentiated
   coefficients or g-computed risk ratios/differences/PAF for
   `boot.glm()`).
3. Estimate the BCa bias-correction (`z0`) and acceleration (`a`)
   constants (only used when `ci.type = "bca"`) — from closed-form or
   one-step (infinitesimal jackknife) leave-one-out formulas, avoiding
   full refits; this jackknife is shared between `boot.method`s, since
   its role (estimating the curvature of the sampling distribution from
   each observation's influence) doesn't depend on how the *main*
   bootstrap distribution was generated.
4. Build the confidence interval from the bootstrap distribution —
   the BCa interval (default), or the plain empirical
   `[alpha/2, 1-alpha/2]` quantiles if `ci.type = "percentile"`. For
   quantities that are *squares* of a signed statistic (`boot.lm()`'s
   `"partial.eta2"`/`"eta2"`), the interval is first built on the
   signed (non-squared) statistic and then converted to the squared
   scale by taking the smallest/largest absolute value attained within
   that signed interval as the new lower/upper bound before squaring —
   not by naively squaring the signed endpoints, which would
   misrepresent an interval that straddles zero. `pred.r.squared` is
   *not* treated this way: it's the textbook `1 - PRESS/TSS`, already on
   its natural signed scale, so its interval is built directly on that
   statistic with no transform. For `boot.glm()`'s `"OR"` and the
   conditional form of `"RR"`, the *same* coefficient-scale interval
   used for `"coef"` is instead just exponentiated (valid because
   exponentiation is monotonic, so it commutes with taking quantiles).
5. Obtain the p-value by **analytically inverting** the interval's
   construction at the natural null value for the scale in question (0
   for most statistics; 1 for a risk/odds ratio), rather than searching
   over confidence levels. For `"OR"`/conditional `"RR"`, the p-value is
   the one already computed on the coefficient scale and is not
   recomputed.

`boot.t.test()` is internally just a call into `boot.lm()`'s own
engines (see above), applied to an equivalent regression (see the
Usage example above) -- so it inherits `boot.method`, `wild.dist`, and
`ci.type` with no separate implementation, and agrees with `boot.lm()`
exactly on the same data.

Every function calls `set.seed()` once at the start (using the `seed`
argument, default `123`) before any resampling, so results are
reproducible across repeated calls by default.

See `?merya` for details.

## License

MIT © Your Name
