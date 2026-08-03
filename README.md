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
- **Built for speed.** Resampling is vectorized wherever a closed-form
  statistic allows it, and the BCa transform is inverted *analytically*
  for p-values.
- **All bootstrap-based functions default to `R = 10000` replicates and
  to a fixed `seed = 123`** for reproducibility. Pass a different
  integer for a different reproducible draw, or `seed = NULL` to use
  whatever random state is currently active instead of reseeding.
- **`boot.lm()` effect sizes.** Besides slope coefficients (the
  default), `boot.lm(..., effect = "partial.cor" | "partial.eta2" |
  "eta2")` reports, per predictor, its partial correlation, partial
  eta-squared, or eta-squared, each with a BCa confidence interval and
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
  every eligible predictor (binomial family). All come with a BCa
  confidence interval and CI-inversion p-value.
- **`boot.ci()`** takes an object you already fit with `t.test()`,
  `cor.test()`, `lm()`, or `glm()` and substitutes in BCa/CI-inversion
  inference, using the exact same internal code as the matching
  `boot.*()` function above -- so `boot.ci(lm(...))` and `boot.lm(...)`
  agree exactly under the shared default seed.

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

## linear regression
fit <- boot.lm(mpg ~ wt + hp, data = mtcars)
summary(fit)

## ... or, for each predictor, an effect size instead of its slope
summary(boot.lm(mpg ~ wt + hp, data = mtcars, effect = "partial.cor"))
summary(boot.lm(mpg ~ wt + hp, data = mtcars, effect = "partial.eta2"))
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

## reproducibility: same seed (default 123) -> same bootstrap result;
## a different seed, or seed = NULL, changes it
boot.lm(mpg ~ wt + hp, data = mtcars)                # seed = 123
boot.lm(mpg ~ wt + hp, data = mtcars)                # identical result
boot.lm(mpg ~ wt + hp, data = mtcars, seed = 456)    # different draw
boot.lm(mpg ~ wt + hp, data = mtcars, seed = NULL)   # uses current RNG state

## boot.ci(): apply BCa/CI-inversion inference to an object you already fit
fit <- lm(mpg ~ wt + hp, data = mtcars)
summary(boot.ci(fit))          # same result boot.lm() would give
boot.ci(t.test(rnorm(30, 1)))  # works for t.test()/cor.test() output too
```

## How it works

`boot.t.test()`, `boot.cor.test()`, `boot.lm()`, and `boot.glm()`:

1. Draw bootstrap resamples via a single `n x R` index matrix (not an
   R-length loop).
2. Compute the statistic of interest for all replicates at once, using
   closed-form vectorized formulas where possible (means, mean
   differences, Pearson/Spearman correlation, OLS/GLM coefficients,
   partial/semi-partial correlations and leave-one-out predictions for
   `boot.lm()`'s effect sizes and `pred.r.squared`, exponentiated
   coefficients or g-computed risk ratios/differences/PAF for
   `boot.glm()`).
3. Estimate the BCa bias-correction (`z0`) and acceleration (`a`)
   constants — the latter from closed-form or one-step (infinitesimal
   jackknife) leave-one-out formulas, avoiding full refits.
4. Build the BCa confidence interval from the bootstrap distribution.
   For quantities that are *squares* of a signed statistic
   (`boot.lm()`'s `"partial.eta2"`/`"eta2"`), the interval is first built
   on the signed (non-squared) statistic and then converted to the
   squared scale by taking the smallest/largest absolute value attained
   within that signed interval as the new lower/upper bound before
   squaring — not by naively squaring the signed endpoints, which would
   misrepresent an interval that straddles zero. `pred.r.squared` is
   *not* treated this way: it's the textbook `1 - PRESS/TSS`, already on
   its natural signed scale (a model worse than the mean legitimately
   scores negative), so its BCa interval is built directly on that
   statistic with no transform. For `boot.glm()`'s `"OR"` and the
   conditional form of `"RR"`, the *same* coefficient-scale BCa interval
   used for `"coef"` is instead just exponentiated (valid because
   exponentiation is monotonic, so it commutes with taking quantiles).
5. Obtain the p-value by **analytically inverting** the BCa
   transformation at the natural null value for the scale in question
   (0 for most statistics; 1 for a risk/odds ratio), rather than
   searching over confidence levels. For `"OR"`/conditional `"RR"`, the
   p-value is the one already computed on the coefficient scale
   (equivalent, since exponentiation doesn't change which side of the
   null a replicate falls on) and is not recomputed.

Every function calls `set.seed()` once at the start (using the `seed`
argument, default `123`) before any resampling, so results are
reproducible across repeated calls by default.

See `?merya` for details.

## License

MIT © Your Name
