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
- **Base R only.** No dependencies beyond `stats` and `utils`, which
  ship with every R installation.
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
  CI-inversion p-value.
- **`boot.glm()` marginal risk ratio/difference.** For a binomial model,
  `boot.glm(..., effect = "rr" | "rd", exposure = "<predictor>")` reports
  the marginal risk ratio or risk difference for a given binary
  predictor via g-computation, with a BCa confidence interval and
  CI-inversion p-value.
- **`boot.ci()`** takes an object you already fit with `t.test()`,
  `cor.test()`, `lm()`, or `glm()` and substitutes in BCa/CI-inversion
  inference, using the exact same internal code as the matching
  `boot.*()` function above -- so `boot.ci(lm(...))` and `boot.lm(...)`
  agree exactly under the shared default seed.

## Installation

```r
# from GitHub
# install.packages("remotes")
remotes::install_github("yourusername/merya")

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

## logistic regression
gfit <- boot.glm(am ~ wt + hp, data = mtcars, family = binomial())
summary(gfit)

## ... or the marginal risk ratio / risk difference for a binary
## predictor, via g-computation (binomial family only)
d <- mtcars; d$vs <- factor(d$vs)
rr_fit <- boot.glm(am ~ vs + wt, data = d, family = binomial(),
                    effect = "rr", exposure = "vs")
summary(rr_fit)

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
   partial/semi-partial correlations for `boot.lm()`'s effect sizes,
   g-computed risk ratios/differences for `boot.glm()`).
3. Estimate the BCa bias-correction (`z0`) and acceleration (`a`)
   constants — the latter from closed-form or one-step (infinitesimal
   jackknife) leave-one-out formulas, avoiding full refits.
4. Build the BCa confidence interval from the bootstrap distribution.
   For `boot.lm()`'s `"partial.eta2"`/`"eta2"`, the interval is first
   built on the signed (non-squared) statistic and then converted to
   the squared scale by taking the smallest/largest absolute value
   attained within that signed interval as the new lower/upper bound
   before squaring — not by naively squaring the signed endpoints,
   which would misrepresent an interval that straddles zero.
5. Obtain the p-value by **analytically inverting** the BCa
   transformation at the null value (0 for most statistics; 1 for a
   risk ratio), rather than searching over confidence levels.

Every function calls `set.seed()` once at the start (using the `seed`
argument, default `123`) before any resampling, so results are
reproducible across repeated calls by default.

See `?merya` for details.

## License

MIT © Your Name
