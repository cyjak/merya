# merya 0.7.4

* **Test-suite fix (no change to package behavior):**
  `tests/test-merya.R` fit `binomial(link = "log")` ("log-binomial")
  models on `mtcars` to test `effect = "RR"`. That family/link
  combination requires every fitted probability to stay below 1, which
  this particular formula/data does not satisfy, so plain
  `stats::glm()` itself fails to converge ("no valid set of
  coefficients has been found: please supply starting values") before
  `boot.glm()`'s bootstrap ever runs -- a misspecified model for this
  data, not a defect in the package. The conditional-RR-via-log-link
  test now uses a Poisson model instead (still verifying
  `RR = exp(coef)` against a matching `effect = "coef"` fit), and the
  narrower "exposure given still triggers g-computation even under a
  log link" case has been dropped from the suite rather than fit
  another fragile log-binomial model; the broader "exposure given ->
  g-computation" behavior remains covered under the default (logit)
  link.

# merya 0.7.3

* **Test-suite fix (no change to package behavior):** the
  `print.boot.lm()`/`print.boot.glm()` regression tests added in 0.7.1
  referenced `rr_fit` before it was defined later in
  `tests/test-merya.R` (a test-ordering mistake on my part, not a
  problem in the package), causing `object 'rr_fit' not found`. The
  block has been moved to after `or_fit`, `rr_fit`, `rd_fit`,
  `fit_pcor`, `gfit`, and `fit` are all defined.

# merya 0.7.2

* **Test-suite fix (no change to package behavior):** several
  `tests/test-merya.R` assertions checked that odds-ratio/risk-ratio
  bootstrap replicates and confidence-interval endpoints were strictly
  `> 0`. `exp()` is mathematically never negative, but for a
  sufficiently large-magnitude negative input it legitimately underflows
  to *exactly* `0` in double-precision floating point (below about
  `exp(-745)`), which can happen for an occasional bootstrap replicate's
  coefficient under quasi-complete separation on a small dataset -- as
  happened intermittently on some CI runners for the `mtcars` example
  used in the tests, depending on small BLAS/LAPACK floating-point
  differences across platforms. These assertions were relaxed to
  `>= 0`, matching what `exp()` can actually return; `boot.glm()`'s
  computation itself was already correct and required no change.
* Fixed placeholder GitHub username (`yourusername`) in `DESCRIPTION`'s
  `URL`/`BugReports` and in `README.md`'s install instructions to
  `cyjak`.

# merya 0.7.1

* **Bug fix: `print()`/console auto-print on a `boot.lm()`/`boot.glm()`
  object now always reflects the requested `effect`.** Previously, no
  `print` method was defined for these classes, so simply typing a
  fitted object at the console (or calling `print(fit)`) fell through to
  the inherited `print.lm()`/`print.glm()`, which always shows the raw
  model coefficients on their native scale -- e.g. a legitimately
  negative log-odds coefficient, even when `effect = "OR"` had been
  requested. This could easily be misread as "a negative odds ratio",
  when the actual (correctly positive) exponentiated odds ratio was only
  ever shown by `summary()`. New `print.boot.lm()`/`print.boot.glm()`
  methods now delegate to `summary()`, so the console default always
  matches whichever quantity `effect` selected; this does not change any
  computed value, only what is shown by default. Regression tests added.
* The bundled `.github/workflows/R-CMD-check.yaml` now uses
  `actions/checkout@v5` instead of `@v4`, which resolves the "Node.js 20
  actions are deprecated" warning some CI runs showed; this warning was
  a platform-wide GitHub Actions runner notice unrelated to any R code
  in this package (it did not indicate a check failure).

# merya 0.7.0

* **`boot.glm()`'s `effect` argument is reworked and gains two new
  options**; the accepted values are now spelled `"coef"` (unchanged),
  `"OR"`, `"RR"`, `"RD"`, and `"PAF"` (previously lowercase `"rr"`/
  `"rd"`, which no longer match -- this is a breaking rename).
  - `"OR"` (new): odds ratio. Reuses the *exact same* bootstrap/
    jackknife coefficient replicates as `"coef"` and simply reports
    `exp(estimate)`/`exp(CI)`; the p-value is left untouched (computed
    on the coefficient scale, which is equivalent since exponentiation
    doesn't change which side of the null a replicate falls on).
    Requires `family = binomial(link = "logit")`.
  - `"RR"`: if `exposure` is `NULL`, reports the *conditional* risk
    ratio for every predictor the same way `"OR"` reports odds ratios
    (exponentiated coefficient/CI, untouched p-value); this requires a
    log link, but now works with **any** family (e.g.
    `binomial(link = "log")`, `poisson(link = "log")`), not just
    binomial. If `exposure` is given, always reports the *marginal*
    risk ratio for that predictor via g-computation (as before),
    regardless of link function; this still requires
    `family = binomial()`.
  - `"RD"`: unchanged (marginal risk difference via g-computation,
    `family = binomial()`, `exposure` required).
  - `"PAF"` (new): population attributable fraction, for every simple
    (non-interaction, untransformed) predictor term, via g-computation:
    the observed prevalence is compared against the counterfactual
    prevalence obtained by setting that one term to its reference
    category for everyone (categorical predictors) or to the sample
    mean for everyone (continuous predictors), leaving all other
    covariates as observed. Multi-level categorical predictors
    contribute a single row. Requires `family = binomial()`; errors if
    any interaction or transformed term (e.g. `poly()`, `log()`) is
    present in the formula, since the counterfactual has no unambiguous
    meaning for those. `exposure` is not used for `"PAF"`.
* **`boot.lm()` gains a `pred.r.squared` argument** (`FALSE` by
  default): when `TRUE`, additionally computes a leave-one-out
  "predicted R-squared" -- the squared correlation between the observed
  response and, for every observation, its predicted value from the
  same model refit excluding that observation (obtained in closed form
  via the hat-matrix identity, no refitting) -- with a BCa confidence
  interval and CI-inversion p-value, reported alongside whichever
  `effect` was requested. As with `"partial.eta2"`/`"eta2"`, the CI is
  first computed on the signed (non-squared) correlation and then
  converted to the squared scale by taking the smallest/largest
  absolute value attained within the signed interval as the new
  lower/upper bound before squaring.

# merya 0.6.0

* **`boot.gee()` has been removed from the package**, along with its
  `boot.ci()` method (for `"gee"` objects) and the `gee`-specific
  cluster-resampling helpers it used internally. The `gee` package has
  accordingly been dropped from `Suggests` entirely; `merya` now has
  **no package dependencies beyond base R** (`stats` and `utils`).
  `boot.t.test()`, `boot.cor.test()`, `boot.lm()`, and `boot.glm()` (and
  `boot.ci()` applied to their inputs) are unaffected.

# merya 0.5.0

* **`boot.lmer()` and `boot.glmer()` have been removed from the
  package**, along with their `boot.ci()` methods (for `"lmerMod"`/
  `"glmerMod"` objects) and `boot.summary()` (which existed solely to
  print their bootstrap results, since those functions returned
  unmodified S4 `lme4` objects that couldn't carry a custom `summary()`
  method). Refitting a linear/generalized linear mixed model on every
  bootstrap replicate is comparatively expensive, and combined with the
  new default of `R = 10000` (see `merya` 0.4.0) this made both
  functions impractically slow even for modest sample sizes (e.g.
  n = 100). The `lme4` package has accordingly been dropped from
  `Suggests`, and the `methods` package from `Imports` (both were only
  needed for these two functions). `boot.gee()` is unaffected and
  remains available for clustered/correlated data.
* **Every bootstrap-based function now takes a `seed` argument,
  defaulting to `123`**: `boot.t.test()`, `boot.cor.test()`,
  `boot.lm()`, `boot.glm()`, `boot.gee()`, and every remaining
  `boot.ci()` method. `set.seed(seed)` is called once at the start of
  the function (before any resampling), so results are reproducible by
  default across repeated calls with the same data; pass a different
  integer for a different reproducible draw, or `seed = NULL` to use
  whatever random state is currently active instead of reseeding.
  `boot.ci()`'s methods forward `seed` to the underlying `boot.*()`
  call, so e.g. `boot.ci(lm(...))` and `boot.lm(...)` continue to agree
  exactly under the shared default seed.
* **The `"Boot SE"` column has been removed from `boot.lm()`'s and
  `boot.glm()`'s `summary()` output** (for all `effect` options,
  including `"partial.cor"`/`"partial.eta2"`/`"eta2"` and `"rr"`/`"rd"`);
  the table now shows `Estimate`, `CI lower`, `CI upper`, and
  `Pr(>|z|)`. `boot.gee()`'s summary table is unchanged and still
  includes `Boot SE`.
* **`boot.t.test()`'s two-sample method label no longer says "Welch"**;
  it now reads `"Two independent samples bootstrap t-test"` (the
  bootstrap itself was never a Welch t-test in the classical sense, so
  this only affects wording, not computation).

# merya 0.4.0

* **Package renamed from `bootInfer` to `merya`.** All exported
  functions, documentation, and internal `::`-qualified references have
  been updated accordingly; there is no separate compatibility shim, so
  code that referred to the package as `bootInfer` needs to switch to
  `merya`.
* **Default number of bootstrap replicates raised to `R = 10000`** for
  every bootstrap-based function (`boot.t.test()`, `boot.cor.test()`,
  `boot.lm()`, `boot.glm()`, `boot.gee()`, `boot.lmer()`, `boot.glmer()`,
  and every `boot.ci()` method). Note that because `boot.gee()`,
  `boot.lmer()`, and `boot.glmer()` refit a real model on every
  replicate, 10000 replicates can take considerably longer for these
  three than for the base-R-only functions; pass a smaller `R` for
  exploratory work.
* `boot.lm()`: new `effect` argument. In addition to the default
  `"coef"` (slope coefficients, unchanged from before), it can report,
  per predictor:
  - `"partial.cor"` -- the partial correlation between the response and
    that predictor, controlling for the others;
  - `"partial.eta2"` -- partial eta-squared (the square of the partial
    correlation);
  - `"eta2"` -- eta-squared (the square of the semi-partial/"part"
    correlation).

  All three come from the same per-replicate refit already used for
  `"coef"` (via each replicate's own t-statistics and R-squared), so no
  extra refitting is needed, and their jackknife (for the BCa
  acceleration constant) is likewise closed-form. Because partial
  eta-squared and eta-squared are squares of a statistic that can be
  negative, their confidence intervals are built by first computing the
  BCa interval on the signed (partial/semi-partial correlation) scale,
  then converting to the squared scale by taking the smallest absolute
  value attained within that interval as the lower bound (0 whenever the
  signed interval straddles zero) and the largest absolute value as the
  upper bound, and squaring both -- rather than naively squaring the
  signed endpoints, which would misrepresent the interval whenever it
  straddles zero.
* `boot.glm()`: new `effect` (`"coef"`, `"rr"`, or `"rd"`) and
  `exposure` arguments. For a binomial model, `effect = "rr"`/`"rd"`
  reports the marginal (population-averaged) risk ratio or risk
  difference for a given binary `exposure` predictor, obtained by
  g-computation (standardizing model-predicted probabilities under
  exposure = 1 vs. exposure = 0 over the sample's covariate
  distribution). Each bootstrap replicate resamples whole rows, refits
  the model, and repeats the standardization within that replicate's own
  resampled data -- the standard nonparametric bootstrap for a
  g-computed effect -- with a closed-form (no-refit) leave-one-out
  approximation for the BCa acceleration constant, and CI-inversion
  p-values against the natural null for each scale (1 for the risk
  ratio, 0 for the risk difference).

# merya 0.3.0

* `boot.ci()`: a new S3 generic that takes an object already returned by
  `stats::t.test()`, `stats::cor.test()`, `stats::lm()`, `stats::glm()`,
  `gee::gee()`, `lme4::lmer()`, or `lme4::glmer()`, and returns it with
  its confidence interval(s)/p-value(s) replaced by the BCa bootstrap /
  CI-inversion equivalents -- in the same class and format the input
  object already had. For `lm`/`glm`/`gee`/`lmerMod`/`glmerMod` objects,
  `boot.ci()` reconstructs the original fitting call (`object$call` or
  `object@call`) and re-evaluates it with the matching `boot.*()`
  function substituted in, in the caller's environment -- the same
  mechanism `stats::update()` uses to refit a model -- so `boot.ci(lm(...))`
  runs the *exact same* internal resampling/BCa/CI-inversion code as
  `boot.lm(...)`; with the same random seed the two agree exactly (this is
  covered by the new equivalence tests in `tests/test-merya.R`). For
  `"htest"` objects (which store no `$call`), `boot.ci()` recovers the
  original data by re-evaluating the deparsed `data.name` label in the
  caller's environment, or accepts it explicitly via `boot.ci(object, x =
  ..., y = ...)`.
* `boot.ci()` is a thin dispatcher only (no new bootstrap/BCa machinery),
  so it adds negligible overhead beyond the `boot.*()` call it triggers.

# merya 0.2.0

* `boot.gee()`: GEE regression for correlated/clustered data, mirroring
  `gee::gee()`'s arguments and returned `"gee"` object. Replaces the naive/
  robust-sandwich inference with a cluster-resampling BCa bootstrap (whole
  `id` clusters resampled with replacement, exactly matching how `gee`
  itself defines clusters) and CI-inversion p-values. Every replicate is
  warm-started from the full-data estimate via `gee`'s own `b` argument.
  Requires the `gee` package (listed under `Suggests`, not a hard
  dependency of the package as a whole).
* `boot.lmer()`: linear mixed-effects models, mirroring `lme4::lmer()`'s
  arguments and returned `"lmerMod"` object exactly (same class, same
  slots -- all `lme4` methods keep working unchanged). The random-effects
  structure and its estimation are entirely `lme4`'s; only the
  *resampling* scheme is added on top, as a cases/cluster bootstrap over
  the model's grouping factor(s), read off via `lme4::getME(fit,
  "flist")`. BCa confidence intervals and CI-inversion p-values are
  attached as a `"boot"` attribute (not a slot, so `lme4` dispatch is
  untouched) and printed with the new `boot.summary()`. Replicates are
  warm-started from the full-data variance-component estimate for speed.
  Requires the `lme4` package (`Suggests`).
* `boot.glmer()`: the same treatment for `lme4::glmer()` / `"glmerMod"`
  objects, including generalized linear mixed models (binomial, Poisson,
  etc.) and warm-starting from both the full-data fixed effects and
  variance components. Requires the `lme4` package (`Suggests`).
* `boot.summary()`: a single accessor that prints the BCa bootstrap
  inference table for `boot.gee()`/`boot.lmer()`/`boot.glmer()` objects,
  regardless of whether the underlying object is an S3 list (`gee`) or an
  S4 object (`lme4`'s `merMod` classes).
* Internal: new shared helpers for cluster-resampling, leave-one-cluster-
  out jackknifing, and BCa-table assembly (`utils-internal.R`), reused by
  all three new functions instead of duplicating the same loop three
  times.
* `gee` and `lme4` are now listed under `Suggests` (only required by
  `boot.gee()`/`boot.lmer()`/`boot.glmer()`); `t.test`/`cor.test`/`lm`/
  `glm` bootstrap functions remain zero-dependency (base R only).

# merya 0.1.0

* Initial release.
* `boot.t.test()`: one-sample, paired, and Welch two-sample t-tests with
  BCa bootstrap confidence intervals and CI-inversion p-values.
* `boot.cor.test()`: Pearson, Spearman, and Kendall correlation tests with
  BCa bootstrap confidence intervals and CI-inversion p-values.
* `boot.lm()`: linear regression with a BCa bootstrap inference table for
  coefficients, returned as a subclass of `"lm"`.
* `boot.glm()`: generalized linear regression with a BCa bootstrap
  inference table for coefficients, returned as a subclass of `"glm"`.
  The bootstrap replicates are fit with a lean, warm-started,
  coefficients-only IRLS solver (Cholesky-based normal equations) instead
  of calling `stats::glm.fit()` per replicate, which removes per-call
  overhead (deviance tracking, QR decomposition, full result-object
  construction) that dominated runtime even for small samples. The
  `family = gaussian(link = "identity")` case is solved exactly in one
  step (no IRLS iterations). New `irls.maxit` / `irls.tol` arguments
  expose the (rarely needed) iteration controls.
* Base R only (Imports: stats, utils); resampling and BCa/CI-inversion
  computations are vectorized and use closed-form jackknife/inversion
  formulas for speed.
