# merya 0.9.0

* **Wild bootstrap, now the default resampling method for `boot.lm()`,
  `boot.glm()`, and `boot.t.test()`**, via a new `boot.method` argument
  (`"wild"` or `"case"`). The wild bootstrap keeps the design matrix
  fixed across every replicate and perturbs each observation's own
  residual by an independent random multiplier (mean 0, variance 1)
  instead of resampling rows -- the standard approach for
  heteroskedasticity-robust bootstrap inference. Because the design
  matrix (for `boot.glm()`, the IRLS-weighted design at convergence)
  never changes across replicates, the entire `R`-replicate bootstrap
  for every quantity these functions support reduces to a handful of
  full-matrix operations (one `crossprod()`, one `backsolve()` handling
  all replicates' right-hand sides at once) with **no explicit loop
  over replicates at all**, which is what makes `"wild"` both the new
  default and, in most cases, faster than the classical
  resample-the-rows-with-replacement bootstrap, which remains available
  as `boot.method = "case"` and is unchanged in its own methodology.
  - `boot.lm()`: every effect option (`"coef"`, `"partial.cor"`,
    `"partial.eta2"`, `"eta2"`) and `pred.r.squared` are supported under
    wild bootstrap, derived from one shared set of wild-bootstrap draws
    per call. The leave-one-out predictions needed for
    `pred.r.squared` reuse the closed-form OLS identity
    `yhat_(-i) = yhat_i - h_i*e_i/(1-h_i)` directly, with no
    per-replicate leave-one-out coefficient matrix needed.
  - `boot.glm()`: since a response-scale wild bootstrap cannot generally
    be defined for non-Gaussian families (e.g. a perturbed 0/1 response
    for `binomial()` would fall outside `{0,1}`), the wild bootstrap
    here perturbs the residual of the IRLS *working response* at
    convergence and solves **one** weighted-least-squares step per
    replicate using the converged IRLS weights -- a linearization
    around the full-data fit, well-defined for any family/link. All
    `effect` options (`"coef"`, `"OR"`, `"RR"`, `"RD"`, `"PAF"`) are
    supported; the g-computation effects (`"RR"`/`"RD"` with
    `exposure`, and `"PAF"`) standardize over the fixed, full covariate
    distribution using each replicate's own linearized coefficients,
    fully vectorized. For `"PAF"` specifically, the observed prevalence
    (which has no wild-bootstrap analogue, since only residuals are
    perturbed) is held fixed at the full-sample value for every
    replicate; only the model-based counterfactual prevalence varies.
  - The BCa acceleration constant's jackknife is unchanged and shared
    between both `boot.method`s in every function, since its role
    (estimating the curvature of the sampling distribution from each
    observation's influence) does not depend on how the main bootstrap
    distribution is generated.
  - `wild.dist` selects the multiplier distribution: `"rademacher"`
    (default; \eqn{v=\pm1} with equal probability, the simplest and
    most commonly recommended choice), `"mammen"` (a skewed two-point
    distribution matching the third moment), or `"normal"`
    (\eqn{v \sim N(0,1)}).
* **`ci.type` argument (`"bca"`, the previous and still-default
  behavior, or `"percentile"`) added to every bootstrap-based function**
  (`boot.t.test()`, `boot.cor.test()`, `boot.lm()`, `boot.glm()`),
  applying uniformly regardless of `boot.method`. `"percentile"` uses
  the plain empirical `[alpha/2, 1-alpha/2]` quantiles of the bootstrap
  distribution (no bias-correction or acceleration adjustment), with a
  correspondingly simpler CI-inversion p-value (the smallest two-sided
  alpha at which the null value would just sit on the interval's
  boundary).
* **`boot.t.test()` is reformulated internally as a regression**, and
  now agrees *exactly* (not just approximately) with `boot.lm()` on the
  same data: a one-sample test is an intercept-only regression on `x`
  (paired: on the differences `x - y`), and a two-sample test is a
  regression of the combined sample on a 0/1 group indicator, testing
  that indicator's coefficient (algebraically identical to the
  difference in means). Both cases are handed directly to `boot.lm()`'s
  own internal engines (the same bootstrap/wild-bootstrap coefficient
  routines and the same closed-form jackknife), so
  `boot.t.test(x, y)` and `boot.lm(c(x, y) ~ group)`'s group
  coefficient are numerically identical given the same
  `boot.method`/`wild.dist`/`ci.type`/`seed`. **This is a deliberate
  methodology change, not just an internal refactor**: the two-sample
  case's resampling is now unstratified (the *combined* sample is
  resampled/perturbed together, exactly as `boot.lm()`'s case-resampling
  bootstrap already did), rather than resampling each group separately
  at its own fixed size as previous versions did; this is what makes
  the exact `boot.lm()` agreement possible, and is also the only
  sensible option for the wild bootstrap, which has no natural
  "resample within group" analogue since it never resamples rows.
  `var.equal` remains an ignored, interface-compatibility-only argument,
  as before.

# merya 0.8.2

* **Test-suite fix (no change to package behavior):** the
  single-predictor `sr == pr` regression test added in 0.8.1 to verify
  the semi-partial correlation fix compared `fit$boot$estimate` between
  an `effect = "partial.cor"` fit and an `effect = "eta2"` fit directly.
  That's not an apples-to-apples comparison: `fit$boot$estimate` holds
  the *reported* quantity, which for `"eta2"` is already squared
  (`eta2 = sr^2`, mirroring how `"partial.eta2"` reports `pr^2`), while
  `"partial.cor"` reports the raw, unsquared statistic -- so the test
  was comparing `pr` against `sr^2` and failing even though the
  underlying fix is correct (`sr == pr` exactly, verified independently
  via `fit$boot$coefficients`, which holds each side's *raw*,
  pre-squaring per-replicate values and does not have this mismatch).
  The test now squares the `"partial.cor"` side before comparing,
  matching the same convention already used by the pre-existing
  `"partial.eta2" == "partial.cor"^2` test just above it in the suite.

# merya 0.8.1

* **Bug fix: `boot.lm(..., effect = "eta2")` (and, downstream, its BCa
  confidence interval and p-value) used an incorrect formula for the
  semi-partial ("part") correlation.** It computed `sr <- pr *
  sqrt(1 - R2)` (partial correlation rescaled by `sqrt(1 - R^2)`); the
  correct closed-form expression (Cohen & Cohen, *Applied Multiple
  Regression/Correlation Analysis*), derivable from the general-linear-
  model identity `F = t^2` for dropping a single predictor, is
  `sr <- t * sqrt((1 - R2) / df_resid)`. The two formulas coincide only
  in the trivial `t -> 0` limit and diverge increasingly as the effect
  gets stronger, so `eta2` (and its confidence interval) was previously
  understated for anything beyond a weak effect. The clearest symptom:
  with a single predictor, there is nothing left to partial out, so
  `sr` must equal `pr` exactly -- the old formula did not satisfy this
  (`pr * sqrt(1-R2) != pr` unless `R2 = 0`); the corrected one does, by
  construction. Fixed in all three places this quantity is computed:
  the per-replicate bootstrap loop, the closed-form leave-one-out
  jackknife used for the BCa acceleration constant, and the full-sample
  point estimate; `"partial.cor"`/`"partial.eta2"` (which do not
  involve this formula) are unaffected. Two correctness checks were
  added to the test suite: the single-predictor `sr == pr` identity
  above (which the bug failed and the fix satisfies, exactly, across
  the entire bootstrap distribution, not just the point estimate), and
  a direct ground-truth check with two correlated predictors comparing
  the reported semi-partial correlation squared against
  `R2_full - R2_reduced` computed by actually refitting the model
  without that predictor.

# merya 0.8.0

* **`boot.lm()`'s `pred.r.squared` now uses the textbook "predicted
  R-squared" definition, `1 - PRESS/TSS`, instead of a squared
  leave-one-out correlation.** This is a deliberate behavior change, not
  a bug fix: the previous implementation bootstrapped the (signed)
  correlation between the observed response and leave-one-out
  predictions and then squared it using the same sign-preserving
  transform as `"partial.eta2"`/`"eta2"`. That transform is only valid
  when the population-level quantity is known to be non-negative (true
  for eta-squared, a variance decomposition), which is *not* true of
  predictive performance -- a model that predicts worse than just the
  response mean has genuinely, meaningfully negative predictive value,
  and squaring silently converted that into a small positive number
  (e.g. a correlation of -0.20 with 95% CI [-0.5, 0.1] became a
  "predicted R-squared" of 0.04 with CI [0.01, 0.25], reading as weak
  evidence of usefulness when the evidence actually pointed the other
  way). `1 - PRESS/TSS` is already on its natural, correctly signed,
  unbounded-below scale, so it is now bootstrapped and given a BCa
  interval directly, with no post-hoc transform at all. `PRESS` is
  computed from the same closed-form leave-one-out predictions as
  before (no refitting), and the BCa acceleration constant again uses a
  fast closed-form "delete-one" approximation rather than a full nested
  double-jackknife. As before, predicted R-squared cannot exceed the
  model's ordinary (in-sample) R-squared. The output structure
  (`fit$boot$pred.r.squared$estimate`/`conf.int`/`p.value`) is
  unchanged; only the values themselves and their computation differ.

# merya 0.7.6

* **Bug fix: `boot.glm(..., effect = "PAF")` now works without a `data`
  argument**, i.e. with formula variables taken from the calling
  environment, exactly like `glm()` and every other `effect` option
  already allowed. Previously, an internal helper
  (`.paf_recipes()`) referenced the raw `data` argument directly to
  classify each predictor as categorical or continuous; since `data`
  has no default, simply touching it when it hadn't been supplied threw
  `argument "data" is missing, with no default` (surfacing as a
  confusing secondary error during `summary()`/error-condition
  dispatch, e.g. `l'argument "data" est manquant, avec aucune valeur
  par défaut` in French R sessions). The helper now uses the model
  frame that `boot.glm()` already builds internally (which is always
  available and correctly resolved against `data`/the calling
  environment either way), with an updated, syntax-based check for
  "simple, untransformed variable name" that still correctly rejects
  interactions and transformed terms like `poly()`/`log()`. Regression
  test added.

# merya 0.7.5

* **Fixed the same misspecified-model issue as 0.7.4, this time in
  documentation, not tests.** `boot.glm()`'s `\examples{}` block (in
  both the roxygen source and `man/boot.glm.Rd`, which `R CMD check`
  executes as part of "checking examples") and `README.md` still fit a
  `binomial(link = "log")` model on `am ~ wt + hp` / `mtcars` to
  illustrate conditional risk ratios. As previously diagnosed, this
  combination fails to converge in plain `stats::glm()` on this
  particular data (predicted probabilities would need to exceed 1),
  independent of `boot.glm()`. All three now use the same
  `poisson(link = "log")` example already used in the test suite
  (`carb ~ wt`), which converges reliably and still demonstrates that
  the conditional-RR path works for any family given a log link, not
  just binomial.
* **Fixed `R CMD check` NOTE: "Namespace in Imports field not imported
  from: 'utils'".** `utils` was declared in `DESCRIPTION`'s `Imports`
  but never actually used anywhere in the package's `R/` code (a
  leftover from an earlier development stage); it has been removed.
  The package's only dependency is now `stats`. (`tests/test-merya.R`'s
  use of `capture.output()`, from `utils`, is unaffected: test scripts
  run with all default packages already attached, so they don't need a
  package `Imports` declaration the way namespaced package code does.)

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
