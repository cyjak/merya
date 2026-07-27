#' BCa Bootstrap Confidence Intervals for an Already-Fitted Object
#'
#' A single generic entry point that takes an object already returned by
#' \code{\link[stats]{t.test}}, \code{\link[stats]{cor.test}},
#' \code{\link[stats]{lm}}, or \code{\link[stats]{glm}}, and returns the
#' same object with its confidence interval(s) replaced by a
#' bias-corrected and accelerated (BCa) bootstrap interval and its
#' p-value(s) replaced by the CI-inversion p-value -- i.e. exactly what
#' the corresponding
#' \code{\link{boot.t.test}}/\code{\link{boot.cor.test}}/\code{\link{boot.lm}}/
#' \code{\link{boot.glm}} function would have produced had it been used
#' to fit the model in the first place.
#'
#' @param object a fitted object of class \code{"htest"} (from
#'   \code{\link[stats]{t.test}} or \code{\link[stats]{cor.test}}),
#'   \code{"lm"}, or \code{"glm"}.
#' @param ... passed on to the underlying \code{boot.*} method; see the
#'   method sections below.
#'
#' @return An object of the \strong{same class, and in the same format},
#'   as \code{object}: a \code{"htest"} list for \code{t.test}/
#'   \code{cor.test} input, or an \code{"lm"}/\code{"glm"} object for
#'   \code{lm}/\code{glm} input.
#'
#' @section How it works: For \code{"lm"} and \code{"glm"} objects,
#'   \code{boot.ci()} reads the original fitting call off the object
#'   (\code{object$call}), substitutes in the corresponding \code{boot.*}
#'   function, and re-evaluates it in the environment that called
#'   \code{boot.ci()} -- the same mechanism \code{\link[stats]{update}}
#'   uses to refit a model. This means \code{boot.ci(lm(formula, data = d))}
#'   runs the \emph{exact same} resampling/BCa/CI-inversion code as
#'   \code{boot.lm(formula, data = d)} (down to using the same internal
#'   helper functions), so with the same seed the two are equivalent; it
#'   also means the original variables (e.g. \code{d} above) must still
#'   be reachable from wherever \code{boot.ci()} is called, exactly as
#'   they would need to be for \code{update(fit)} to work.
#'
#'   \code{"htest"} objects from \code{t.test()}/\code{cor.test()} are a
#'   partial exception: base R's \code{t.test}/\code{cor.test} do not
#'   store a \code{$call}, only a deparsed \code{data.name} label (e.g.
#'   \code{"x and y"}). \code{boot.ci()} re-parses and re-evaluates that
#'   label in the calling environment by default, which works whenever
#'   \code{x}/\code{y} were plain variable names (the overwhelmingly
#'   common case); pass the data explicitly via \code{boot.ci(object, x =
#'   ..., y = ...)} to sidestep this whenever \code{x}/\code{y} were
#'   one-off expressions (e.g. \code{t.test(rnorm(10))}).
#'
#'   No new resampling machinery is introduced anywhere in this function:
#'   \code{boot.ci()} is a thin (and therefore fast) dispatcher, so its
#'   own overhead is negligible next to the cost of the bootstrap it
#'   triggers.
#'
#' @section Method-specific arguments:
#' \describe{
#'   \item{\code{htest}}{\code{x}, \code{y} (optional, to override
#'     automatic recovery from \code{data.name}), \code{R} (default
#'     \code{10000}), \code{conf.level} (default: the original test's
#'     level if recoverable, else \code{0.95}), \code{seed} (default
#'     \code{123}).}
#'   \item{\code{lm}, \code{glm}}{\code{R} (default \code{10000}),
#'     \code{conf.level} (default \code{0.95}), \code{seed} (default
#'     \code{123}); further named arguments (e.g. \code{irls.maxit},
#'     \code{irls.tol}, \code{effect}, \code{exposure} for \code{glm})
#'     are forwarded to \code{\link{boot.lm}}/\code{\link{boot.glm}}.}
#' }
#'
#' Every \code{seed} argument above is passed to
#' \code{\link[base]{set.seed}} at the start of the underlying
#' \code{boot.*} call; pass \code{seed = NULL} to use whatever random
#' state is currently active instead of reseeding.
#'
#' @examples
#' ## lm: boot.ci(lm(...)) and boot.lm(...) use the same machinery, and
#' ## agree exactly under the shared default seed (123)
#' fit <- lm(mpg ~ wt + hp, data = mtcars)
#' b1 <- boot.ci(fit, R = 300)
#' b2 <- boot.lm(mpg ~ wt + hp, data = mtcars, R = 300)
#' summary(b1)
#'
#' ## htest
#' x <- rnorm(30, mean = 1)
#' boot.ci(t.test(x))
#'
#' @seealso \code{\link{boot.t.test}}, \code{\link{boot.cor.test}},
#'   \code{\link{boot.lm}}, \code{\link{boot.glm}}
#' @export
boot.ci <- function(object, ...) UseMethod("boot.ci")

#' @export
#' @method boot.ci default
boot.ci.default <- function(object, ...) {
  stop("boot.ci() does not know how to handle an object of class ",
       paste(class(object), collapse = "/"), ". Supported input classes ",
       "are: \"htest\" (from stats::t.test()/stats::cor.test()), \"lm\", ",
       "and \"glm\".", call. = FALSE)
}

#' Take the reconstructed call for a boot.* re-fit, point it at `fn`, and
#' merge in the arguments that boot.ci() always adds/overrides
#' (`R`, `conf.level`, `seed`, and any further named arguments in
#' `extra`). Shared by every "redispatch to boot.*" method below so that
#' logic (and its one error message) isn't repeated twice.
#' @keywords internal
#' @noRd
.boot_ci_redispatch <- function(cl, fn, R, conf.level, seed, extra = list()) {
  if (is.null(cl)) {
    stop("boot.ci() requires the fitted object to carry its original ",
         "fitting call (object$call), which this object does not have.",
         call. = FALSE)
  }
  cl[[1L]] <- fn
  cl$R <- R
  cl$conf.level <- conf.level
  cl$seed <- seed
  if (length(extra)) cl[names(extra)] <- extra
  cl
}

#' @export
#' @method boot.ci htest
boot.ci.htest <- function(object, x, y, R = 10000L, conf.level = NULL,
                           seed = 123, ...) {
  meth <- object$method
  is_t   <- grepl("t-test", meth, ignore.case = TRUE)
  is_cor <- grepl("correlation", meth, ignore.case = TRUE)
  if (!is_t && !is_cor) {
    stop("boot.ci() only supports \"htest\" objects returned by ",
         "stats::t.test() or stats::cor.test() (got method: \"", meth,
         "\").", call. = FALSE)
  }

  if (is.null(conf.level)) {
    conf.level <- attr(object$conf.int, "conf.level")
    if (is.null(conf.level)) conf.level <- 0.95
  }

  ## htest objects carry no $call and no raw data, only a deparsed label
  ## (data.name); recover x/(y) from it unless supplied explicitly. This
  ## mirrors exactly how t.test()/cor.test() built data.name in the first
  ## place (deparse(substitute(x))), so it only works when x/y are plain,
  ## still-reachable variables -- pass them explicitly to sidestep this.
  have_x <- !missing(x)
  have_y <- !missing(y)
  if (!have_x) {
    parts <- strsplit(object$data.name, " and ", fixed = TRUE)[[1]]
    x <- tryCatch(eval(parse(text = parts[1])[[1]], parent.frame()),
                  error = function(e) NULL)
    if (is.null(x)) {
      stop("boot.ci() could not recover the original data for this ",
           "\"htest\" object from its data.name (\"", object$data.name,
           "\"); pass the data explicitly, e.g. ",
           "boot.ci(object, x = ..., y = ...).", call. = FALSE)
    }
    if (!have_y && length(parts) > 1L) {
      y <- tryCatch(eval(parse(text = parts[2])[[1]], parent.frame()),
                    error = function(e) NULL)
      have_y <- !is.null(y)
    }
  }

  if (is_t) {
    paired <- grepl("Paired", meth)
    mu <- object$null.value
    mu <- if (is.null(mu)) 0 else as.numeric(mu)
    boot.t.test(x, y = if (have_y) y else NULL,
                alternative = object$alternative, mu = mu, paired = paired,
                conf.level = conf.level, R = R, seed = seed)
  } else {
    if (!have_y) {
      stop("boot.ci() could not recover the second variable for this ",
           "correlation test; pass it explicitly via ",
           "boot.ci(object, x = ..., y = ...).", call. = FALSE)
    }
    cmethod <- if (grepl("Kendall", meth)) "kendall" else
      if (grepl("Spearman", meth)) "spearman" else "pearson"
    boot.cor.test(x, y, method = cmethod, alternative = object$alternative,
                   conf.level = conf.level, R = R, seed = seed)
  }
}

#' @export
#' @method boot.ci lm
boot.ci.lm <- function(object, R = 10000L, conf.level = 0.95, seed = 123, ...) {
  cl <- .boot_ci_redispatch(object$call, quote(merya::boot.lm),
                             R, conf.level, seed, list(...))
  eval(cl, envir = parent.frame())
}

#' @export
#' @method boot.ci glm
boot.ci.glm <- function(object, R = 10000L, conf.level = 0.95, seed = 123, ...) {
  cl <- .boot_ci_redispatch(object$call, quote(merya::boot.glm),
                             R, conf.level, seed, list(...))
  eval(cl, envir = parent.frame())
}
