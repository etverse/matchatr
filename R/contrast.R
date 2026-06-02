#' Contrast estimates from a fitted design
#'
#' @description
#' The second step of the two-step API (mirroring `causatr::contrast()` and the
#' etverse convention): given a [matcha()] fit, compute the requested causal
#' contrast on the chosen scale and return a `matchatr_result`. `matcha()`
#' resolves the sampling design and analysis; `contrast()` turns the fitted
#' estimates into a reported effect with uncertainty.
#'
#' @details
#' Producing a contrast requires fitted estimates. A fit whose estimation engine
#' has not run carries none — its `model` slot is `NULL` — so `contrast()`
#' validates its arguments and then aborts with a classed
#' `matchatr_not_estimated` error rather than returning an empty result.
#'
#' @param fit A `matchatr_fit` object returned by [matcha()].
#' @param type Character contrast scale: `"difference"` (default), `"ratio"`,
#'   or `"or"` (odds ratio).
#' @param ci_method Character confidence-interval method: `"sandwich"`
#'   (default) or `"bootstrap"`.
#' @param ... Reserved for estimator-specific contrast arguments.
#'
#' @returns A `matchatr_result` object carrying the estimates, the pairwise
#'   contrasts, and their variance-covariance matrix.
#'
#' @examples
#' df <- data.frame(case = c(1, 0, 1, 0), x = c(1, 0, 1, 0))
#' fit <- matcha(df, outcome = "case", exposure = "x", design = unmatched_cc())
#' # A design / dispatch fit carries no estimates yet, so contrast() errors:
#' try(contrast(fit))
#'
#' @seealso [matcha()]
#' @export
contrast <- function(
  fit,
  type = c("difference", "ratio", "or"),
  ci_method = c("sandwich", "bootstrap"),
  ...
) {
  if (!inherits(fit, "matchatr_fit")) {
    rlang::abort(
      "`fit` must be a `matchatr_fit` object returned by `matcha()`.",
      class = c("matchatr_bad_input", "matchatr_error")
    )
  }
  # Validate the contrast scale / CI method here so the public signature is
  # fixed regardless of which estimation engine eventually fills in the body.
  type <- match.arg(type)
  ci_method <- match.arg(ci_method)

  # `fit$model` is populated by an estimation engine; until one has run there
  # are no fitted estimates to contrast. Estimation engines replace this guard
  # with the standardisation / contrast computation that returns a
  # `matchatr_result` (via `new_matchatr_result()`).
  if (is.null(fit$model)) {
    rlang::abort(
      c(
        "This `matchatr_fit` carries no estimates to contrast.",
        i = paste0(
          "Its `",
          fit$engine,
          "` estimation engine has not produced a ",
          "fitted model."
        )
      ),
      class = c("matchatr_not_estimated", "matchatr_error")
    )
  }

  # Unreachable while every fit has `model = NULL`; an estimation engine
  # supplies the contrast computation that returns here.
  rlang::abort(
    "Internal error: `contrast()` reached estimation with no engine wired.",
    class = c("matchatr_not_estimated", "matchatr_error")
  )
}
