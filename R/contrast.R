#' Contrast estimates from a fitted design
#'
#' @description
#' The second step of the two-step API (mirroring `causatr::contrast()` and the
#' etverse convention): given a [matcha()] fit, compute the requested effect on
#' the chosen scale and return a `matchatr_result`. `matcha()` resolves the
#' sampling design and runs the analysis; `contrast()` turns the fitted
#' estimates into a reported effect with uncertainty.
#'
#' @details
#' What is identifiable depends on the design. From an unmatched case-control
#' sample only the conditional odds ratio (`type = "or"`) is identified: under
#' separate case / control sampling the marginal outcome frequency is fixed, so
#' a marginal risk difference (`type = "difference"`) or risk ratio
#' (`type = "ratio"`) requires the source-population prevalence q0 and a
#' case-control-weighted estimator. Requesting an unidentified estimand aborts
#' with the classed `matchatr_unidentified_estimand` condition.
#'
#' A fit whose engine has no wired estimator carries no estimates — its `model`
#' slot is `NULL` — so `contrast()` validates its arguments and then aborts with
#' `matchatr_not_estimated`.
#'
#' @param fit A `matchatr_fit` object returned by [matcha()].
#' @param type Character contrast scale: `"difference"`, `"ratio"`, or `"or"`
#'   (odds ratio). When omitted, it defaults to the estimand the design
#'   identifies — `"or"` for the classical odds-ratio engines (unmatched
#'   case-control logistic), `"difference"` otherwise.
#' @param ci_method Character variance source for the interval: `"model"`
#'   (information-matrix Wald, the default), `"sandwich"` (Huber-White robust),
#'   or `"bootstrap"`.
#' @param conf_level Numeric confidence level for the interval, a single number
#'   strictly in (0, 1). Defaults to 0.95.
#' @param ... Reserved for estimator-specific contrast arguments.
#'
#' @returns A `matchatr_result` object carrying the estimates, the contrasts,
#'   and their variance-covariance matrix.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   case = rep(c(1, 0), each = 100),
#'   x = rbinom(200, 1, 0.4),
#'   age = rnorm(200, 50, 10)
#' )
#' fit <- matcha(df, outcome = "case", exposure = "x",
#'               design = unmatched_cc(), confounders = ~ age)
#' # The conditional odds ratio is identified:
#' contrast(fit, type = "or")
#' # The risk difference is not (no prevalence q0):
#' try(contrast(fit, type = "difference"))
#'
#' @seealso [matcha()], [tidy.matchatr_fit()]
#' @export
contrast <- function(
  fit,
  type = c("difference", "ratio", "or"),
  ci_method = c("model", "sandwich", "bootstrap"),
  conf_level = 0.95,
  ...
) {
  if (!inherits(fit, "matchatr_fit")) {
    rlang::abort(
      "`fit` must be a `matchatr_fit` object returned by `matcha()`.",
      class = c("matchatr_bad_input", "matchatr_error")
    )
  }
  # Validate the contrast scale / CI method / confidence level up front so the
  # public signature is fixed regardless of which engine fills in the body.
  # When the caller does not name a `type`, default to the estimand the design
  # identifies (the OR for an odds-ratio-only engine) rather than the generic
  # risk difference, which such an engine would have to reject.
  if (missing(type)) {
    type <- default_contrast_type(fit$engine)
  } else {
    type <- match.arg(type)
  }
  ci_method <- match.arg(ci_method)
  check_conf_level(conf_level)

  # An engine with no wired estimator leaves `model = NULL`; there is nothing to
  # contrast until its estimation layer lands.
  if (is.null(fit$model)) {
    rlang::abort(
      c(
        "This `matchatr_fit` carries no estimates to contrast.",
        i = paste0(
          "Its `",
          fit$engine,
          "` estimation engine has not produced a fitted model."
        )
      ),
      class = c("matchatr_not_estimated", "matchatr_error")
    )
  }

  call <- match.call()
  switch(
    fit$engine,
    glm_logistic = contrast_logistic(
      fit,
      type = type,
      ci_method = ci_method,
      conf_level = conf_level,
      call = call
    ),
    # Defensive: a fitted model whose engine has no contrast assembly wired.
    rlang::abort(
      paste0(
        "No contrast is wired for the `",
        fit$engine,
        "` engine."
      ),
      class = c("matchatr_not_estimated", "matchatr_error")
    )
  )
}
