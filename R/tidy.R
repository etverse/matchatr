#' @importFrom generics tidy
#' @export
generics::tidy

#' Tidy a matchatr fit into a coefficient / odds-ratio table
#'
#' @description
#' Returns the fitted unmatched case-control logistic regression as a tidy
#' coefficient table — one row per model term — on the log-odds scale by
#' default, or the odds-ratio scale with `exponentiate = TRUE`.
#'
#' @details
#' Standard errors come from the model information matrix, or the Huber-White
#' sandwich when `robust = TRUE`. The Wald statistic and p-value are always on
#' the coefficient (log-odds) scale; with `exponentiate = TRUE` the estimate and
#' the confidence bounds are exponentiated while `std.error` stays on the
#' log-odds scale (the broom convention). The intercept row is included but is
#' not an interpretable baseline risk: under separate case / control sampling it
#' is offset by the log sampling-fraction ratio (Prentice & Pyke, 1979).
#'
#' @param x A `matchatr_fit` whose `model` is a fitted binomial `glm`.
#' @param conf.int Logical; add `conf.low` / `conf.high` Wald bounds. Default
#'   `TRUE`.
#' @param conf.level Numeric confidence level in (0, 1). Default 0.95.
#' @param exponentiate Logical; report `estimate` and confidence bounds on the
#'   odds-ratio scale. Default `FALSE`.
#' @param robust Logical; use the Huber-White sandwich standard error instead of
#'   the model information matrix. Default `FALSE`.
#' @param ... Unused; present for generic consistency.
#' @returns A `data.table` with columns `term`, `estimate`, `std.error`,
#'   `statistic`, `p.value`, and (when `conf.int`) `conf.low`, `conf.high`.
#' @examples
#' set.seed(1)
#' df <- data.frame(case = rep(c(1, 0), each = 100), x = rbinom(200, 1, 0.4))
#' fit <- matcha(df, outcome = "case", exposure = "x", design = unmatched_cc())
#' tidy(fit, exponentiate = TRUE)
#' @seealso [matcha()], [contrast()], [summary.matchatr_fit()]
#' @family tidiers
#' @export
tidy.matchatr_fit <- function(
  x,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = FALSE,
  robust = FALSE,
  ...
) {
  require_estimated(x)
  check_conf_level(conf.level)

  model <- x$model
  beta <- stats::coef(model)
  vcov_mat <- if (isTRUE(robust)) {
    sandwich::sandwich(model)
  } else {
    stats::vcov(model)
  }
  se <- sqrt(diag(vcov_mat))
  # Wald z-statistic and two-sided p-value on the log-odds scale.
  z_stat <- unname(beta / se)
  z <- stats::qnorm(1 - (1 - conf.level) / 2)

  estimate <- unname(beta)
  std_error <- unname(se)
  conf_low <- estimate - z * std_error
  conf_high <- estimate + z * std_error
  if (isTRUE(exponentiate)) {
    # Exponentiate the point estimate and bounds; std.error stays on the
    # log-odds scale per broom's convention so it is not silently delta-method
    # transformed.
    estimate <- exp(estimate)
    conf_low <- exp(conf_low)
    conf_high <- exp(conf_high)
  }

  # Assemble in one call (the package uses qualified data.table calls, never
  # `:=`, so it does not register a data.table import).
  cols <- list(
    term = names(beta),
    estimate = estimate,
    std.error = std_error,
    statistic = z_stat,
    p.value = 2 * stats::pnorm(-abs(z_stat))
  )
  if (isTRUE(conf.int)) {
    cols$conf.low <- conf_low
    cols$conf.high <- conf_high
  }
  data.table::as.data.table(cols)
}

#' Tidy a matchatr contrast result
#'
#' Returns the contrasts carried by a `matchatr_result` as
#' a tidy `data.table`, one row per reported contrast, with the contrast scale
#' recorded in a `type` column.
#'
#' @param x A `matchatr_result` object returned by [contrast()].
#' @param ... Unused; present for generic consistency.
#' @returns A `data.table` with columns `term`, `estimate`, `std.error`,
#'   `type`, `conf.low`, `conf.high`.
#' @examples
#' set.seed(1)
#' df <- data.frame(case = rep(c(1, 0), each = 100), x = rbinom(200, 1, 0.4))
#' fit <- matcha(df, outcome = "case", exposure = "x", design = unmatched_cc())
#' tidy(contrast(fit, type = "or"))
#' @seealso [contrast()]
#' @family tidiers
#' @export
tidy.matchatr_result <- function(x, ...) {
  out <- data.table::data.table(
    term = x$contrasts$comparison,
    estimate = x$contrasts$estimate,
    std.error = x$contrasts$se,
    type = x$type,
    conf.low = x$contrasts$ci_lower,
    conf.high = x$contrasts$ci_upper
  )
  out[]
}

#' Abort if a fit carries no estimated model
#'
#' Shared guard for the fit tidier / summary: an engine with no wired estimator
#' leaves `model = NULL`, so any method that reads the fitted model must reject
#' such a fit with the classed `matchatr_not_estimated` condition.
#'
#' @param x A `matchatr_fit` object.
#' @param call Caller environment surfaced in the error.
#' @returns `NULL` invisibly; aborts with class `matchatr_not_estimated` when
#'   `x$model` is `NULL`.
#' @family validators
#' @noRd
require_estimated <- function(x, call = rlang::caller_env()) {
  if (is.null(x$model)) {
    rlang::abort(
      c(
        "This `matchatr_fit` carries no estimated model.",
        i = paste0(
          "Its `",
          x$engine,
          "` estimation engine has not produced a fitted model."
        )
      ),
      class = c("matchatr_not_estimated", "matchatr_error"),
      call = call
    )
  }
  invisible(NULL)
}
