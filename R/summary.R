#' Summarise a matchatr fit
#'
#' @description
#' Prints the odds-ratio table of a fitted unmatched case-control logistic
#' regression — one row per model term, on the odds-ratio scale with a Wald
#' confidence interval — under a short header naming the design, estimator, and
#' variance source.
#'
#' @details
#' The intercept is shown but flagged as non-interpretable: it is not a baseline
#' risk in a case-control sample. Standard errors are model-based by default, or
#' Huber-White robust when `robust = TRUE`.
#'
#' @param object A `matchatr_fit` whose `model` is a fitted binomial `glm`.
#' @param conf.level Numeric confidence level in (0, 1). Default 0.95.
#' @param robust Logical; use the Huber-White sandwich standard error. Default
#'   `FALSE`.
#' @param ... Unused; present for S3 consistency.
#' @returns Invisibly, the odds-ratio `data.table` (as returned by
#'   [tidy.matchatr_fit()] with `exponentiate = TRUE`).
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   case = rep(c(1, 0), each = 100),
#'   x = rbinom(200, 1, 0.4),
#'   age = rnorm(200, 50, 10)
#' )
#' fit <- matcha(df, outcome = "case", exposure = "x",
#'               design = unmatched_cc(), confounders = ~ age)
#' summary(fit)
#' @seealso [matcha()], [tidy.matchatr_fit()], [contrast()]
#' @export
summary.matchatr_fit <- function(
  object,
  conf.level = 0.95,
  robust = FALSE,
  ...
) {
  require_estimated(object)
  check_conf_level(conf.level)
  or_table <- tidy(
    object,
    conf.int = TRUE,
    conf.level = conf.level,
    exponentiate = TRUE,
    robust = robust
  )

  cat("<matchatr_fit> summary\n")
  cat(" Design:     ", design_label(object$design$type), "\n", sep = "")
  cat(
    " Estimator:  ",
    object$estimator,
    "  (engine: ",
    object$engine,
    ")\n",
    sep = ""
  )
  se_label <- if (isTRUE(robust)) "robust (sandwich)" else "model-based"
  cat(" Std. error: ", se_label, "\n", sep = "")
  # The nested case-control (risk-set) design reports hazard ratios; every other
  # design tidied through this path reports odds ratios.
  scale_label <- if (identical(object$design$type, "nested_cc")) {
    "hazard ratios"
  } else {
    "odds ratios"
  }
  cat(
    "\nConditional ",
    scale_label,
    " (",
    format(100 * conf.level),
    "% Wald CI):\n",
    sep = ""
  )
  print(or_table)
  cat(
    "\nNote: the intercept is not an interpretable baseline risk in a ",
    "case-control sample.\n",
    sep = ""
  )
  invisible(or_table)
}
