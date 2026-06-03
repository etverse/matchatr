#' Fit an unmatched polytomous (multinomial) logistic regression
#'
#' Fits a baseline-category multinomial logistic model to a case-control sample
#' with three or more outcome groups (multiple disease subtypes, or several
#' control groups), via [nnet::multinom()]. Each non-reference outcome
#' level's equation is a binary logistic contrast against the shared reference,
#' so the exposure coefficient in that equation is the log odds ratio for that
#' subtype versus the reference (Begg & Gray, 1984; Dubin & Pasternack, 1986).
#'
#' @details
#' The model is `outcome ~ exposure + confounders` with the outcome releveled so
#' the reference group is the baseline (first) level — [nnet::multinom()]
#' contrasts every other level against that baseline. The reference is resolved
#' upstream in [matcha()] (`resolve_polytomous_outcome()`), which also rejects a
#' two-group outcome (`matchatr_bad_outcome`); here the outcome column is already
#' a reference-first factor with unused levels dropped.
#'
#' Only the slope coefficients are interpretable odds ratios: as in the binary
#' case-control logistic, separate case / control sampling offsets each
#' equation's intercept by the log sampling-fraction ratio, so the intercepts
#' are not baseline risks. Rows with a missing outcome, exposure, or confounder
#' are dropped by [nnet::multinom()]'s default `na.action`; a
#' `matchatr_dropped_rows` warning reports how many. A constant exposure is
#' rejected (`matchatr_unestimable_exposure`) because — unlike [stats::glm()] —
#' [nnet::multinom()] does not alias a collinear column to `NA`, so it would
#' otherwise return a silent, degenerate fit.
#'
#' @param fit A `matchatr_fit` whose `engine` resolved to `"multinom"`, carrying
#'   the analysis `data` (with `outcome` a reference-first factor), the
#'   `outcome` / `exposure` column names, and the `confounders` formula (or
#'   `NULL`).
#' @returns The fitted [nnet::multinom()] model object.
#' @family estimators
#' @seealso [matcha()], [contrast()]
#' @noRd
fit_polytomous <- function(fit) {
  # nnet::multinom does not alias a collinear predictor to NA (it returns a
  # degenerate fit with duplicated coefficients), so the constant-exposure case
  # that glm catches by aliasing must be rejected explicitly here.
  xcol <- fit$data[[fit$exposure]]
  if (length(unique(stats::na.omit(xcol))) < 2L) {
    rlang::abort(
      c(
        paste0("Exposure `", fit$exposure, "` is constant."),
        i = "A polytomous odds ratio needs the exposure to vary."
      ),
      class = c("matchatr_unestimable_exposure", "matchatr_error"),
      call = fit$call
    )
  }

  conf_terms <- if (is.null(fit$confounders)) {
    character(0)
  } else {
    attr(stats::terms(fit$confounders), "term.labels")
  }
  # outcome ~ exposure + confounders, with the exposure entered parametrically so
  # its per-subtype coefficient identifies the conditional log OR.
  model_formula <- stats::reformulate(
    termlabels = c(fit$exposure, conf_terms),
    response = fit$outcome
  )
  # `trace = FALSE` silences multinom's iteration log; the outcome is already a
  # reference-first factor (set in matcha()) so the baseline equation is correct.
  model <- nnet::multinom(model_formula, data = fit$data, trace = FALSE)

  # multinom has no nobs() method; the fitted residual matrix has one row per
  # complete case used. Surface listwise deletion so a missing-data problem is
  # not mistaken for the full sample (mirrors the glm engine's report).
  n_used <- nrow(model$residuals)
  n_dropped <- nrow(fit$data) - n_used
  if (n_dropped > 0L) {
    rlang::warn(
      c(
        paste0(
          n_dropped,
          " row(s) with missing values were dropped from the fit."
        ),
        i = "The odds ratios are estimated on the complete cases only."
      ),
      class = c("matchatr_dropped_rows", "matchatr_warning")
    )
  }
  model
}

#' Assemble the per-subtype odds-ratio contrast from a multinomial fit
#'
#' Turns a fitted polytomous logistic regression into a `matchatr_result`
#' reporting the exposure's odds ratio for each non-reference outcome group
#' versus the shared reference, each with a Wald confidence interval from the
#' model information matrix. The risk difference and risk ratio are rejected
#' (`matchatr_unidentified_estimand`): under separate case / control sampling the
#' marginal group frequencies are fixed by design, so only the (non-collapsible)
#' conditional ORs are identified without the source-population prevalences.
#'
#' @details
#' For a binary or continuous exposure each non-reference level contributes one
#' OR row; for an unordered factor exposure, each (subtype, exposure-level) pair
#' is a row, the exposure level taken against the factor's own reference. The
#' interval is Wald on the log-odds scale, exponentiated, so it is asymmetric on
#' the OR scale (`estimate +/- z * se` does not reproduce the bounds; the
#' OR-scale `se` is the delta-method value OR * SE(log OR)). The variance is the
#' multinomial information matrix ([nnet::multinom()]'s [stats::vcov()]); the
#' robust-sandwich and bootstrap interval methods do not apply, so they abort
#' with `matchatr_unsupported_variance`.
#'
#' @param fit A `matchatr_fit` whose `model` is a fitted [nnet::multinom()].
#' @param type Character contrast scale; only `"or"` is computed,
#'   `"difference"` / `"ratio"` abort with `matchatr_unidentified_estimand`.
#' @param ci_method Character variance source; only `"model"` (the information
#'   matrix) applies — `"sandwich"` / `"bootstrap"` abort with
#'   `matchatr_unsupported_variance`.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param call Caller environment surfaced in any error.
#' @returns A `matchatr_result` carrying the per-subtype log-OR estimates and the
#'   OR-scale contrasts for the exposure term(s).
#' @family estimators
#' @seealso [contrast()], `fit_polytomous()`
#' @noRd
contrast_polytomous <- function(
  fit,
  type,
  ci_method,
  conf_level,
  call = rlang::caller_env()
) {
  reject_unidentified_rd_rr(type, call = call)
  # The multinomial OR has one likelihood-based variance (the information
  # matrix); the model-vs-robust and bootstrap interval choices do not apply.
  if (!identical(ci_method, "model")) {
    rlang::abort(
      c(
        paste0(
          "`ci_method = \"",
          ci_method,
          "\"` is not available for the polytomous estimator."
        ),
        i = "It reports the multinomial information-matrix interval; use `ci_method = \"model\"` (the default)."
      ),
      class = c("matchatr_unsupported_variance", "matchatr_error"),
      call = call
    )
  }

  model <- fit$model
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  tab <- multinom_exposure_or(
    model,
    fit$exposure,
    conf_level = conf_level,
    call = call
  )

  estimates <- data.table::data.table(
    term = tab$comparison,
    estimate = tab$log_or,
    se = tab$se,
    ci_lower = tab$log_or - z * tab$se,
    ci_upper = tab$log_or + z * tab$se
  )
  contrasts <- data.table::data.table(
    comparison = tab$comparison,
    estimate = exp(tab$log_or), # OR
    se = exp(tab$log_or) * tab$se, # delta-method SE on the OR scale
    ci_lower = exp(tab$log_or - z * tab$se),
    ci_upper = exp(tab$log_or + z * tab$se)
  )
  vcov_exp <- tab$vcov
  dimnames(vcov_exp) <- list(tab$comparison, tab$comparison)

  new_matchatr_result(
    estimates = estimates,
    contrasts = contrasts,
    type = "or",
    estimand = "subtype OR",
    ci_method = ci_method,
    # The reference is the baseline outcome group every subtype OR is taken
    # against (multinom's first outcome level).
    reference = model$lev[1],
    n = nrow(model$residuals),
    estimator = fit$estimator,
    engine = fit$engine,
    vcov = vcov_exp,
    call = call
  )
}

#' Per-subtype exposure odds ratios from a multinomial fit
#'
#' Locates the exposure coefficient(s) in each non-reference equation of a
#' [nnet::multinom()] fit and returns their log odds ratios, standard errors, and
#' the variance sub-matrix. The exposure columns are found by term *position*
#' (via `term_assign()`), collision-free even when a factor exposure yields
#' several columns; the variance entries are read by the `level:predictor` names
#' [nnet::multinom()] gives its [stats::vcov()], so the lookup never relies on the
#' coefficient ordering.
#'
#' @param model A fitted [nnet::multinom()] object (3+ outcome levels, so its
#'   coefficients form a matrix with one row per non-reference level).
#' @param exposure Character scalar exposure column name; it must enter the model
#'   as a parametric main-effect term.
#' @param conf_level Numeric confidence level (unused for the point/SE, accepted
#'   for signature symmetry with the contrast caller).
#' @param call Caller environment surfaced in any error.
#' @returns A list with `comparison` (character `"<subtype>: <term>"` labels),
#'   `log_or` (numeric log odds ratios), `se` (their standard errors), and
#'   `vcov` (their variance sub-matrix). Aborts with `matchatr_bad_input` when
#'   the exposure is not a parametric term, or `matchatr_unestimable_exposure`
#'   when an exposure coefficient is `NA`.
#' @family estimators
#' @noRd
multinom_exposure_or <- function(
  model,
  exposure,
  conf_level = 0.95,
  call = rlang::caller_env()
) {
  cf <- stats::coef(model) # matrix: rows = non-reference levels, cols = predictors
  vc <- stats::vcov(model) # rows/cols named "<level>:<predictor>"
  subtypes <- rownames(cf)
  predictors <- colnames(cf)

  # Exposure columns by term position (collision-free), reusing the shared
  # parametric term map; multinom exposes model.matrix()/terms() like a glm.
  ta <- term_assign(model)
  pos <- match(exposure, ta$labels)
  if (is.na(pos)) {
    rlang::abort(
      paste0(
        "Exposure `",
        exposure,
        "` is not a parametric main-effect term in the fitted model."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  exp_cols <- which(ta$assign == pos)

  # One OR per (subtype, exposure-column); the variance entry is matched by the
  # "<subtype>:<predictor>" vcov name, so it is robust to the coefficient
  # ordering rather than assuming a row-major layout.
  comparison <- character(0)
  log_or <- numeric(0)
  vnames <- character(0)
  for (k in seq_along(subtypes)) {
    for (j in exp_cols) {
      comparison <- c(comparison, paste0(subtypes[k], ": ", predictors[j]))
      log_or <- c(log_or, cf[k, j])
      vnames <- c(vnames, paste0(subtypes[k], ":", predictors[j]))
    }
  }
  # A degenerate (e.g. perfectly separated) equation can leave a coefficient NA;
  # refuse rather than report a silent NA odds ratio.
  if (anyNA(log_or)) {
    rlang::abort(
      c(
        paste0("Exposure `", exposure, "` has no estimable effect."),
        i = "It is constant or collinear with the confounders, so its odds ratio is not identified."
      ),
      class = c("matchatr_unestimable_exposure", "matchatr_error"),
      call = call
    )
  }
  vcov_sel <- vc[vnames, vnames, drop = FALSE]
  list(
    comparison = comparison,
    log_or = unname(log_or),
    se = unname(sqrt(diag(vcov_sel))),
    vcov = vcov_sel
  )
}

#' Tidy a multinomial fit into a per-equation coefficient table
#'
#' Builds the tidy coefficient table for a [nnet::multinom()] fit: one row per
#' (outcome level, model term), with a `y.level` column naming the non-reference
#' outcome group each equation contrasts against the reference (the broom
#' convention). Used by [tidy.matchatr_fit()] when the fitted engine is
#' polytomous, where the coefficients form a matrix rather than the single vector
#' the glm tidier assumes.
#'
#' @param model A fitted [nnet::multinom()] object.
#' @param conf.int Logical; add `conf.low` / `conf.high` Wald bounds.
#' @param conf.level Numeric confidence level in (0, 1).
#' @param exponentiate Logical; report `estimate` and bounds on the odds-ratio
#'   scale (`std.error` stays on the log-odds scale, per broom).
#' @returns A `data.table` with columns `y.level`, `term`, `estimate`,
#'   `std.error`, `statistic`, `p.value`, and (when `conf.int`) `conf.low`,
#'   `conf.high`.
#' @family tidiers
#' @noRd
tidy_multinom <- function(
  model,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = FALSE
) {
  cf <- stats::coef(model)
  vc <- stats::vcov(model)
  subtypes <- rownames(cf)
  predictors <- colnames(cf)
  z <- stats::qnorm(1 - (1 - conf.level) / 2)

  # Walk the coefficient matrix row (subtype) by row, reading each entry's SE
  # from the "<subtype>:<predictor>" vcov name so the table never assumes a
  # particular coefficient ordering.
  y_level <- character(0)
  term <- character(0)
  estimate <- numeric(0)
  std_error <- numeric(0)
  for (k in seq_along(subtypes)) {
    for (j in seq_along(predictors)) {
      y_level <- c(y_level, subtypes[k])
      term <- c(term, predictors[j])
      estimate <- c(estimate, cf[k, j])
      vname <- paste0(subtypes[k], ":", predictors[j])
      std_error <- c(std_error, sqrt(vc[vname, vname]))
    }
  }
  statistic <- estimate / std_error
  conf_low <- estimate - z * std_error
  conf_high <- estimate + z * std_error
  if (isTRUE(exponentiate)) {
    estimate <- exp(estimate)
    conf_low <- exp(conf_low)
    conf_high <- exp(conf_high)
  }

  cols <- list(
    y.level = y_level,
    term = term,
    estimate = estimate,
    std.error = std_error,
    statistic = statistic,
    p.value = 2 * stats::pnorm(-abs(statistic))
  )
  if (isTRUE(conf.int)) {
    cols$conf.low <- conf_low
    cols$conf.high <- conf_high
  }
  data.table::as.data.table(cols)
}
