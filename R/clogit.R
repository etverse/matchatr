#' Fit a conditional partial likelihood for a matched or nested case-control design
#'
#' Fits the conditional maximum-likelihood (CMLE) effect for a matched or nested
#' case-control sample via [survival::clogit()]: the model is
#' `outcome ~ exposure + confounders + strata(set)`, with each matched set (or
#' sampled risk set) as a stratum. Conditioning on the per-stratum totals removes
#' the matching / risk-set nuisance parameters, so only the exposure / adjustment
#' effects are reported; the conditioning variables are controlled implicitly and
#' have no estimable coefficient. The two designs share this engine because a
#' matched set and a risk-set-sampled set are the same stratum construction; they
#' differ only in the estimand the contrast reports (a conditional odds ratio for
#' the matched design, a hazard ratio for the nested design — OR = HR exactly
#' under risk-set sampling, Prentice & Breslow 1978).
#'
#' @details
#' The conditional likelihood for a stratum with one case and M controls is
#'
#'   prod exp(x_case . beta) / sum_j exp(x_j . beta),
#'
#' which for 1:1 matching reduces to expit{(x_case - x_control) . beta}.
#' `survival::clogit` is the Cox partial likelihood with each matched set as a
#' stratum, so this is exactly the CMLE. Unconditional logistic regression on
#' matched-set indicators is never used: for 1:1 matching its MLE converges to
#' the squared odds ratio in large samples (Pike et al. 1980; Breslow & Day
#' 1980).
#'
#' The confounder formula's terms — transforms (`poly(age, 2)`), interactions
#' (`age:smoke`) — are carried through verbatim for adjustment of non-matching
#' covariates. Several matching columns are crossed into one `strata()` term
#' (frequency matching on, e.g., age group and sex). Matched sets with no case
#' or no control carry no information and are dropped by `clogit`; the
#' `matcha()` entry point already warns about them
#' (`matchatr_uninformative_stratum`). Rows with a missing outcome, exposure, or
#' confounder are dropped by the default `na.action`; a `matchatr_dropped_rows`
#' warning reports how many.
#'
#' For a nested case-control sample the sampled risk set is the stratum and the
#' design's `time` column defines how controls were drawn (incidence-density
#' sampling); the conditional likelihood reads the risk-set membership straight
#' from `strata`, so `time` is not entered in the model here — it feeds the later
#' inclusion-weight / weighted-Cox designs. Risk-set reuse (and hence a
#' cluster-robust variance) belongs to those designs, not this classical analysis
#' where each sampled set is an independent stratum.
#'
#' @param fit A `matchatr_fit` whose `engine` resolved to `"clogit"`, carrying
#'   the analysis `data`, the `outcome` / `exposure` column names, the
#'   `confounders` formula (or `NULL`), and the design's matched-set / risk-set
#'   `strata`.
#' @returns The fitted `survival::clogit` object (a `clogit` / `coxph`), or
#'   `NULL` when the design is neither matched nor nested case-control (the
#'   conditional engine has nothing wired for the other designs yet).
#' @family estimators
#' @seealso [matcha()], [contrast()], [survival::clogit()]
#' @noRd
fit_clogit <- function(fit) {
  # The conditional partial likelihood serves both the matched case-control and
  # the nested case-control (risk-set-sampled) designs: a matched set and a
  # sampled risk set are the same stratum construction. Inclusion weighting and
  # the weighted-Cox / case-cohort designs are wired by their own engines, so
  # this engine stays unestimated (NULL) for everything else.
  if (!fit$design$type %in% c("matched_cc", "nested_cc")) {
    return(NULL)
  }

  conf_terms <- if (is.null(fit$confounders)) {
    character(0)
  } else {
    attr(stats::terms(fit$confounders), "term.labels")
  }
  # `strata()` must be a bare, un-namespaced term: `survival::clogit` detects it
  # as a formula special by the name "strata", and a namespaced `survival::strata`
  # would not be recognised as the special. Several matching columns cross into
  # one stratifying factor (frequency matching).
  strata_term <- paste0(
    "strata(",
    paste(fit$design$strata, collapse = ", "),
    ")"
  )
  # With an effect modifier the exposure enters crossed with the modifier
  # (`exposure * modifier` = exposure + modifier + exposure:modifier); the
  # interaction coefficients carry the per-level shift in the exposure log OR
  # (the stratum-specific contrast). Without one the exposure is a plain main
  # effect whose coefficient(s) identify the conditional log OR.
  exposure_term <- if (is.null(fit$effect_modifier)) {
    fit$exposure
  } else {
    paste0(fit$exposure, " * ", fit$effect_modifier)
  }
  # outcome ~ exposure (* modifier) + confounders + strata(set). reformulate()
  # preserves the confounder transforms / interactions.
  model_formula <- stats::reformulate(
    termlabels = c(exposure_term, conf_terms, strata_term),
    response = fit$outcome
  )
  # Fit on a copy whose modifier is coerced to a factor with its unused levels
  # dropped: per-level odds ratios need discrete levels (and the model's
  # `xlevels`), and an empty/unused factor level would otherwise contribute an
  # all-zero interaction column aliased to NA, which would wrongly mark the
  # whole stratum-specific OR unestimable. droplevels() keeps the order of the
  # remaining declared levels, so a user-set reference level is preserved.
  fit_data <- fit$data
  em <- fit$effect_modifier
  if (!is.null(em)) {
    fit_data[[em]] <- droplevels(as.factor(fit_data[[em]]))
  }
  # `clogit` rewrites its own call to an unqualified `coxph(Surv(...) ~ ... +
  # strata(...))` and evaluates it in this frame, so those three survival names
  # must resolve here (the imports below) even though the entry point is called
  # qualified.
  model <- survival::clogit(model_formula, data = fit_data)

  # clogit's default na.action silently drops rows with a missing outcome,
  # exposure, or confounder. `model$n` is the rows actually used (coxph's
  # nobs() counts events, not rows), so the dropped count is taken against it.
  n_dropped <- nrow(fit$data) - model$n
  if (n_dropped > 0L) {
    rlang::warn(
      c(
        paste0(
          n_dropped,
          " row(s) with missing values were dropped from the fit."
        ),
        i = "The odds ratio is estimated on the complete cases only."
      ),
      class = c("matchatr_dropped_rows", "matchatr_warning")
    )
  }
  model
}

#' Assemble the conditional contrast from a clogit fit
#'
#' Turns a fitted matched or nested case-control conditional partial likelihood
#' into a `matchatr_result` reporting the exposure's effect with a Wald
#' confidence interval. The design fixes the estimand: a matched case-control
#' design reports the conditional **odds ratio** (`type = "or"`), a nested
#' case-control (risk-set-sampled) design the **hazard ratio** (`type = "hr"`).
#' Reuses the shared `conditional_or_result()` assembly (the same exp(beta) +
#' Wald-interval layer as the unmatched logistic engine); the math is identical
#' because OR = HR exactly under risk-set sampling (Prentice & Breslow 1978).
#'
#' @details
#' The variance is the inverse partial-likelihood information matrix that
#' `survival::clogit` returns; the Wald interval is on the log scale and
#' exponentiated. The risk difference and risk ratio are rejected as unidentified
#' from a case-control sample (no source-population prevalence q0), shared with
#' the unmatched engine. Each design identifies exactly one conditional scale, so
#' the off-design request (an odds ratio from a risk-set design, or a hazard ratio
#' from a matched design) names an estimand the design does not target and aborts
#' with `matchatr_unidentified_estimand`. The conditional fit needs no sandwich; a
#' cluster-robust variance (relevant only when controls are reused, handled by the
#' inclusion-weight designs) is not offered here, so `ci_method = "sandwich"` /
#' `"bootstrap"` abort with `matchatr_unsupported_variance`.
#'
#' @param fit A `matchatr_fit` whose `model` is a fitted `survival::clogit`.
#' @param type Character contrast scale; the design's native scale (`"or"` for
#'   matched, `"hr"` for nested) is computed, while the off-design scale and
#'   `"difference"` / `"ratio"` abort with `matchatr_unidentified_estimand`.
#' @param ci_method Character variance source; only `"model"` (the
#'   partial-likelihood information matrix) applies.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param call Caller environment surfaced in any error.
#' @returns A `matchatr_result` carrying the conditional odds ratio(s) (matched)
#'   or hazard ratio(s) (nested) for the exposure term.
#' @family estimators
#' @seealso [contrast()], `fit_clogit()`, `conditional_or_result()`
#' @noRd
contrast_clogit <- function(
  fit,
  type,
  ci_method,
  conf_level,
  call = rlang::caller_env()
) {
  # The conditional partial likelihood gives exp(beta); the sampling design fixes
  # its meaning. A matched case-control design targets the conditional odds
  # ratio; a risk-set-sampled nested case-control design targets the hazard ratio
  # (OR = HR exactly, no rare-disease assumption; Prentice & Breslow 1978). Each
  # design therefore identifies exactly one conditional scale.
  is_ncc <- identical(fit$design$type, "nested_cc")
  native_scale <- if (is_ncc) "hr" else "or"
  estimand <- if (is_ncc) "hazard ratio" else "conditional OR"

  # A marginal risk difference / ratio needs q0; not identified here for either
  # design. (This passes "or" / "hr" through untouched.)
  reject_unidentified_rd_rr(type, call = call)
  # The off-design conditional scale -- an odds ratio asked of a risk-set design,
  # or a hazard ratio asked of a matched design -- names an estimand the design
  # does not target, even though exp(beta) is numerically the same value.
  if (!identical(type, native_scale)) {
    reject_offdesign_conditional_scale(is_ncc = is_ncc, call = call)
  }
  # The conditional fit reports the partial-likelihood information-matrix
  # interval. A cluster-robust sandwich (for reused controls) and the bootstrap
  # belong to the inclusion-weight designs, not the conditional CMLE.
  if (!identical(ci_method, "model")) {
    rlang::abort(
      c(
        paste0(
          "`ci_method = \"",
          ci_method,
          "\"` is not available for the conditional logistic estimator."
        ),
        i = "It reports the partial-likelihood information-matrix interval; use `ci_method = \"model\"` (the default)."
      ),
      class = c("matchatr_unsupported_variance", "matchatr_error"),
      call = call
    )
  }

  # With an effect modifier the contrast is one exposure effect per modifier
  # level (the stratum-specific OR / HR), assembled from the joint
  # partial-likelihood variance; otherwise it is the single conditional effect
  # shared with the unmatched logistic engine. The design's native scale labels
  # both (the modifier gate keys on the clogit engine, which the nested design
  # also uses).
  if (!is.null(fit$effect_modifier)) {
    return(stratum_specific_or_result(
      fit,
      model = fit$model,
      conf_level = conf_level,
      ci_method = ci_method,
      type = native_scale,
      estimand = paste("stratum-specific", estimand),
      # coxph's nobs() counts events; the analysis size is the rows used (`$n`).
      n = fit$model$n,
      call = call
    ))
  }

  conditional_or_result(
    fit,
    model = fit$model,
    robust = FALSE,
    ci_method = ci_method,
    conf_level = conf_level,
    type = native_scale,
    estimand = estimand,
    # coxph's nobs() counts events; the analysis size is the rows used (`$n`).
    n = fit$model$n,
    call = call
  )
}

#' Reject an off-design conditional-scale request
#'
#' The conditional partial likelihood gives exp(beta), whose meaning is fixed by
#' the sampling design: a matched case-control design targets the conditional
#' odds ratio, a risk-set-sampled nested case-control design the hazard ratio
#' (OR = HR exactly under risk-set sampling; Prentice & Breslow 1978). Asking for
#' the *other* scale names an estimand the design does not target — a hazard
#' ratio has no meaning without risk-set sampling, and a risk-set estimate is
#' reported as a hazard ratio, not an odds ratio — so the request is rejected
#' rather than silently mislabelled.
#'
#' @param is_ncc Logical; `TRUE` when the fit is a nested case-control design
#'   (which asked for `"or"`), `FALSE` for a matched design (which asked for
#'   `"hr"`). Selects the design-specific message.
#' @param call Caller environment surfaced in the error.
#' @returns Never returns; always aborts with class
#'   `matchatr_unidentified_estimand`.
#' @family estimators
#' @seealso `contrast_clogit()`
#' @noRd
reject_offdesign_conditional_scale <- function(
  is_ncc,
  call = rlang::caller_env()
) {
  msg <- if (is_ncc) {
    c(
      "A nested case-control design is reported on the hazard-ratio scale.",
      i = paste0(
        "Risk-set (incidence-density) sampling identifies the hazard ratio ",
        "(OR = HR exactly; Prentice & Breslow 1978). Use `type = \"hr\"` ",
        "(the default)."
      )
    )
  } else {
    c(
      "A matched case-control design does not identify a hazard ratio.",
      i = paste0(
        "The hazard ratio needs risk-set (incidence-density) sampling -- a ",
        "nested case-control design (`nested_cc()`). Report the conditional ",
        "odds ratio with `type = \"or\"` (the default)."
      )
    )
  }
  rlang::abort(
    msg,
    class = c("matchatr_unidentified_estimand", "matchatr_error"),
    call = call
  )
}
