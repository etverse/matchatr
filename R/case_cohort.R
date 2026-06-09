#' Fit a case-cohort Cox pseudo-likelihood
#'
#' Fits the Prentice, Self-Prentice, or Lin-Ying pseudo-likelihood hazard
#' ratio for a case-cohort sample via [survival::cch()]. The case-cohort
#' design draws a random subcohort from the full cohort at baseline and
#' augments it with every case; controls are reused across failure times, so
#' the pseudo-likelihood's dependent score factors mean the variance does NOT
#' come from the naive information matrix. `survival::cch` returns the correct
#' asymptotic (Self-Prentice) or robust (Lin-Ying) variance; matchatr passes
#' it through unchanged.
#'
#' @details
#' The three unstratified methods differ in which "risk set" backs each failure
#' time:
#' - `"Prentice"` uses the subcohort members at risk plus the case if the case
#'   is outside the subcohort. Score-unbiased.
#' - `"SelfPrentice"` uses only the subcohort members at risk. Same asymptotic
#'   variance as Prentice (Self & Prentice 1988).
#' - `"LinYing"` uses the full subcohort plus all failures. The robust variance
#'   estimator (Lin & Ying 1993) is returned by default; the conservative
#'   design-based SE is available via `robust = TRUE` inside `cch`.
#'
#' `matcha()` is called with the FULL cohort; `fit_cch()` subsets internally
#' to cases plus subcohort members before passing to `survival::cch`. The
#' cohort size used in the pseudo-likelihood denominator is `nrow(data)` (the
#' full cohort passed to `matcha()`).
#'
#' Missing values in the outcome, exposure, or confounders are handled by
#' `survival::cch`'s default `na.action`; a `matchatr_dropped_rows` warning
#' reports how many rows were dropped.
#'
#' @param fit A `matchatr_fit` whose `engine` resolved to `"cch"`, carrying
#'   the analysis `data` (full cohort), the `outcome` / `exposure` column
#'   names, the `confounders` formula (or `NULL`), and the design's `subcohort`,
#'   `time`, `method`, and (optionally) `id` slots.
#' @returns A fitted `survival::cch` object.
#' @family estimators
#' @seealso [matcha()], [contrast()], [survival::cch()]
#' @noRd
fit_cch <- function(fit) {
  dt <- fit$data
  outcome_col <- fit$outcome
  subcohort_col <- fit$design$subcohort
  time_col <- fit$design$time
  method <- fit$design$method %||% "Prentice"
  id_col <- fit$design$id

  n_cohort <- nrow(dt)

  # Subset to cases + subcohort members; censored non-subcohort subjects
  # contribute nothing to the pseudo-likelihood and survival::cch rejects them.
  sc_vec <- dt[[subcohort_col]]
  # subcohort column may be logical or numeric 0/1; treat any non-zero as TRUE.
  is_sc <- if (is.logical(sc_vec)) sc_vec else as.logical(sc_vec != 0L)
  y01 <- dt[[outcome_col]] # already resolved to integer 0/1 by matcha()
  in_sample <- is_sc | (y01 == 1L)
  subset_dt <- dt[in_sample, ]

  # Subject IDs for survival::cch: needed to correctly pair subjects that appear
  # as both subcohort member and case. Use the named column when supplied; fall
  # back to the original row indices in the full cohort.
  if (!is.null(id_col)) {
    id_vec <- subset_dt[[id_col]]
  } else {
    id_vec <- which(in_sample)
  }

  # Build the formula RHS from exposure + confounders.
  conf_terms <- if (is.null(fit$confounders)) {
    character(0)
  } else {
    attr(stats::terms(fit$confounders), "term.labels")
  }
  rhs_terms <- c(fit$exposure, conf_terms)

  # Build the Surv response as a string; survival::cch expects a right-censored
  # Surv object. reformulate() accepts a character string for the response.
  model_formula <- stats::reformulate(
    termlabels = rhs_terms,
    response = paste0("survival::Surv(", time_col, ", ", outcome_col, ")")
  )

  model <- survival::cch(
    formula = model_formula,
    data = subset_dt,
    subcoh = stats::as.formula(paste0("~", subcohort_col)),
    id = id_vec,
    cohort.size = n_cohort,
    method = method
  )

  # survival::cch may silently drop rows with missing values; warn if so.
  n_dropped <- length(in_sample[in_sample]) - model$n
  if (n_dropped > 0L) {
    rlang::warn(
      c(
        paste0(
          n_dropped,
          " row(s) with missing values were dropped from the fit."
        ),
        i = "The hazard ratio is estimated on the complete cases only."
      ),
      class = c("matchatr_dropped_rows", "matchatr_warning")
    )
  }
  model
}

#' Assemble the case-cohort hazard-ratio contrast from a cch fit
#'
#' Turns a fitted case-cohort pseudo-likelihood into a `matchatr_result`
#' reporting the exposure's hazard ratio with a Wald confidence interval.
#' The variance returned by [survival::cch()] is the correct asymptotic
#' variance for the method used (Prentice / Self-Prentice share the
#' Self-Prentice correction; Lin-Ying uses its own variance). The naive
#' information-matrix SE is never used.
#'
#' @details
#' The interval is Wald on the log scale and exponentiated, so it is
#' asymmetric on the HR scale: `estimate +/- z * se` does not reproduce
#' `ci_lower` / `ci_upper`. The risk difference and risk ratio are rejected as
#' unidentified from a case-cohort sample (no source-population prevalence q0).
#' Only `ci_method = "model"` applies (the pseudo-likelihood asymptotic
#' variance); sandwich and bootstrap variants of the pseudo-likelihood are not
#' offered here.
#'
#' @param fit A `matchatr_fit` whose `model` is a fitted `survival::cch`.
#' @param type Character contrast scale; `"hr"` (the default for this design)
#'   is computed, while `"or"`, `"difference"`, and `"ratio"` abort with
#'   `matchatr_unidentified_estimand`.
#' @param ci_method Character variance source; only `"model"` (the cch
#'   asymptotic variance) applies.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param call Caller environment surfaced in any error.
#' @returns A `matchatr_result` carrying the hazard ratio(s) for the exposure
#'   term, with the cch asymptotic variance.
#' @family estimators
#' @seealso [contrast()], `fit_cch()`
#' @noRd
contrast_cch <- function(
  fit,
  type,
  ci_method,
  conf_level,
  call = rlang::caller_env()
) {
  # A marginal risk difference / ratio needs q0; not identified here.
  reject_unidentified_rd_rr(type, call = call)
  # OR is rejected: a case-cohort design with risk-set pseudo-likelihood targets
  # the hazard ratio. Requesting an odds ratio confuses the estimand.
  if (identical(type, "or")) {
    rlang::abort(
      c(
        "A case-cohort design is reported on the hazard-ratio scale.",
        i = paste0(
          "The Prentice / Self-Prentice / Lin-Ying pseudo-likelihood identifies ",
          "the hazard ratio. Use `type = \"hr\"` (the default)."
        )
      ),
      class = c("matchatr_unidentified_estimand", "matchatr_error"),
      call = call
    )
  }
  if (!identical(ci_method, "model")) {
    rlang::abort(
      c(
        paste0(
          "`ci_method = \"",
          ci_method,
          "\"` is not available for the case-cohort estimator."
        ),
        i = paste0(
          "It reports the pseudo-likelihood asymptotic variance; ",
          "use `ci_method = \"model\"` (the default)."
        )
      ),
      class = c("matchatr_unsupported_variance", "matchatr_error"),
      call = call
    )
  }

  model <- fit$model
  beta_vec <- stats::coef(model)
  vcov_mat <- stats::vcov(model)
  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  # Locate the exposure coefficient(s) by building a model.matrix from the
  # user-facing formula against the data; the column names in that matrix
  # exactly match the cch coefficient names (both use standard R contrasts).
  idx_names <- cch_exposure_coef_names(fit)
  if (length(idx_names) == 0L) {
    rlang::abort(
      c(
        paste0("Exposure `", fit$exposure, "` has no estimable effect."),
        i = paste0(
          "It is constant or collinear with the confounders, so ",
          "its hazard ratio is not identified."
        )
      ),
      class = c("matchatr_unestimable_exposure", "matchatr_error"),
      call = call
    )
  }
  # Guard: coefficient present in cch output (a constant or perfectly collinear
  # exposure would be absent from cch's result, unlike glm which aliases to NA).
  missing_coefs <- setdiff(idx_names, names(beta_vec))
  if (length(missing_coefs) > 0L) {
    rlang::abort(
      c(
        paste0("Exposure `", fit$exposure, "` has no estimable effect."),
        i = paste0(
          "It is constant or collinear with the confounders, so ",
          "its hazard ratio is not identified."
        )
      ),
      class = c("matchatr_unestimable_exposure", "matchatr_error"),
      call = call
    )
  }

  b <- unname(beta_vec[idx_names])
  v <- vcov_mat[idx_names, idx_names, drop = FALSE]
  s <- unname(sqrt(diag(v)))
  log_lower <- b - z_crit * s
  log_upper <- b + z_crit * s

  # For a factor exposure, record the reference (baseline) level.
  exposure_col <- fit$data[[fit$exposure]]
  reference <- if (is.factor(exposure_col)) {
    xl <- model$xlevels[[fit$exposure]]
    if (is.null(xl)) levels(droplevels(exposure_col))[1L] else xl[1L]
  } else {
    NULL
  }

  estimates <- data.table::data.table(
    term = idx_names,
    estimate = b,
    se = s,
    ci_lower = log_lower,
    ci_upper = log_upper
  )
  contrasts <- data.table::data.table(
    comparison = idx_names,
    estimate = exp(b),
    se = exp(b) * s, # delta-method SE on the exponentiated scale
    ci_lower = exp(log_lower),
    ci_upper = exp(log_upper)
  )

  new_matchatr_result(
    estimates = estimates,
    contrasts = contrasts,
    type = "hr",
    estimand = "hazard ratio",
    ci_method = ci_method,
    reference = reference,
    n = model$n,
    estimator = fit$estimator,
    engine = fit$engine,
    vcov = v,
    call = call
  )
}

#' Coefficient names for the exposure term in a cch fit
#'
#' Builds a model matrix from the user-facing RHS formula (exposure +
#' confounders) against the fit's data to obtain the standard R coefficient
#' names for the exposure term. This is needed because `survival::cch`
#' rewrites its internal formula and `assign` structure in a non-standard form,
#' so the usual `term_assign()` approach does not apply; the column names from
#' a fresh `model.matrix()` call match the `cch` coefficient names exactly.
#'
#' @param fit A `matchatr_fit` with `engine = "cch"`.
#' @returns A character vector of coefficient names corresponding to the
#'   exposure term (length 1 for binary/continuous, length k-1 for a factor
#'   with k levels). Returns `character(0)` if the exposure is not found.
#' @family estimators
#' @noRd
cch_exposure_coef_names <- function(fit) {
  conf_terms <- if (is.null(fit$confounders)) {
    character(0)
  } else {
    attr(stats::terms(fit$confounders), "term.labels")
  }
  rhs_terms <- c(fit$exposure, conf_terms)
  rhs_formula <- stats::reformulate(rhs_terms)

  # Build the model matrix from the full data (the cch subset shares the same
  # column types / factor levels); suppress the intercept via [-1] to avoid
  # issues when assign == 0.
  mm <- tryCatch(
    stats::model.matrix(rhs_formula, data = fit$data),
    error = function(e) NULL
  )
  if (is.null(mm)) {
    return(character(0))
  }

  assign_attr <- attr(mm, "assign")
  term_labels <- attr(stats::terms(rhs_formula), "term.labels")
  term_pos <- match(fit$exposure, term_labels)
  if (is.na(term_pos)) {
    return(character(0))
  }
  # Position 0 is the intercept; positions >= 1 are the model terms in order.
  coef_col_idx <- which(assign_attr == term_pos)
  colnames(mm)[coef_col_idx]
}
