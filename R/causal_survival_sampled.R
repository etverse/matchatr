#' Fit a design-weighted marginal causal-survival model for a sampled cohort
#'
#' Estimates marginal causal-survival contrasts (absolute risk, risk difference,
#' risk ratio, restricted mean survival time) from a nested case-control or
#' case-cohort sample by g-computation on the design's own weighted Cox model:
#' the design-faithful weighted partial likelihood gives the log hazard ratio,
#' an inverse-probability-weighted Breslow baseline gives the cumulative hazard,
#' and the subject-specific absolute risk F(t | x, W) is standardized over the
#' cohort covariate distribution. The result is the full-cohort marginal effect,
#' not the sampled subset's.
#'
#' @details
#' The fit step delegates to the design's classical weighted-Cox estimator —
#' `fit_cch()` (case-cohort Prentice / Self-Prentice / Borgan pseudo-likelihood)
#' or `fit_ipw_cox()` (nested case-control Samuelsen IPW weighted Cox) — both of
#' which handle the sampled risk-set structure correctly (the weights enter the
#' partial likelihood as the design intends, not as constant person-period
#' weights). The g-computation contrast is assembled later by
#' `contrast_surv_gcomp()`, which marginalizes the engine's `absolute_risk()`
#' over the cohort covariate distribution.
#'
#' Standardizing F(t | x, W) over the covariate distribution is what makes the
#' constant inclusion weights valid here: the marginalization is a
#' Horvitz-Thompson population mean (each sampled subject's weight is its inverse
#' inclusion probability), which is unbiased for the cohort mean of the
#' (correctly fitted) risk function — unlike feeding the same constant weights
#' into a discrete-time hazard model, where they misrepresent the time-varying
#' risk sets.
#'
#' @param fit A `matchatr_fit` whose `engine` resolved to `"surv_gcomp"`, with a
#'   `nested_cc()` or `case_cohort()` design, a binary `outcome`, a binary
#'   `exposure`, and the design `time` column (plus the `ipw_weight` /
#'   `.cohort_row` columns from `sample_ncc(incl_prob = TRUE)` for a nested
#'   case-control design).
#' @returns The fitted weighted-Cox model (a `survival::cch` for a case-cohort
#'   design, a `survival::coxph` for a nested case-control design), stored in the
#'   `matchatr_fit`'s `model` slot; [contrast()] turns it into the marginal
#'   survival contrast.
#' @family causal survival
#' @seealso [matcha()], [contrast()], [absolute_risk()], `fit_cch()`,
#'   `fit_ipw_cox()`
#' @noRd
fit_surv_gcomp <- function(fit) {
  fit <- surv_gcomp_recode_exposure(fit)
  switch(
    fit$design$type,
    case_cohort = fit_cch(fit),
    nested_cc = fit_ipw_cox(fit),
    # Defensive: the dispatch table only routes case_cohort / nested_cc here.
    rlang::abort(
      paste0(
        "Design-weighted causal survival is not available for design `",
        fit$design$type,
        "`."
      ),
      class = c("matchatr_bad_design", "matchatr_error")
    )
  )
}

#' Recode the exposure to 0/1 for a design-weighted survival fit
#'
#' Replaces the exposure column with its 0/1 coding so the treat-all
#' (X = 1) / treat-none (X = 0) marginalization sets the column to a value the
#' fitted model matrix understands, and the coefficient is on the binary scale.
#' A non-binary exposure has no binary average-treatment-effect contrast and is
#' rejected here, so the guard fires at `matcha()` time.
#'
#' @param fit A `matchatr_fit` with a `surv_gcomp` engine.
#' @returns The `matchatr_fit` with a `data.table` copy of `data` whose exposure
#'   column is the 0/1 coding; aborts `matchatr_bad_input` for a non-binary
#'   exposure.
#' @family causal survival
#' @noRd
surv_gcomp_recode_exposure <- function(fit) {
  x01 <- resolve_binary_exposure(
    fit$data,
    fit$exposure,
    estimator_label = "design-weighted causal survival",
    alternative = "a conditional estimator (e.g. estimator = \"cch\" or \"ipw_cox\")"
  )
  dt <- data.table::copy(data.table::as.data.table(fit$data))
  dt[[fit$exposure]] <- x01
  fit$data <- dt
  fit
}

#' Underlying weighted-Cox engine for a design-weighted survival fit
#'
#' The `surv_gcomp` engine standardizes the design's classical weighted-Cox
#' absolute risk; this maps the sampling design to that engine so the
#' `absolute_risk()` dispatch (which keys on the engine) reaches the matching
#' inverse-probability Breslow baseline.
#'
#' @param design_type Character scalar design type.
#' @returns The underlying engine key (`"cch"` or `"ipw_cox"`).
#' @family causal survival
#' @noRd
surv_gcomp_underlying_engine <- function(design_type) {
  switch(
    design_type,
    case_cohort = "cch",
    nested_cc = "ipw_cox",
    rlang::abort(
      paste0("No weighted-Cox engine for design `", design_type, "`."),
      class = c("matchatr_bad_design", "matchatr_error")
    )
  )
}

#' Standardization sample and inclusion weights for the marginalization
#'
#' Returns the covariate rows over which the subject-specific absolute risk is
#' averaged to obtain the marginal risk, together with each row's inclusion
#' weight (so the average is the Horvitz-Thompson estimate of the cohort mean).
#' A case-cohort design marginalizes over its subcohort — a random sample of the
#' cohort, weighted by the inverse subcohort sampling fraction `N_s / n_s`,
#' stratum-specific when sampled within strata. A nested case-control design
#' marginalizes over the deduplicated risk-set sample weighted by the Samuelsen
#' inclusion weights `1/π_j` (cases at weight 1).
#'
#' @param fit A `matchatr_fit` with a `nested_cc()` or `case_cohort()` design.
#' @returns A list with `newdata` (a data frame of covariate rows), `weights`
#'   (the inclusion weights aligned to `newdata`), and `engine` (the underlying
#'   weighted-Cox engine for the `absolute_risk()` dispatch).
#' @family causal survival
#' @seealso `contrast_surv_gcomp()`, `subcohort_std_weights()`
#' @noRd
surv_gcomp_std_sample <- function(fit) {
  design_type <- fit$design$type
  under <- surv_gcomp_underlying_engine(design_type)
  dt <- data.table::as.data.table(fit$data)

  if (identical(design_type, "nested_cc")) {
    require_ipw_ncc_columns(fit, "surv_gcomp")
    samp <- data.table::as.data.table(ncc_ipw_analysis_data(fit))
    return(list(
      newdata = as.data.frame(samp),
      weights = as.numeric(samp[["ipw_weight"]]),
      engine = under
    ))
  }

  if (is.null(fit$design$time)) {
    rlang::abort(
      c(
        "The design-weighted causal-survival fit requires a `time` column in `case_cohort()`.",
        i = "Use `case_cohort(subcohort = \"sub\", time = \"t\")` with the follow-up time column."
      ),
      class = c("matchatr_bad_design", "matchatr_error")
    )
  }
  sc <- dt[[fit$design$subcohort]]
  is_sc <- if (is.logical(sc)) sc else as.logical(sc != 0L)
  # The subcohort is the random cohort sample; weight each member by the inverse
  # of its (stratum-specific) subcohort sampling fraction so the average is the
  # cohort covariate mean.
  w_sub <- subcohort_std_weights(dt, fit$design, is_sc)
  list(
    newdata = as.data.frame(dt[is_sc, , drop = FALSE]),
    weights = w_sub[is_sc],
    engine = under
  )
}

#' Inverse subcohort-sampling-fraction weights for subcohort members
#'
#' Each subcohort member represents `N_s / n_s` cohort subjects in its sampling
#' stratum (the whole cohort when the subcohort is a simple random sample), so
#' that is its weight in the Horvitz-Thompson covariate-distribution average.
#'
#' @param dt A `data.table` of the full cohort.
#' @param design The `matchatr_design` (its `subcohort` and `stratum` slots).
#' @param is_sc Logical vector flagging subcohort membership, aligned to `dt`.
#' @returns A numeric vector of standardization weights, one per row of `dt`
#'   (the value is only consumed for subcohort members).
#' @family causal survival
#' @seealso `surv_gcomp_std_sample()`
#' @noRd
subcohort_std_weights <- function(dt, design, is_sc) {
  n <- nrow(dt)
  stratum_cols <- design$stratum
  if (is.null(stratum_cols)) {
    return(rep(n / sum(is_sc), n))
  }
  strat <- if (length(stratum_cols) == 1L) {
    as.character(dt[[stratum_cols]])
  } else {
    do.call(paste, c(dt[, stratum_cols, with = FALSE], sep = ":"))
  }
  big_n_s <- table(strat)
  small_n_s <- table(strat[is_sc])
  as.numeric(big_n_s[strat]) / as.numeric(small_n_s[strat])
}
