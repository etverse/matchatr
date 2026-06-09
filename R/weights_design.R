#' Compute working-model inclusion probability weights for an NCC sample
#'
#' @description
#' Replaces or sets the `ipw_weight` column of a nested case-control dataset
#' using a **working model** for the probability that a subject is selected as
#' a control at each event time. Two estimation strategies are available:
#'
#' * `"glm"` — fits a logistic regression of the binary selection indicator
#'   (`selected ~ selection_formula`) across all (subject, event-time) pairs
#'   where the subject is in the risk set; uses [stats::glm()].
#' * `"gam"` — same but uses [mgcv::gam()] so that smooth terms
#'   (e.g. `s(risk_time)`) are available in `selection_formula`.
#'
#' Both methods require the **full Phase-1 cohort** to reconstruct the risk set
#' at every event time and determine which subjects were eligible but not
#' selected. Omitting `cohort` (or supplying a cohort without the `time`
#' column) aborts with `matchatr_missing_phase1`.
#'
#' @details
#' The working-model inclusion probability for subject j is:
#'
#'   π_j = 1 − ∏_{k: j ∈ R(t_k)} (1 − p̂_jk)
#'
#' where the product runs over all event times t_k at which j was in the risk
#' set R(t_k) (excluding the case's own failure time), and p̂_jk is the
#' predicted probability from the fitted selection model. Cases (cohort events)
#' are always included with weight 1.
#'
#' The default `selection_formula = NULL` maps to `~ risk_time`, a time-only
#' logistic model matching the simplest GLM specification of Borgan, Samuelsen
#' & Aastveit (2003, *Lifetime Data Analysis* 9(2)). Richer models can include
#' cohort covariates: `~ risk_time + age + sex`.
#'
#' **Population-stratum matching caveat.** If `ncc` was generated with
#' [sample_ncc()] using a `match` argument (e.g. `match = ~ sex`), controls
#' were drawn only from subjects in the case's matching stratum. The default
#' working model does not condition on the matching variable, which causes all
#' at-risk subjects (including those from other strata) to be treated as
#' "eligible but not selected". Add the matching variable to `selection_formula`
#' (e.g. `selection_formula = ~ risk_time + sex`) to partially account for this.
#'
#' @param ncc A `data.table` or `data.frame` from [sample_ncc()] with
#'   `incl_prob = TRUE`. Must contain the columns `.cohort_row` (integer cohort
#'   row index), `ipw_weight`, `set` (matched-set id), `case` (per-set 0/1
#'   indicator), and `risk_time` (set's event time).
#' @param cohort A `data.frame` or `data.table` — the **full Phase-1 cohort**
#'   from which `ncc` was drawn. Must contain the `time` column and any
#'   covariates referenced in `selection_formula`. Required for working-model
#'   methods; `NULL` is rejected with `matchatr_missing_phase1`.
#' @param method Character scalar; the weight estimation strategy. One of
#'   `"glm"` (logistic regression via [stats::glm()]) or `"gam"` (generalised
#'   additive model via [mgcv::gam()]).
#' @param selection_formula `NULL` (the default) or a one-sided formula naming
#'   the predictors for the selection model. The response is always `selected`
#'   (the binary indicator of being chosen as a control). Default `NULL` maps
#'   to `~ risk_time`. Smooth terms (e.g. `s(risk_time)`) are supported only
#'   for `method = "gam"`.
#' @param time A single character string naming the **event/follow-up time**
#'   column in `cohort`. Determines which cohort subjects are at risk at each
#'   event time: a subject is eligible if `cohort[[time]] >= risk_time_k` (and
#'   `cohort[[entry]] < risk_time_k` when `entry` is supplied).
#' @param entry `NULL` (everyone enters at the origin) or a single character
#'   string naming a **delayed-entry** column in `cohort`.
#'
#' @returns A copy of `ncc` (as a `data.table`) with the `ipw_weight` column
#'   set to the working-model inverse inclusion probabilities: 1 for cases,
#'   1/π_j (>= 1) for sampled controls. The `.cohort_row` column is preserved
#'   unchanged.
#'
#' @examples
#' set.seed(1)
#' cohort <- data.frame(
#'   id = 1:500,
#'   t  = rexp(500, 0.1),
#'   d  = rbinom(500, 1, 0.15),
#'   x  = rbinom(500, 1, 0.4)
#' )
#' ncc <- sample_ncc(cohort, time = "t", event = "d", m = 2, incl_prob = TRUE)
#'
#' # Replace Samuelsen KM weights with GLM working-model weights
#' ncc_glm <- compute_ncc_weights(ncc, cohort = cohort, method = "glm", time = "t")
#'
#' # Fit IPW Cox using the GLM weights
#' fit <- matcha(ncc_glm, outcome = "d", exposure = "x",
#'               design = nested_cc(strata = "set", time = "t"),
#'               estimator = "ipw_cox")
#' contrast(fit)
#'
#' @family sampling
#' @seealso [sample_ncc()], [matcha()], [nested_cc()]
#' @export
compute_ncc_weights <- function(
  ncc,
  cohort,
  method = c("glm", "gam"),
  selection_formula = NULL,
  time,
  entry = NULL
) {
  call <- rlang::current_env()
  method <- match.arg(method)

  # Phase-1 cohort is required: the working model needs the full risk set at
  # every event time, which cannot be reconstructed from the NCC alone.
  if (is.null(cohort)) {
    rlang::abort(
      c(
        paste0(
          "Working-model weight `method = \"",
          method,
          "\"` requires the full Phase-1 cohort data."
        ),
        i = paste0(
          "Supply the original cohort as `cohort`. ",
          "It must contain the `time` column and any covariates in ",
          "`selection_formula`."
        )
      ),
      class = c("matchatr_missing_phase1", "matchatr_error"),
      call = call
    )
  }

  if (!is.data.frame(ncc)) {
    rlang::abort(
      "`ncc` must be a data.frame or data.table from `sample_ncc(incl_prob = TRUE)`.",
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  if (!is.data.frame(cohort)) {
    rlang::abort(
      "`cohort` must be a data.frame or data.table.",
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }

  # `time` is required to determine who is at risk at each event time.
  check_string(time, call = call)
  if (!time %in% names(cohort)) {
    rlang::abort(
      c(
        paste0(
          "Column `",
          time,
          "` not found in `cohort`."
        ),
        i = paste0(
          "Working-model weights need the cohort's follow-up time column to ",
          "reconstruct risk sets. Set `time` to its name."
        )
      ),
      class = c("matchatr_missing_phase1", "matchatr_error"),
      call = call
    )
  }

  if (!is.null(entry)) {
    check_string(entry, call = call)
    if (!entry %in% names(cohort)) {
      rlang::abort(
        paste0("`entry` column `", entry, "` not found in `cohort`."),
        class = c("matchatr_bad_input", "matchatr_error"),
        call = call
      )
    }
  }

  # `.cohort_row` is required: it maps each NCC row back to the cohort row that
  # generated it, which is how we identify controls in the augmented selection
  # dataset. Call sample_ncc(incl_prob = TRUE) to get this column.
  if (!".cohort_row" %in% names(ncc)) {
    rlang::abort(
      c(
        "Working-model weights require a `.cohort_row` column in `ncc`.",
        i = paste0(
          "Generate the NCC sample with ",
          "`sample_ncc(..., incl_prob = TRUE)` ",
          "to attach cohort row indices."
        )
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }

  required_cols <- c(".cohort_row", "case", "set", "risk_time")
  missing_cols <- setdiff(required_cols, names(ncc))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      c(
        paste0(
          "`ncc` is missing required column(s): ",
          paste0("`", missing_cols, "`", collapse = ", "),
          "."
        ),
        i = "Supply an NCC dataset from `sample_ncc(incl_prob = TRUE)`."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }

  if (method == "gam") {
    rlang::check_installed(
      "mgcv",
      reason = 'for working-model NCC weights with `method = "gam"`'
    )
  }

  cohort_df <- as.data.frame(cohort)
  ncc_df <- as.data.frame(ncc)
  n_cohort <- nrow(cohort_df)

  # Build the augmented dataset: one row per (eligible j, event k) pair.
  aug <- build_ncc_selection_dataset(
    ncc_df = ncc_df,
    cohort_df = cohort_df,
    time_col = time,
    entry_col = entry
  )

  # Construct the selection model formula.
  if (is.null(selection_formula)) {
    # Default: time-only working model (Borgan, Samuelsen & Aastveit 2003).
    fit_formula <- selected ~ risk_time
  } else {
    check_formula(selection_formula, call = call)
    fit_formula <- stats::update(
      stats::as.formula(selection_formula),
      selected ~ .
    )
  }

  # Fit the selection model and extract predicted probabilities.
  p_hat <- if (identical(method, "glm")) {
    glm_fit <- stats::glm(
      fit_formula,
      data = aug,
      family = stats::binomial()
    )
    stats::fitted(glm_fit)
  } else {
    gam_fit <- mgcv::gam(
      fit_formula,
      data = aug,
      family = stats::binomial()
    )
    stats::fitted(gam_fit)
  }

  # Compute per-cohort-subject inclusion probabilities via the product formula.
  pi_j <- working_model_inclusion_probs(
    cohort_rows = aug$.cohort_row,
    p_hat = p_hat,
    n_cohort = n_cohort
  )

  # Cases are always selected (probability 1); force their weight to 1 to
  # prevent 1/0 if the formula assigns them pi = 0 (they are excluded from
  # the selection pool at their own failure time).
  case_cohort_rows <- unique(ncc_df$.cohort_row[ncc_df$case == 1L])
  pi_j[case_cohort_rows] <- 1.0

  # Map inclusion probabilities to each NCC row via .cohort_row index.
  ncc_out <- ncc_df
  ncc_out$ipw_weight <- (1 / pi_j)[ncc_out$.cohort_row]
  # Enforce case weight = 1 even if a case appeared as a control in another set.
  ncc_out$ipw_weight[as.logical(ncc_out$case)] <- 1.0

  data.table::as.data.table(ncc_out)
}

#' Build the augmented selection dataset for working-model weight estimation
#'
#' For each event time t_k in the NCC sample, identifies the full risk set
#' R(t_k) (all cohort subjects at risk, excluding the case), marks which were
#' selected as controls, and stacks the result into a single data.frame with
#' cohort covariates merged. The resulting dataset is what the logistic /
#' GAM selection model is fitted on.
#'
#' @param ncc_df A data.frame with `.cohort_row`, `case`, `set`, `risk_time`.
#' @param cohort_df A data.frame (Phase-1 cohort) with the time column.
#' @param time_col Character scalar; name of the time column in `cohort_df`.
#' @param entry_col `NULL` or character scalar; name of the delayed-entry
#'   column in `cohort_df`.
#' @returns A data.frame with columns `.cohort_row`, `risk_time`, `selected`
#'   (1 = sampled control, 0 = eligible but not selected), plus all cohort
#'   columns available as covariates. Row-names are reset.
#' @family sampling
#' @noRd
build_ncc_selection_dataset <- function(
  ncc_df,
  cohort_df,
  time_col,
  entry_col
) {
  cohort_time <- cohort_df[[time_col]]
  cohort_entry <- if (!is.null(entry_col)) cohort_df[[entry_col]] else NULL

  sets <- sort(unique(ncc_df$set))
  aug_parts <- vector("list", length(sets))

  for (ki in seq_along(sets)) {
    sk <- sets[ki]
    set_mask <- ncc_df$set == sk
    # Failure time anchoring this risk set.
    risk_time_k <- ncc_df$risk_time[which(set_mask)[1L]]
    case_row_k <- ncc_df$.cohort_row[set_mask & ncc_df$case == 1L]
    ctrl_rows_k <- ncc_df$.cohort_row[set_mask & ncc_df$case == 0L]

    # Eligible pool: at risk at the case's failure time, excluding the case.
    # This mirrors the eligibility condition in eligible_controls().
    at_risk <- cohort_time >= risk_time_k
    if (!is.null(cohort_entry)) {
      at_risk <- at_risk & (cohort_entry < risk_time_k)
    }
    # The case is never in its own control pool.
    at_risk[case_row_k] <- FALSE
    elig_idx <- which(at_risk)

    if (length(elig_idx) == 0L) {
      next
    }

    aug_parts[[ki]] <- data.frame(
      .cohort_row = elig_idx,
      risk_time = risk_time_k,
      # Selected = 1 if this cohort subject was drawn as a control at t_k.
      selected = as.integer(elig_idx %in% ctrl_rows_k),
      stringsAsFactors = FALSE
    )
  }

  aug_parts <- aug_parts[!vapply(aug_parts, is.null, logical(1L))]
  aug <- do.call(rbind, aug_parts)
  rownames(aug) <- NULL

  # Merge all cohort columns so selection_formula can reference any of them.
  cohort_extra <- cohort_df[aug$.cohort_row, , drop = FALSE]
  rownames(cohort_extra) <- NULL
  cbind(aug, cohort_extra)
}

#' Compute inclusion probabilities from working-model fitted probabilities
#'
#' Applies the product formula to convert per-(subject, event) predicted
#' selection probabilities into per-subject inclusion probabilities:
#'
#'   π_j = 1 − ∏_{k: j ∈ R(t_k)} (1 − p̂_jk)
#'
#' Numerically stable via log-sum: log(1 − π_j) = Σ log(1 − p̂_jk), and
#' `π_j = −expm1(Σ log(1 − p̂_jk))`. Predicted probabilities are clamped below
#' 1 − ε to avoid log(0).
#'
#' @param cohort_rows Integer vector of cohort row indices aligned to `p_hat`.
#' @param p_hat Numeric vector of fitted probabilities P(selected | j, t_k).
#' @param n_cohort Integer total cohort size.
#' @returns Numeric vector of length `n_cohort` with inclusion probabilities
#'   π_j. Subjects never in any risk set (π_j = 0) are not in the NCC sample;
#'   the returned 0 is safe because they are excluded by index.
#' @family sampling
#' @noRd
working_model_inclusion_probs <- function(cohort_rows, p_hat, n_cohort) {
  # Clamp to avoid log(1 - 1) = log(0) when a predicted probability is exactly 1.
  p_clamped <- pmin(p_hat, 1 - .Machine$double.eps)
  log_contrib <- log1p(-p_clamped)

  # Accumulate log(1 - p_jk) per cohort subject using tapply.
  by_subject <- tapply(log_contrib, cohort_rows, sum)
  idx <- as.integer(names(by_subject))

  log_surv_fac <- numeric(n_cohort)
  log_surv_fac[idx] <- as.numeric(by_subject)

  # π_j = 1 - exp(log_surv_fac) = -expm1(log_surv_fac)
  -expm1(log_surv_fac)
}
