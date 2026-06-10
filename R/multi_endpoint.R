#' Reuse a nested case-control control set for a second endpoint
#'
#' @description
#' Prepares a nested case-control (NCC) sample drawn for one (primary) endpoint
#' so it can be reused to fit a weighted Cox model for a *different* (secondary)
#' endpoint, sharing the single sampled control set. The sampled controls keep
#' their primary Samuelsen inclusion weights 1/π_j (the inclusion probability is
#' a property of the sampling, not of the endpoint analysed), and the secondary
#' endpoint's cases that were not already drawn into the sample are augmented
#' from the full Phase-1 cohort with weight 1 (every case of an analysed endpoint
#' is ascertained with probability 1). The result feeds straight into
#' `matcha(design = nested_cc(...), estimator = "ipw_cox")` with the secondary
#' endpoint as `outcome`.
#'
#' @details
#' Nested case-control studies often record several endpoints in the same cohort.
#' Classical (matched) analysis (Phase 5) ties each control to its case's failure
#' time, so a control set sampled for one endpoint cannot serve another. The IPW
#' reformulation breaks the matching: the union of cases and controls is treated
#' as a biased cohort subsample weighted by inverse inclusion probability, which
#' makes control reuse across endpoints legitimate (Samuelsen 1997; Saarela,
#' Kulathinal, Arjas & Läärä 2008; Støer & Samuelsen 2012).
#'
#' The reuse semantics are:
#'
#'   * **Controls** keep their primary inclusion weight 1/π_j. π_j is the
#'     probability the subject was drawn as a control during the *primary*
#'     risk-set sampling; it does not depend on which endpoint is later analysed,
#'     so reusing the primary weight is the theoretically correct
#'     inverse-sampling-probability weight.
#'   * **Secondary-endpoint cases** are ascertained with probability 1, hence
#'     weight 1. Those already present in the sample (drawn as a primary case or
#'     as a control) keep their row; those absent are augmented from the cohort,
#'     each as its own singleton risk set (a case with no sampled controls).
#'   * **Primary-endpoint cases** are competing events for the secondary
#'     analysis. They were ascertained by the primary sampling, so they also keep
#'     weight 1; the `ipw_cox` engine enforces this when it builds the
#'     deduplicated analysis sample.
#'
#' When the primary sample already contains every secondary-endpoint case — for
#' example when the NCC was drawn on the union "any-failure" event so all
#' endpoints' cases were ascertained at once — no augmentation is needed and the
#' input is returned unchanged. In that shared-control case the secondary
#' analysis can also be run directly with
#' `matcha(ncc, outcome = "<secondary>", estimator = "ipw_cox")`.
#'
#' This function does not recompute inclusion probabilities; the secondary
#' endpoint reuses the primary π_j attached by `sample_ncc(incl_prob = TRUE)`.
#' The robust Lin-Wei sandwich variance the `ipw_cox` engine reports treats the
#' weights as fixed (the standard, slightly conservative treatment).
#'
#' @param ncc A `data.table` or `data.frame` from [sample_ncc()] with
#'   `incl_prob = TRUE`. Must contain `.cohort_row`, `ipw_weight`, `set`, `case`,
#'   and `risk_time`, plus the cohort columns (carried by `sample_ncc()`) and a
#'   binary column for the secondary `event`.
#' @param cohort A `data.frame` or `data.table` — the **full Phase-1 cohort**
#'   from which `ncc` was drawn. Used to recover the secondary endpoint's cases
#'   that were not sampled. `NULL` is rejected with `matchatr_missing_phase1`.
#' @param time A single character string naming the cohort follow-up / event-time
#'   column (the same column passed to `nested_cc(time = ...)`). Augmented cases
#'   take their `risk_time` from it. Must be present in `cohort`.
#' @param event A single character string naming the secondary endpoint's binary
#'   indicator column in `cohort` (logical, two-level factor, or numeric 0/1)
#'   with at least one case. The same column must be present in `ncc` (it is, as
#'   a carried cohort column) and is passed as `outcome` to [matcha()].
#'
#' @returns A `data.table` with the same columns as `ncc`: the original NCC rows
#'   followed by one augmented row per previously unsampled secondary-endpoint
#'   case (weight 1, `case = 1`, its own `set` id beyond the sampled sets). The
#'   secondary cases' `ipw_weight` is set to 1 throughout. Ready for
#'   `matcha(design = nested_cc(strata = "set", time = time),
#'   estimator = "ipw_cox")` with `outcome = event`.
#'
#' @examples
#' set.seed(1)
#' t1 <- rexp(600, 0.05)
#' t2 <- rexp(600, 0.05)
#' tau <- 6
#' tt <- pmin(t1, t2, tau)
#' cause <- ifelse(tt >= tau, 0L, ifelse(t1 < t2, 1L, 2L))
#' cohort <- data.frame(
#'   id = 1:600,
#'   t  = tt,
#'   d1 = as.integer(cause == 1L),   # primary endpoint
#'   d2 = as.integer(cause == 2L),   # secondary endpoint
#'   x  = rbinom(600, 1, 0.4)
#' )
#' ncc <- sample_ncc(cohort, time = "t", event = "d1", m = 3, incl_prob = TRUE)
#'
#' # Reuse the one control set for the secondary endpoint
#' ncc2 <- reuse_ncc_endpoint(ncc, cohort = cohort, time = "t", event = "d2")
#' fit2 <- matcha(ncc2, outcome = "d2", exposure = "x",
#'                design = nested_cc(strata = "set", time = "t"),
#'                estimator = "ipw_cox")
#' contrast(fit2)
#'
#' @family sampling
#' @seealso [sample_ncc()], [compute_ncc_weights()], [matcha()], [nested_cc()]
#' @export
reuse_ncc_endpoint <- function(ncc, cohort, time, event) {
  call <- rlang::current_env()

  # The Phase-1 cohort is required: a secondary-endpoint case that was not drawn
  # into the primary sample exists only in the cohort, and its event must be
  # ascertained for the reuse analysis to be unbiased.
  if (is.null(cohort)) {
    rlang::abort(
      c(
        "Reusing a control set for a second endpoint requires the full Phase-1 cohort.",
        i = paste0(
          "Supply the original cohort as `cohort` so the secondary endpoint's ",
          "cases that were not sampled can be ascertained."
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
  check_string(time, call = call)
  check_string(event, call = call)

  # `.cohort_row` maps each NCC row back to its cohort subject; without it the
  # already-sampled secondary cases cannot be told apart from the augmented ones.
  if (!".cohort_row" %in% names(ncc)) {
    rlang::abort(
      c(
        "Reusing a control set requires a `.cohort_row` column in `ncc`.",
        i = "Generate the NCC sample with `sample_ncc(..., incl_prob = TRUE)`."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  required_cols <- c(".cohort_row", "ipw_weight", "set", "case", "risk_time")
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

  # `time` anchors the augmented cases on the cohort time scale.
  if (!time %in% names(cohort)) {
    rlang::abort(
      c(
        paste0("Column `", time, "` not found in `cohort`."),
        i = "Set `time` to the cohort follow-up / event-time column name."
      ),
      class = c("matchatr_missing_phase1", "matchatr_error"),
      call = call
    )
  }
  if (!event %in% names(cohort)) {
    rlang::abort(
      c(
        paste0(
          "Secondary endpoint column `",
          event,
          "` not found in `cohort`."
        ),
        i = "Set `event` to the binary column of the endpoint to reuse the control set for."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }

  cohort_df <- as.data.frame(cohort)
  ncc_df <- as.data.frame(ncc)

  # Resolve the secondary endpoint to a 0/1 indicator; reusing for an endpoint
  # with no cases is meaningless (`resolve_event_indicator` enforces >= 1 case).
  d2 <- resolve_event_indicator(cohort_df, event, call = call)
  d2_rows <- which(d2 == 1L)

  # Secondary cases already represented in the sample keep their rows; the rest
  # are augmented. Their inclusion weight is 1 (ascertained), set on every row
  # of the subject so a secondary case drawn as a control is not left at 1/π_j.
  in_sample <- unique(ncc_df[[".cohort_row"]])
  aug_rows <- setdiff(d2_rows, in_sample)
  ncc_df[["ipw_weight"]][ncc_df[[".cohort_row"]] %in% d2_rows] <- 1.0

  if (length(aug_rows) == 0L) {
    # Shared-control case: every secondary case was already ascertained (e.g. the
    # NCC was drawn on the union event), so there is nothing to augment.
    return(data.table::as.data.table(ncc_df))
  }

  # The NCC carries the cohort columns plus this bookkeeping; the augmented rows
  # need the same column set, pulled from the cohort.
  book <- c("set", "case", "risk_time", ".cohort_row", "ipw_weight")
  cohort_cols <- setdiff(names(ncc_df), book)
  missing_in_cohort <- setdiff(cohort_cols, names(cohort_df))
  if (length(missing_in_cohort) > 0L) {
    rlang::abort(
      c(
        paste0(
          "`ncc` carries column(s) absent from `cohort`: ",
          paste0("`", missing_in_cohort, "`", collapse = ", "),
          "."
        ),
        i = "Augmented secondary cases are pulled from `cohort`; its columns must cover the NCC's cohort columns."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }

  # Each augmented case is its own singleton risk set (a case with no sampled
  # control), beyond the existing set ids. The `ipw_cox` engine treats it as an
  # ordinary weight-1 observation; the `set` id is bookkeeping it does not read.
  aug <- cohort_df[aug_rows, cohort_cols, drop = FALSE]
  max_set <- as.integer(max(ncc_df[["set"]], na.rm = TRUE))
  aug[["set"]] <- max_set + seq_along(aug_rows)
  aug[["case"]] <- 1L
  aug[["risk_time"]] <- cohort_df[[time]][aug_rows]
  aug[[".cohort_row"]] <- aug_rows
  aug[["ipw_weight"]] <- 1.0

  aug <- aug[, names(ncc_df), drop = FALSE]
  out <- rbind(ncc_df, aug)
  rownames(out) <- NULL
  data.table::as.data.table(out)
}
