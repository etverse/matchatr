#' Draw a nested case-control sample from a cohort by risk-set sampling
#'
#' @description
#' Generates a nested case-control (NCC) dataset from a full cohort by risk-set
#' (incidence-density) sampling: at each event time the failing subject (the
#' case) is matched to `m` controls drawn at random from the subjects still at
#' risk at that instant. Each case and its sampled controls form one matched set
#' (a sampled risk set), which is the stratum the conditional partial likelihood
#' conditions on. The result feeds straight into
#' `matcha(design = nested_cc(strata = "set", time = "risk_time"))`, whose
#' conditional logistic fit reports the hazard ratio (OR = HR exactly under
#' risk-set sampling; Prentice & Breslow 1978).
#'
#' @details
#' The risk set at a failure time `tc` is the set of subjects under observation
#' then: those who have entered follow-up and not yet left, i.e.
#' `entry < tc <= time` (with `entry` taken as the origin when not supplied, so
#' the condition reduces to `time >= tc`). The case is excluded from its own
#' control pool, and `min(m, n_eligible)` controls are sampled without
#' replacement; a late failure time with fewer than `m` eligible controls yields
#' a smaller set rather than an error, mirroring real NCC sampling. A subject
#' sampled as a control may itself fail later in the cohort — it serves as a
#' control before its own event — so the per-set `case` indicator is distinct
#' from the cohort-wide `event` column.
#'
#' Additional matching (`match = ~ s1 + s2`) restricts each case's control pool
#' to subjects sharing the case's values on the named population-stratum
#' variables (e.g. sex, birth cohort). Sampling uses the ambient random-number
#' stream, so wrap the call in [withr::with_seed()] or precede it with
#' [set.seed()] for reproducibility.
#'
#' The sampler is implemented natively (base R + data.table) so it is always
#' available and deterministically seedable; [Epi::ccwc()] is an equivalent
#' external implementation used as a cross-check in the package's tests.
#'
#' A case whose risk set contains no eligible control carries no
#' conditional-likelihood information and signals a sampling failure (a
#' misspecified time origin/scale, an entry/exit mismatch, or over-fine `match`
#' strata), so it aborts with `matchatr_empty_risk_set` rather than silently
#' producing a singleton set.
#'
#' @param cohort A data.frame or data.table with one row per subject. Copied,
#'   never mutated.
#' @param time A single character string naming the exit / event-time column (the
#'   time scale on which risk sets are formed). Must be numeric.
#' @param event A single character string naming the cohort event indicator
#'   (logical, a two-level factor, or numeric 0/1); at least one event must
#'   occur.
#' @param m A single whole number >= 1, the number of controls sampled per case
#'   (default 1).
#' @param match `NULL` (no additional matching) or a one-sided formula naming
#'   population-stratum column(s) controls must share with the case
#'   (e.g. `~ sex + birth_cohort`).
#' @param entry `NULL` (everyone enters at the time origin) or a single character
#'   string naming a delayed-entry / left-truncation column. Must be numeric when
#'   supplied; a subject is at risk at `tc` only if `entry < tc`.
#'
#' @param incl_prob Logical; if `TRUE`, compute Samuelsen (1997) Kaplan-Meier
#'   inclusion probabilities for each sampled subject and append two columns to
#'   the output: `.cohort_row` (integer row index of the subject in `cohort`)
#'   and `ipw_weight` (inverse inclusion probability: 1 for cases, 1/pi_j for
#'   sampled controls). The weights are needed by the IPW nested case-control
#'   analysis (`estimator = "ipw_cox"`). The inclusion probability for control
#'   j is pi_j = 1 - prod over event times where j was eligible of
#'   (1 - m_i / n_elig_i), where m_i is the controls sampled and n_elig_i the
#'   eligible pool size at event i (Samuelsen 1997, Biometrika). Computation
#'   is O(n x K) where n is the cohort size and K is the number of events.
#'   Default `FALSE`.
#'
#' @returns A `data.table` with one row per sampled subject: the selected rows of
#'   `cohort` (all original columns) plus `set` (integer matched-set id), `case`
#'   (per-set 0/1 indicator, 1 for the case), and `risk_time` (the set's failure
#'   time). When `incl_prob = TRUE`, also `.cohort_row` (integer, 1-indexed row
#'   index in the original `cohort`) and `ipw_weight` (numeric, 1/pi_j).
#'   Aborts with `matchatr_empty_risk_set` when any case has no eligible control.
#'
#' @examples
#' # A small cohort with event times; draw 2 controls per case.
#' cohort <- data.frame(
#'   id = 1:8,
#'   t  = c(2, 5, 1, 8, 3, 9, 4, 7),
#'   d  = c(1, 0, 1, 0, 1, 0, 0, 0),
#'   x  = c(1, 0, 1, 0, 0, 1, 0, 1)
#' )
#' set.seed(1)
#' ncc <- sample_ncc(cohort, time = "t", event = "d", m = 2)
#' ncc
#'
#' # Analyse it: each sampled risk set is a stratum -> hazard ratio.
#' fit <- matcha(ncc, outcome = "case", exposure = "x",
#'               design = nested_cc(strata = "set", time = "risk_time"),
#'               estimator = "clogit")
#' contrast(fit)
#'
#' # With Samuelsen KM inclusion probabilities for IPW analysis.
#' set.seed(1)
#' ncc_ipw <- sample_ncc(cohort, time = "t", event = "d", m = 2, incl_prob = TRUE)
#' ncc_ipw[, c("id", "case", "set", "ipw_weight", ".cohort_row")]
#'
#' @family sampling
#' @seealso [nested_cc()], [matcha()], [Epi::ccwc()]
#' @export
sample_ncc <- function(
  cohort,
  time,
  event,
  m = 1,
  match = NULL,
  entry = NULL,
  incl_prob = FALSE
) {
  call <- rlang::current_env()
  if (!is.data.frame(cohort)) {
    rlang::abort(
      "`cohort` must be a data.frame or data.table.",
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  check_string(time, call = call)
  check_string(event, call = call)
  check_sample_m(m, call = call)
  if (!is.null(entry)) {
    check_string(entry, call = call)
  }
  match_vars <- character(0)
  if (!is.null(match)) {
    check_formula(match, call = call)
    match_vars <- all.vars(match)
    if (length(match_vars) == 0L) {
      rlang::abort(
        "`match` must name at least one column (e.g. `~ sex`).",
        class = c("matchatr_bad_input", "matchatr_error"),
        call = call
      )
    }
  }
  check_unique_colnames(cohort, call = call)
  # Every referenced column must exist before sampling reads it.
  check_cols_exist(cohort, time, arg = "time", call = call)
  check_cols_exist(cohort, event, arg = "event", call = call)
  if (!is.null(entry)) {
    check_cols_exist(cohort, entry, arg = "entry", call = call)
  }
  if (length(match_vars) > 0L) {
    check_cols_exist(cohort, match_vars, arg = "match", call = call)
  }
  # The output columns are appended; a cohort that already carries them would
  # collide (and `[[` would silently resolve to the wrong one downstream).
  added <- c("set", "case", "risk_time")
  clash <- intersect(added, names(cohort))
  if (length(clash) > 0L) {
    rlang::abort(
      c(
        paste0(
          "`cohort` already has column(s) ",
          paste0("`", clash, "`", collapse = ", "),
          ", which `sample_ncc()` appends to the output."
        ),
        i = "Rename them so the sampled set id / case indicator / risk time are unambiguous."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }

  # A base data.frame copy drives the row subsetting: data.table's `[`
  # row-selection is gated on the calling namespace being "data.table-aware"
  # (cedta), which is brittle under load_all, so the sampler reads columns with
  # `[[`, assembles the result by base row indexing, and converts to data.table
  # only at the return. Adding the output columns to the fresh subset never
  # touches `cohort`, so the user's data is not mutated.
  df <- as.data.frame(cohort)
  tvec <- df[[time]]
  if (!is.numeric(tvec)) {
    rlang::abort(
      paste0("`time` column `", time, "` must be numeric."),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  entryvec <- NULL
  if (!is.null(entry)) {
    entryvec <- df[[entry]]
    if (!is.numeric(entryvec)) {
      rlang::abort(
        paste0("`entry` column `", entry, "` must be numeric."),
        class = c("matchatr_bad_input", "matchatr_error"),
        call = call
      )
    }
  }
  evec <- resolve_event_indicator(df, event, call = call)

  n <- nrow(df)
  rid <- seq_len(n)
  # Cases anchor the risk sets; a case needs a finite failure time to define one.
  case_rows <- rid[which(evec == 1L)]
  if (anyNA(tvec[case_rows]) || any(!is.finite(tvec[case_rows]))) {
    rlang::abort(
      paste0(
        "Some case rows have a missing or non-finite `",
        time,
        "`, so their risk set cannot be formed."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  # Order cases by failure time (ties broken by row) so the set index is stable
  # and independent of the input row order.
  case_rows <- case_rows[order(tvec[case_rows], case_rows)]

  # A single per-row stratum code (the match columns crossed) lets the
  # eligibility test be one vector comparison per case. A missing matching value
  # leaves the stratum undefined -- it must not silently merge with the literal
  # string "NA" (which a pasted character key would do), nor match anyone -- so
  # NA in a match column is rejected up front.
  match_key <- NULL
  if (length(match_vars) > 0L) {
    na_cols <- match_vars[vapply(
      match_vars,
      function(v) anyNA(df[[v]]),
      logical(1)
    )]
    if (length(na_cols) > 0L) {
      rlang::abort(
        c(
          paste0(
            "`match` column(s) ",
            paste0("`", na_cols, "`", collapse = ", "),
            " contain missing values, so the matching stratum is undefined."
          ),
          i = "Drop or impute the missing strata before sampling."
        ),
        class = c("matchatr_bad_input", "matchatr_error"),
        call = call
      )
    }
    # interaction() crosses the match columns into one factor; its integer codes
    # compare exactly, with no separator a value could collide on (unlike a
    # pasted string key). unname() + do.call spreads the columns as the factors
    # interaction() expects.
    match_key <- as.integer(do.call(
      interaction,
      c(unname(df[match_vars]), list(drop = TRUE))
    ))
  }

  # Pass 1 -- eligibility only, no random draws. Computing every case's eligible
  # pool first means the empty-risk-set abort fires BEFORE any RNG is consumed,
  # so a failed call leaves the random stream untouched.
  elig_list <- vector("list", length(case_rows))
  for (k in seq_along(case_rows)) {
    elig_list[[k]] <- eligible_controls(
      case_row = case_rows[k],
      tvec = tvec,
      entryvec = entryvec,
      match_key = match_key
    )
  }
  empty <- case_rows[lengths(elig_list) == 0L]
  if (length(empty) > 0L) {
    reject_empty_risk_set(empty, tvec, call = call)
  }

  # Pass 2 -- sample min(m, n_eligible) controls per case, accumulating the
  # selected rows and their per-set labels for a single base subset at the end. A
  # late failure time with fewer than m eligible controls keeps all of them (a
  # smaller set), which is correct, not an error.
  all_rows <- integer(0)
  set_id <- integer(0)
  case_ind <- integer(0)
  risk_t <- numeric(0)
  for (k in seq_along(case_rows)) {
    ci <- case_rows[k]
    elig <- elig_list[[k]]
    take <- if (length(elig) > m) elig[sample.int(length(elig), m)] else elig
    members <- c(ci, take)
    all_rows <- c(all_rows, members)
    set_id <- c(set_id, rep.int(k, length(members)))
    case_ind <- c(case_ind, as.integer(members == ci))
    risk_t <- c(risk_t, rep.int(tvec[ci], length(members)))
  }

  out <- df[all_rows, , drop = FALSE]
  out[["set"]] <- set_id
  out[["case"]] <- case_ind
  out[["risk_time"]] <- risk_t
  rownames(out) <- NULL

  if (isTRUE(incl_prob)) {
    added_ipw <- c(".cohort_row", "ipw_weight")
    clash_ipw <- intersect(added_ipw, names(df))
    if (length(clash_ipw) > 0L) {
      rlang::abort(
        c(
          paste0(
            "`cohort` already has column(s) ",
            paste0("`", clash_ipw, "`", collapse = ", "),
            ", which `sample_ncc(incl_prob = TRUE)` appends to the output."
          ),
          i = "Rename them before calling `sample_ncc()`."
        ),
        class = c("matchatr_bad_input", "matchatr_error"),
        call = call
      )
    }
    # Samuelsen (1997) KM inclusion probability for each cohort subject:
    # π_j = 1 − ∏_{i: j ∈ elig(t_i)} (1 − m_i / |elig(t_i)|)
    # where the product runs over ALL event times where j was eligible, and
    # m_i = min(m, |elig(t_i)|) is the controls actually sampled at event i.
    ipw <- samuelsen_km_weights(
      n = n,
      case_rows = case_rows,
      elig_list = elig_list,
      m_requested = m
    )
    out[[".cohort_row"]] <- all_rows
    # ipw already contains 1/π_j with cases forced to 1; no further override.
    out[["ipw_weight"]] <- ipw[all_rows]
  }

  data.table::as.data.table(out)
}

#' Compute the eligible control rows for one case's risk set
#'
#' The risk set at the case's failure time `tc` is everyone under observation
#' then -- `entry < tc <= time` -- restricted (when matching) to the case's
#' population stratum, with the case itself removed. `which()` drops the `NA`s a
#' missing time / entry comparison produces, so a subject with unknown timing is
#' simply not eligible.
#'
#' @param case_row Integer row index of the case in the cohort.
#' @param tvec Numeric vector of exit / event times.
#' @param entryvec `NULL` (no delayed entry) or a numeric vector of entry times.
#' @param match_key `NULL` (no matching) or an integer vector of per-row crossed
#'   stratum codes (from `interaction()`), compared exactly to the case's code.
#' @returns An integer vector of eligible control row indices (possibly empty).
#' @family sampling
#' @seealso [sample_ncc()]
#' @noRd
eligible_controls <- function(case_row, tvec, entryvec, match_key) {
  tc <- tvec[case_row]
  at_risk <- tvec >= tc
  if (!is.null(entryvec)) {
    # Delayed entry: a subject is only at risk after entering follow-up.
    at_risk <- at_risk & (entryvec < tc)
  }
  if (!is.null(match_key)) {
    at_risk <- at_risk & (match_key == match_key[case_row])
  }
  at_risk[case_row] <- FALSE
  which(at_risk)
}

#' Samuelsen KM inclusion probabilities for a nested case-control sample
#'
#' Computes the Kaplan-Meier inclusion probability pi_j for each subject in the
#' full cohort using the Samuelsen (1997) formula:
#'
#'   pi_j = 1 - prod_{i: j in elig(t_i)} (1 - m_i / n_elig_i)
#'
#' where the product runs over ALL event times t_i where subject j was eligible
#' (in the risk set), m_i = min(m_requested, |elig(t_i)|) is the number of
#' controls actually sampled, and |elig(t_i)| is the eligible pool size at t_i.
#' Cases (in case_rows) are not eligible for their own risk set; they may be
#' eligible at earlier events. The returned vector has length n (cohort
#' dimension); non-sampled subjects with π_j = 0 indicate that they were never
#' eligible as a control (no event time placed them in any risk set).
#'
#' @param n Integer cohort size.
#' @param case_rows Integer vector of case row indices in the cohort.
#' @param elig_list List with one element per event (aligned with case_rows);
#'   each element is an integer vector of eligible control row indices.
#' @param m_requested Integer; the nominal number of controls per case.
#' @returns Numeric vector of length `n` with IPW weights (1/π_j). Cohort
#'   cases (in `case_rows`) get weight 1; sampled controls get 1/π_j >= 1;
#'   never-eligible subjects get weight 0 (they cannot be in the NCC sample).
#' @family sampling
#' @seealso [sample_ncc()]
#' @noRd
samuelsen_km_weights <- function(n, case_rows, elig_list, m_requested) {
  n_elig <- lengths(elig_list)
  m_actual <- pmin(m_requested, n_elig)

  # Accumulate log(1 - m_i / |elig(t_i)|) for each cohort subject across all
  # event times where they appear in the eligible pool.
  log_surv_fac <- numeric(n) # starts at 0 (product = 1 -> π = 0)
  for (k in seq_along(case_rows)) {
    ne <- n_elig[k]
    if (ne == 0L) {
      next
    }
    # log(1 - m_actual / n_elig): when m_actual == n_elig the factor is 0 (log = -Inf),
    # meaning any subject in this pool was certainly sampled (contributes π -> 1).
    log_fac <- log1p(-m_actual[k] / ne)
    log_surv_fac[elig_list[[k]]] <- log_surv_fac[elig_list[[k]]] + log_fac
  }
  # π_j = 1 - exp(sum of log-factors); -expm1(x) = 1 - exp(x) is numerically
  # stable near zero. A subject eligible at every event where all eligible
  # controls were sampled gets log_surv_fac = -Inf -> π = 1.
  pi_j <- -expm1(log_surv_fac)

  # Cohort cases are always included in the NCC sample (as the failing subject
  # of their risk set), so their IPW weight is 1 regardless of what the control-
  # sampling probability formula gives. This also avoids 1/0 for cases whose
  # pi from the formula is 0 (never eligible as a control before their event).
  pi_j[case_rows] <- 1.0

  # Return inverse probabilities (weights for the Cox model). Non-sampled
  # subjects with pi = 0 get Inf but are never in all_rows, so they never
  # appear in the NCC output (1/0 is safe to return as-is).
  1 / pi_j
}

#' Draw a counter-matched nested case-control sample from a cohort
#'
#' @description
#' Generates a counter-matched nested case-control (NCC) dataset from a full
#' cohort: at each event time the case is matched to `m` controls drawn
#' exclusively from the *opposite* surrogate stratum, concentrating study
#' resources in subjects whose surrogate exposure differs from the case. The
#' analysis weights (log-sampling-weights) for the Langholz-Borgan (1995)
#' weighted partial likelihood are appended as a `log_w` column: the case
#' represents its entire surrogate stratum (log-weight = log(n_same + 1)) and
#' each control represents the opposite stratum divided by the controls drawn
#' (log-weight = log(n_other / m)). The result feeds straight into
#' `matcha(design = counter_matched(strata = "set", time = "risk_time",
#' weights = "log_w"))`, whose `survival::coxph` fit reports the hazard ratio.
#'
#' @details
#' Controls are drawn without replacement from subjects at risk at the case's
#' failure time (`entry < tc <= time`) in the *opposite* surrogate stratum.
#' When fewer than `m` eligible controls exist (late failure times or narrow
#' strata), the smaller available set is returned rather than an error. A case
#' with no eligible control in the opposite stratum — which means the entire
#' at-risk population shares the case's surrogate value — aborts with
#' `matchatr_empty_risk_set`, signalling a sampling-design failure.
#'
#' Additional matching (`match = ~ s1 + s2`) restricts each risk set to
#' subjects sharing the case's population-stratum values on the named
#' variables. The eligibility pool is first restricted by the match stratum,
#' and the surrogate split is then applied within that pool.
#'
#' The surrogate column must be binary: logical, numeric 0/1, or a two-level
#' factor. Missing values abort with `matchatr_bad_input`.
#'
#' The log-weights are the key difference from unweighted `sample_ncc()`: with
#' pure counter-matching the unweighted clogit is biased because the controls
#' were drawn non-uniformly (opposite-stratum only), so the weighted partial
#' likelihood (coxph + offset) must be used for analysis.
#'
#' @param cohort A data.frame or data.table with one row per subject.
#' @param time A single character string naming the exit / event-time column.
#'   Must be numeric.
#' @param event A single character string naming the event indicator (logical,
#'   two-level factor, or numeric 0/1); at least one event must occur.
#' @param surrogate A single character string naming the binary surrogate
#'   exposure column. Controls are drawn from subjects whose surrogate value
#'   differs from the case's. Must be logical, numeric 0/1, or a two-level
#'   factor; no missing values.
#' @param m A single whole number >= 1, the number of controls sampled per
#'   case from the opposite surrogate stratum (default 1).
#' @param match `NULL` (no additional matching) or a one-sided formula naming
#'   population-stratum column(s) controls must share with the case
#'   (e.g. `~ sex + birth_cohort`).
#' @param entry `NULL` (everyone enters at the time origin) or a single
#'   character string naming a delayed-entry / left-truncation column.
#'
#' @returns A `data.table` with one row per sampled subject: the selected rows
#'   of `cohort` (all original columns) plus `set` (integer matched-set id),
#'   `case` (per-set 0/1 indicator), `risk_time` (the set's failure time), and
#'   `log_w` (the log-sampling-weight for the weighted partial likelihood).
#'   Aborts with `matchatr_empty_risk_set` when any case has no eligible
#'   control in the opposite surrogate stratum.
#'
#' @examples
#' set.seed(1)
#' tt <- rexp(200, 0.1)
#' cohort <- data.frame(
#'   id  = 1:200,
#'   t   = pmin(tt, 15),
#'   d   = as.integer(tt <= 15),
#'   x   = rbinom(200, 1, 0.4),
#'   z   = rbinom(200, 1, 0.5)   # binary surrogate
#' )
#' ncc_cm <- sample_ncc_counter_matched(cohort, time = "t", event = "d",
#'                                      surrogate = "z", m = 1L)
#' # Analyse: weighted partial likelihood identifies the HR.
#' fit <- matcha(ncc_cm, outcome = "case", exposure = "x",
#'               design = counter_matched(strata = "set", time = "risk_time",
#'                                        weights = "log_w"),
#'               estimator = "weighted_cox")
#' contrast(fit)
#'
#' @family sampling
#' @seealso [counter_matched()], [matcha()], [sample_ncc()]
#' @export
sample_ncc_counter_matched <- function(
  cohort,
  time,
  event,
  surrogate,
  m = 1L,
  match = NULL,
  entry = NULL
) {
  call <- rlang::current_env()
  if (!is.data.frame(cohort)) {
    rlang::abort(
      "`cohort` must be a data.frame or data.table.",
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  check_string(time, call = call)
  check_string(event, call = call)
  check_string(surrogate, call = call)
  check_sample_m(m, call = call)
  if (!is.null(entry)) {
    check_string(entry, call = call)
  }
  match_vars <- character(0)
  if (!is.null(match)) {
    check_formula(match, call = call)
    match_vars <- all.vars(match)
    if (length(match_vars) == 0L) {
      rlang::abort(
        "`match` must name at least one column (e.g. `~ sex`).",
        class = c("matchatr_bad_input", "matchatr_error"),
        call = call
      )
    }
  }
  check_unique_colnames(cohort, call = call)
  check_cols_exist(cohort, time, arg = "time", call = call)
  check_cols_exist(cohort, event, arg = "event", call = call)
  check_cols_exist(cohort, surrogate, arg = "surrogate", call = call)
  if (!is.null(entry)) {
    check_cols_exist(cohort, entry, arg = "entry", call = call)
  }
  if (length(match_vars) > 0L) {
    check_cols_exist(cohort, match_vars, arg = "match", call = call)
  }
  added <- c("set", "case", "risk_time", "log_w")
  clash <- intersect(added, names(cohort))
  if (length(clash) > 0L) {
    rlang::abort(
      c(
        paste0(
          "`cohort` already has column(s) ",
          paste0("`", clash, "`", collapse = ", "),
          ", which `sample_ncc_counter_matched()` appends to the output."
        ),
        i = "Rename them so the sampled set id / case indicator / risk time / log-weights are unambiguous."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }

  df <- as.data.frame(cohort)
  tvec <- df[[time]]
  if (!is.numeric(tvec)) {
    rlang::abort(
      paste0("`time` column `", time, "` must be numeric."),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  entryvec <- NULL
  if (!is.null(entry)) {
    entryvec <- df[[entry]]
    if (!is.numeric(entryvec)) {
      rlang::abort(
        paste0("`entry` column `", entry, "` must be numeric."),
        class = c("matchatr_bad_input", "matchatr_error"),
        call = call
      )
    }
  }

  # The surrogate must be binary; coerce to integer 0/1 and validate.
  svec <- resolve_surrogate(df, surrogate, call = call)
  evec <- resolve_event_indicator(df, event, call = call)

  n <- nrow(df)
  rid <- seq_len(n)
  case_rows <- rid[which(evec == 1L)]
  if (anyNA(tvec[case_rows]) || any(!is.finite(tvec[case_rows]))) {
    rlang::abort(
      paste0(
        "Some case rows have a missing or non-finite `",
        time,
        "`, so their risk set cannot be formed."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  case_rows <- case_rows[order(tvec[case_rows], case_rows)]

  match_key <- NULL
  if (length(match_vars) > 0L) {
    na_cols <- match_vars[vapply(
      match_vars,
      function(v) anyNA(df[[v]]),
      logical(1)
    )]
    if (length(na_cols) > 0L) {
      rlang::abort(
        c(
          paste0(
            "`match` column(s) ",
            paste0("`", na_cols, "`", collapse = ", "),
            " contain missing values, so the matching stratum is undefined."
          ),
          i = "Drop or impute the missing strata before sampling."
        ),
        class = c("matchatr_bad_input", "matchatr_error"),
        call = call
      )
    }
    match_key <- as.integer(do.call(
      interaction,
      c(unname(df[match_vars]), list(drop = TRUE))
    ))
  }

  # Pass 1 — eligibility check without random draws. Identify the opposite-
  # stratum eligible pool for each case; abort before any RNG is consumed if
  # any pool is empty.
  elig_opp_list <- vector("list", length(case_rows))
  elig_same_list <- vector("list", length(case_rows))
  for (k in seq_along(case_rows)) {
    ci <- case_rows[k]
    zc <- svec[ci]
    # All at-risk subjects within the population stratum (time / entry / match).
    pool <- eligible_controls(ci, tvec, entryvec, match_key)
    elig_opp_list[[k]] <- pool[svec[pool] != zc]
    elig_same_list[[k]] <- pool[svec[pool] == zc]
  }
  empty <- case_rows[lengths(elig_opp_list) == 0L]
  if (length(empty) > 0L) {
    reject_empty_risk_set(empty, tvec, call = call)
  }

  # Pass 2 — sample min(m, n_opp) controls from the opposite surrogate stratum
  # for each case, compute log-weights, and accumulate output rows.
  all_rows <- integer(0)
  set_id <- integer(0)
  case_ind <- integer(0)
  risk_t <- numeric(0)
  log_w_vec <- numeric(0)

  for (k in seq_along(case_rows)) {
    ci <- case_rows[k]
    opp <- elig_opp_list[[k]]
    n_opp <- length(opp)
    n_same <- length(elig_same_list[[k]])
    # Draw min(m, n_opp) controls from the opposite stratum.
    m_take <- min(m, n_opp)
    ctrl_rows <- if (n_opp > m) opp[sample.int(n_opp, m)] else opp

    members <- c(ci, ctrl_rows)
    all_rows <- c(all_rows, members)
    set_id <- c(set_id, rep.int(k, length(members)))
    case_ind <- c(case_ind, as.integer(members == ci))
    risk_t <- c(risk_t, rep.int(tvec[ci], length(members)))

    # Langholz & Borgan (1995) sampling weights:
    #   case: represents all n_same + 1 at-risk in its surrogate stratum
    #   each control: one of m_take sampled from n_opp opposite-stratum subjects
    log_w_case <- log(n_same + 1L)
    log_w_ctrl <- log(n_opp / m_take)
    log_w_vec <- c(log_w_vec, log_w_case, rep.int(log_w_ctrl, m_take))
  }

  out <- df[all_rows, , drop = FALSE]
  out[["set"]] <- set_id
  out[["case"]] <- case_ind
  out[["risk_time"]] <- risk_t
  out[["log_w"]] <- log_w_vec
  rownames(out) <- NULL
  data.table::as.data.table(out)
}

#' Resolve the surrogate column to a binary integer vector
#'
#' The counter-matching surrogate must be binary (two distinct non-missing
#' values): logical, numeric 0/1, or a two-level factor. Missing values abort.
#' Returns a 0/1 integer vector aligned to the cohort rows.
#'
#' @param data data.frame copy of the cohort.
#' @param surrogate Character scalar column name.
#' @param call Caller environment for the error.
#' @returns An integer vector of 0/1 codes.
#' @family sampling
#' @noRd
resolve_surrogate <- function(data, surrogate, call = rlang::caller_env()) {
  sv <- data[[surrogate]]
  if (anyNA(sv)) {
    rlang::abort(
      paste0(
        "`surrogate` column `",
        surrogate,
        "` contains missing values. ",
        "Counter-matching requires a fully observed binary surrogate."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  if (is.logical(sv)) {
    return(as.integer(sv))
  }
  if (is.numeric(sv)) {
    uv <- sort(unique(sv))
    if (!identical(uv, c(0, 1)) && !identical(uv, 0:1)) {
      rlang::abort(
        paste0(
          "`surrogate` column `",
          surrogate,
          "` must have exactly two distinct values (0 and 1); ",
          "found: ",
          paste(uv, collapse = ", "),
          "."
        ),
        class = c("matchatr_bad_input", "matchatr_error"),
        call = call
      )
    }
    return(as.integer(sv))
  }
  if (is.factor(sv)) {
    lv <- levels(droplevels(sv))
    if (length(lv) != 2L) {
      rlang::abort(
        paste0(
          "`surrogate` column `",
          surrogate,
          "` must be a two-level factor for counter-matching; ",
          "found ",
          length(lv),
          " level(s)."
        ),
        class = c("matchatr_bad_input", "matchatr_error"),
        call = call
      )
    }
    return(as.integer(sv) - 1L)
  }
  rlang::abort(
    paste0(
      "`surrogate` column `",
      surrogate,
      "` must be logical, numeric 0/1, or a two-level factor."
    ),
    class = c("matchatr_bad_input", "matchatr_error"),
    call = call
  )
}
