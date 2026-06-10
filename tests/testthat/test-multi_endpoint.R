# Multiple endpoints from one reused NCC control set (PHASE_7 Chunk 4). Tests
# the two reuse modes of the IPW weighted Cox:
#   (A) combined-event sampling: sample_ncc() on the union "any-failure" event
#       ascertains every endpoint's cases at once; each cause-specific endpoint
#       is then analysed directly via matcha(outcome = "<cause>", "ipw_cox").
#   (B) cohort-augmented reuse: sample_ncc() on a primary endpoint, then
#       reuse_ncc_endpoint() augments the sample with the secondary endpoint's
#       unsampled cohort cases before matcha(outcome = "<secondary>").
#
# Oracles:
#   - Mode (A): multipleNCC::wpl() reproduces the Samuelsen weighted partial
#     likelihood exactly for each endpoint (its $coefficients exposes the
#     endpoint coded 2 in samplestat). The Surv status is ignored by wpl; the
#     endpoint and weights come from samplestat. Exact agreement is expected.
#   - Mode (B): an independent reconstruction — KMprob() for the d1-only
#     inclusion probabilities + survival::coxph() for the weighted fit — matches
#     to machine precision, validating the augmentation and weight assembly.
#   - Truth-based: a competing-risks cohort with known cause-specific Cox log-HRs;
#     each reuse mode recovers the full-cohort cause-specific HR within a 3.5-SE
#     band.
# No Python oracle is added: as in Chunk 1, the weighted Cox robust sandwich has
# no convenient statsmodels equivalent and is validated exactly against
# survival::coxph and the canonical multipleNCC::wpl.

# ---- DGP: competing risks with known cause-specific Cox log-HRs -------------

# Two latent failure times with cause-specific hazards depending on the binary
# exposure and a continuous confounder; the realised cause is whichever fires
# first before administrative censoring at `tau`. The cause-specific Cox log-HR
# for the exposure is `beta1_x` (cause 1) and `beta2_x` (cause 2).
make_competing_cohort <- function(
  n = 4000L,
  beta1_x = log(2),
  beta2_x = log(0.5),
  seed = 20L
) {
  withr::with_seed(seed, {
    x <- stats::rbinom(n, 1L, 0.4)
    z <- stats::rnorm(n)
    r1 <- 0.05 * exp(beta1_x * x + 0.3 * z)
    r2 <- 0.05 * exp(beta2_x * x + 0.2 * z)
    t1 <- stats::rexp(n, r1)
    t2 <- stats::rexp(n, r2)
    tau <- 6
    tt <- pmin(t1, t2, tau)
    cause <- ifelse(tt >= tau, 0L, ifelse(t1 < t2, 1L, 2L))
    cohort <- data.frame(
      id = seq_len(n),
      t = tt,
      d_any = as.integer(cause != 0L),
      d1 = as.integer(cause == 1L),
      d2 = as.integer(cause == 2L),
      x = x,
      z = z
    )
    attr(cohort, "truth") <- c(beta1_x = beta1_x, beta2_x = beta2_x)
    cohort
  })
}

# samplestat for multipleNCC::wpl() with two endpoints: 0 = not sampled, 1 =
# sampled control, 2 = the analysed endpoint (reachable via $coefficients), 3 =
# the other endpoint (used as additional controls + extra event times for KM).
multi_samplestat <- function(cohort, ctrl_rows, analysed_col, other_col) {
  n <- nrow(cohort)
  ss <- rep(0L, n)
  ss[ctrl_rows] <- 1L
  ss[cohort[[other_col]] == 1L] <- 3L
  ss[cohort[[analysed_col]] == 1L] <- 2L
  ss
}

ipw_hr <- function(ncc, outcome, cohort_time = "t") {
  fit <- matcha(
    ncc,
    outcome = outcome,
    exposure = "x",
    design = nested_cc(strata = "set", time = cohort_time),
    confounders = ~z,
    estimator = "ipw_cox"
  )
  res <- contrast(fit)
  c(log_hr = res$estimates$estimate[1], se = res$estimates$se[1])
}

# ---- structural: reuse_ncc_endpoint() --------------------------------------

test_that("reuse_ncc_endpoint augments unsampled secondary cases with weight 1", {
  cohort <- make_competing_cohort(n = 1500L, seed = 3L)
  ncc <- withr::with_seed(
    7L,
    sample_ncc(cohort, "t", "d1", m = 3L, incl_prob = TRUE)
  )
  reused <- reuse_ncc_endpoint(ncc, cohort = cohort, time = "t", event = "d2")

  # Same columns as the input NCC, in the same order.
  expect_identical(names(reused), names(ncc))
  expect_s3_class(reused, "data.table")

  # Augmented rows = secondary cases not already in the sample.
  d2_rows <- which(cohort$d2 == 1L)
  in_sample <- unique(ncc$.cohort_row)
  expect_equal(nrow(reused) - nrow(ncc), length(setdiff(d2_rows, in_sample)))

  # Every secondary-endpoint subject carries weight 1 in the reused data.
  d2_mask <- reused$.cohort_row %in% d2_rows
  expect_equal(reused$ipw_weight[d2_mask], rep(1, sum(d2_mask)))

  # Every cohort secondary case is represented exactly once.
  expect_setequal(reused$.cohort_row[reused$d2 == 1L], d2_rows)
})

test_that("reuse_ncc_endpoint on a combined-event NCC needs no augmentation", {
  cohort <- make_competing_cohort(n = 1500L, seed = 4L)
  # Sampling on the union event ascertains every endpoint's cases already.
  ncc <- withr::with_seed(
    7L,
    sample_ncc(cohort, "t", "d_any", m = 3L, incl_prob = TRUE)
  )
  reused <- reuse_ncc_endpoint(ncc, cohort = cohort, time = "t", event = "d2")
  # No rows added; the input is returned unchanged (as a data.table).
  expect_equal(nrow(reused), nrow(ncc))
  expect_setequal(reused$.cohort_row, ncc$.cohort_row)
})

# ---- mode (A): combined-event reuse, exact multipleNCC::wpl oracle ---------

test_that("combined-event reuse matches multipleNCC::wpl exactly for each endpoint", {
  skip_if_not_installed("multipleNCC")
  cohort <- make_competing_cohort(n = 2500L, seed = 20L)
  ncc <- withr::with_seed(
    7L,
    sample_ncc(cohort, "t", "d_any", m = 3L, incl_prob = TRUE)
  )

  est_d1 <- ipw_hr(ncc, "d1")
  est_d2 <- ipw_hr(ncc, "d2")

  ctrl_rows <- unique(ncc$.cohort_row[ncc$case == 0L])
  ss1 <- multi_samplestat(
    cohort,
    ctrl_rows,
    analysed_col = "d1",
    other_col = "d2"
  )
  ss2 <- multi_samplestat(
    cohort,
    ctrl_rows,
    analysed_col = "d2",
    other_col = "d1"
  )
  w1 <- multipleNCC::wpl(
    survival::Surv(t, d1) ~ x + z,
    data = cohort,
    samplestat = ss1,
    m = 3L,
    weight.method = "KM",
    variance = "robust"
  )
  w2 <- multipleNCC::wpl(
    survival::Surv(t, d2) ~ x + z,
    data = cohort,
    samplestat = ss2,
    m = 3L,
    weight.method = "KM",
    variance = "robust"
  )

  # wpl ascertains every endpoint's cases at weight 1 and computes π over all
  # event times — identical to the combined-event NCC, so agreement is exact.
  expect_equal(
    unname(est_d1["log_hr"]),
    unname(w1$coefficients[1]),
    tolerance = 1e-6
  )
  expect_equal(unname(est_d1["se"]), sqrt(w1$var[1, 1]), tolerance = 1e-6)
  expect_equal(
    unname(est_d2["log_hr"]),
    unname(w2$coefficients[1]),
    tolerance = 1e-6
  )
  expect_equal(unname(est_d2["se"]), sqrt(w2$var[1, 1]), tolerance = 1e-6)
})

# ---- mode (A): truth recovery for both endpoints ---------------------------

test_that("combined-event reuse recovers both cause-specific cohort HRs", {
  cohort <- make_competing_cohort(n = 4000L, seed = 21L)
  b1 <- unname(stats::coef(survival::coxph(
    survival::Surv(t, d1) ~ x + z,
    cohort
  ))["x"])
  b2 <- unname(stats::coef(survival::coxph(
    survival::Surv(t, d2) ~ x + z,
    cohort
  ))["x"])

  ncc <- withr::with_seed(
    7L,
    sample_ncc(cohort, "t", "d_any", m = 3L, incl_prob = TRUE)
  )
  est_d1 <- ipw_hr(ncc, "d1")
  est_d2 <- ipw_hr(ncc, "d2")

  # The IPW weighted partial likelihood is unbiased for the cause-specific Cox
  # log-HR; 3.5-SE band (the truth-DGP convention).
  expect_lt(abs(est_d1["log_hr"] - b1), 3.5 * est_d1["se"])
  expect_lt(abs(est_d2["log_hr"] - b2), 3.5 * est_d2["se"])
})

# ---- mode (B): cohort-augmented reuse, exact independent reconstruction ----

test_that("cohort-augmented reuse matches an independent KMprob + coxph fit", {
  skip_if_not_installed("multipleNCC")
  cohort <- make_competing_cohort(n = 2500L, seed = 20L)
  ncc <- withr::with_seed(
    7L,
    sample_ncc(cohort, "t", "d1", m = 3L, incl_prob = TRUE)
  )
  reused <- reuse_ncc_endpoint(ncc, cohort = cohort, time = "t", event = "d2")
  est <- ipw_hr(reused, "d2")

  # Independent reconstruction: d1-only inclusion probabilities from KMprob, the
  # analysis sample assembled by hand (all d1/d2 cases at weight 1, sampled
  # controls at 1/π), fit with survival::coxph(ties = "breslow").
  n <- nrow(cohort)
  ctrl_rows <- unique(ncc$.cohort_row[ncc$case == 0L])
  ss <- rep(0L, n)
  ss[ctrl_rows] <- 1L
  ss[cohort$d1 == 1L] <- 2L
  pi_km <- multipleNCC::KMprob(survtime = cohort$t, samplestat = ss, m = 3)
  keep <- sort(unique(c(
    which(cohort$d1 == 1L),
    which(cohort$d2 == 1L),
    ctrl_rows
  )))
  w <- ifelse(cohort$d1[keep] == 1L | cohort$d2[keep] == 1L, 1, 1 / pi_km[keep])
  samp <- cohort[keep, ]
  ora <- survival::coxph(
    survival::Surv(t, d2) ~ x + z,
    data = samp,
    weights = w,
    robust = TRUE,
    ties = "breslow"
  )

  expect_equal(
    unname(est["log_hr"]),
    unname(stats::coef(ora)["x"]),
    tolerance = 1e-6
  )
  expect_equal(
    unname(est["se"]),
    sqrt(stats::vcov(ora)["x", "x"]),
    tolerance = 1e-6
  )
})

# ---- mode (B): truth recovery ----------------------------------------------

test_that("cohort-augmented reuse recovers the secondary cohort HR", {
  cohort <- make_competing_cohort(n = 4000L, seed = 21L)
  b2 <- unname(stats::coef(survival::coxph(
    survival::Surv(t, d2) ~ x + z,
    cohort
  ))["x"])

  ncc <- withr::with_seed(
    7L,
    sample_ncc(cohort, "t", "d1", m = 3L, incl_prob = TRUE)
  )
  reused <- reuse_ncc_endpoint(ncc, cohort = cohort, time = "t", event = "d2")
  est <- ipw_hr(reused, "d2")

  expect_lt(abs(est["log_hr"] - b2), 3.5 * est["se"])
})

# ---- competing-endpoint cases keep weight 1 (the ncc_ipw_analysis_data fix) -

test_that("a primary case sampled as a control keeps weight 1 in the reused fit", {
  # Construct a case that is sampled as a control before its own (later) event,
  # so deduplication can retain its control row. The generalised analysis-sample
  # builder must still give it weight 1 (it is ascertained), not 1/π.
  cohort <- make_competing_cohort(n = 1200L, seed = 8L)
  ncc <- withr::with_seed(
    7L,
    sample_ncc(cohort, "t", "d_any", m = 3L, incl_prob = TRUE)
  )

  # Identify a d1 case that also appears as a control somewhere in the sample.
  d1_case_rows <- unique(ncc$.cohort_row[ncc$case == 1L & ncc$d1 == 1L])
  ctrl_rows <- unique(ncc$.cohort_row[ncc$case == 0L])
  reused_as_control <- intersect(d1_case_rows, ctrl_rows)
  skip_if(
    length(reused_as_control) == 0L,
    "no primary case reused as a control"
  )

  fit <- matcha(
    ncc,
    outcome = "d2",
    exposure = "x",
    design = nested_cc(strata = "set", time = "t"),
    confounders = ~z,
    estimator = "ipw_cox"
  )
  analysis <- ncc_ipw_analysis_data(fit)
  subj <- reused_as_control[1L]
  # This subject is a d1 case (competing event for the d2 analysis), ascertained
  # with probability 1, so its weight in the deduplicated analysis sample is 1.
  expect_equal(analysis$ipw_weight[analysis$.cohort_row == subj], 1)
})

# ---- rejections -------------------------------------------------------------

test_that("reuse_ncc_endpoint aborts with matchatr_missing_phase1 when cohort = NULL", {
  cohort <- make_competing_cohort(n = 400L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d1", m = 2L, incl_prob = TRUE)
  )
  expect_error(
    reuse_ncc_endpoint(ncc, cohort = NULL, time = "t", event = "d2"),
    class = "matchatr_missing_phase1"
  )
})

test_that("reuse_ncc_endpoint aborts when time column is absent from cohort", {
  cohort <- make_competing_cohort(n = 400L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d1", m = 2L, incl_prob = TRUE)
  )
  bad <- cohort
  names(bad)[names(bad) == "t"] <- "time_renamed"
  expect_error(
    reuse_ncc_endpoint(ncc, cohort = bad, time = "t", event = "d2"),
    class = "matchatr_missing_phase1"
  )
})

test_that("reuse_ncc_endpoint aborts when .cohort_row is absent from ncc", {
  cohort <- make_competing_cohort(n = 400L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d1", m = 2L, incl_prob = TRUE)
  )
  ncc$.cohort_row <- NULL
  expect_error(
    reuse_ncc_endpoint(ncc, cohort = cohort, time = "t", event = "d2"),
    class = "matchatr_bad_input"
  )
})

test_that("reuse_ncc_endpoint aborts when ncc lacks a required bookkeeping column", {
  cohort <- make_competing_cohort(n = 400L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d1", m = 2L, incl_prob = TRUE)
  )
  ncc$set <- NULL
  expect_error(
    reuse_ncc_endpoint(ncc, cohort = cohort, time = "t", event = "d2"),
    class = "matchatr_bad_input"
  )
})

test_that("reuse_ncc_endpoint aborts when the secondary event is absent from cohort", {
  cohort <- make_competing_cohort(n = 400L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d1", m = 2L, incl_prob = TRUE)
  )
  expect_error(
    reuse_ncc_endpoint(
      ncc,
      cohort = cohort,
      time = "t",
      event = "not_a_column"
    ),
    class = "matchatr_bad_input"
  )
})

test_that("reuse_ncc_endpoint aborts when the secondary endpoint has no cases", {
  cohort <- make_competing_cohort(n = 400L, seed = 1L)
  cohort$d_none <- 0L
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d1", m = 2L, incl_prob = TRUE)
  )
  expect_error(
    reuse_ncc_endpoint(ncc, cohort = cohort, time = "t", event = "d_none"),
    class = "matchatr_bad_outcome"
  )
})
