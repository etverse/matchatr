# sample_ncc() draws a nested case-control sample from a cohort by risk-set
# (incidence-density) control sampling: each event anchors a matched set (a
# sampled risk set) containing the case and m controls drawn from those at risk
# at the case's failure time. The result feeds matcha(design = nested_cc()),
# whose conditional partial likelihood reports the hazard ratio. Oracles:
# structural invariants of risk-set sampling (one case per set, controls at
# risk, no reuse within a set); Epi::ccwc as an external definition-of-risk-set
# cross-check; a cohort DGP with a known Cox log-HR (make_ncc_cohort) the sampled
# subsample must recover; and the classical (m+1)/m null efficiency.

# --- structural validity -------------------------------------------------

test_that("a risk-set sample has one case and m controls at risk per set", {
  co <- make_ncc_cohort(n = 1200L, seed = 51L)
  m <- 3L
  ncc <- withr::with_seed(11L, sample_ncc(co, time = "t", event = "d", m = m))

  expect_s3_class(ncc, "data.table")
  # One set per cohort event, exactly one case per set.
  expect_identical(length(unique(ncc$set)), sum(co$d))
  expect_true(all(tapply(ncc$case, ncc$set, sum) == 1L))
  # Every set has between 1 and m controls (a smaller set only at late times).
  set_sizes <- tabulate(ncc$set)
  expect_true(all(set_sizes >= 2L & set_sizes <= m + 1L))
  # Controls are genuinely at risk: their cohort exit time is at or after the
  # set's failure time.
  ctrl <- ncc[ncc$case == 0L, ]
  expect_true(all(ctrl$t >= ctrl$risk_time))
  # The case's risk_time is its own failure time.
  cases <- ncc[ncc$case == 1L, ]
  expect_true(all(cases$t == cases$risk_time))
  # No subject appears twice within one set (sampling without replacement).
  dup_within_set <- tapply(ncc$id, ncc$set, function(ids) {
    anyDuplicated(ids) > 0L
  })
  expect_false(any(dup_within_set))
})

test_that("sample_ncc does not mutate the input cohort", {
  co <- make_ncc_cohort(n = 300L, seed = 51L)
  before <- data.table::copy(co)
  withr::with_seed(1L, sample_ncc(co, time = "t", event = "d", m = 2L))
  expect_identical(co, before)
})

test_that("sample_ncc does not mutate a data.table cohort", {
  # The as.data.frame() copy is the mutation guard; a data.table input is the
  # case where in-place `:=` semantics would otherwise bite.
  co <- data.table::as.data.table(make_ncc_cohort(n = 300L, seed = 51L))
  before <- data.table::copy(co)
  withr::with_seed(1L, sample_ncc(co, time = "t", event = "d", m = 2L))
  expect_identical(co, before)
})

test_that("a control may serve before its own later event", {
  # The defining NCC structure: a subject sampled as a control can itself fail
  # later in the cohort. With a known early case and later events, at least one
  # sampled control must be a future cohort case (d == 1).
  co <- make_ncc_cohort(n = 800L, seed = 51L)
  ncc <- withr::with_seed(3L, sample_ncc(co, time = "t", event = "d", m = 4L))
  ctrl <- ncc[ncc$case == 0L, ]
  expect_true(any(ctrl$d == 1L))
})

# --- tie handling at equal event times -----------------------------------

test_that("a tied co-case is in the risk set (incidence-density semantics)", {
  # Two subjects fail at the same time t = 2. At that instant each is still at
  # risk for the other's event (the risk set is {t >= tc}), so the eligible pool
  # for one case includes the tied co-case. This pins the Langholz/Borgan
  # convention, which the continuous-time cohorts never exercise (and where the
  # behaviour diverges from Epi::ccwc's tied-failure grouping).
  coh <- data.frame(id = 1:5, t = c(2, 2, 5, 8, 9), d = c(1, 1, 0, 0, 0), x = 1:5)
  # Case in row 1 (t = 2): the tied co-case in row 2 (t = 2) is at risk (t >= 2).
  pool <- eligible_controls(1L, tvec = coh$t, entryvec = NULL, match_key = NULL)
  expect_true(2L %in% pool)
  # A full sample over both tied cases forms two sets (one per event), and each
  # set can draw the other tied case as a control.
  ncc <- withr::with_seed(1L, sample_ncc(coh, time = "t", event = "d", m = 3L))
  expect_identical(length(unique(ncc$set)), 2L)
})

# --- Epi::ccwc external oracle (definition of the risk set) ---------------

test_that("the risk-set pool agrees with Epi::ccwc", {
  skip_if_not_installed("Epi")
  co <- make_ncc_cohort(n = 800L, seed = 51L)
  tvec <- co$t

  # Direction 1: every control Epi SAMPLES (at the analysis m) must lie in our
  # eligible pool -- our pool is not too narrow.
  epi <- epi_ccwc_riskset(co, time = "t", event = "d", m = 3L)
  expect_identical(length(unique(epi$set)), sum(co$d))
  by_set <- split(epi, epi$set)
  for (s in by_set) {
    case_row <- s$row[s$case == 1L]
    expect_length(case_row, 1L)
    pool <- eligible_controls(case_row, tvec, entryvec = NULL, match_key = NULL)
    expect_true(all(s$row[s$case == 0L] %in% pool))
    # Epi's set failure time is the case's exit time.
    expect_equal(unique(s$risk_time), tvec[case_row])
  }

  # Direction 2: asking Epi for more controls than any risk set can hold makes it
  # return the FULL eligible set per case, so our pool must equal it EXACTLY --
  # our pool is not too wide either (a one-directional subset check would miss an
  # over-permissive `t > tc` boundary or a missing self-exclusion). Over-
  # requesting makes Epi warn "sets are incomplete" -- that is the intended
  # signal that it returned every available control, so it is suppressed here.
  epi_full <- suppressWarnings(
    epi_ccwc_riskset(co, time = "t", event = "d", m = nrow(co))
  )
  for (s in split(epi_full, epi_full$set)) {
    case_row <- s$row[s$case == 1L]
    pool <- eligible_controls(case_row, tvec, entryvec = NULL, match_key = NULL)
    expect_setequal(s$row[s$case == 0L], pool)
  }
})

# --- truth-based: the sampled NCC recovers the cohort Cox log-HR ----------

test_that("a sampled NCC recovers the cohort Cox log-HR (truth-based)", {
  co <- make_ncc_cohort(beta_x = log(2.2), beta_z = log(1.5))
  truth <- attr(co, "truth")
  ncc <- withr::with_seed(71L, sample_ncc(co, time = "t", event = "d", m = 4L))
  fit <- matcha(
    ncc,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    confounders = ~z,
    estimator = "clogit"
  )
  res <- contrast(fit)
  # SELF-SCALING band of 3.5 reported SEs (the truth-DGP convention): the
  # estimator's sampling SD is the reported SE, so a fixed absolute tolerance
  # would pass only by luck of the seed.
  expect_lt(
    abs(log(res$contrasts$estimate) - unname(truth["beta_x"])),
    3.5 * res$estimates$se
  )
})

test_that("a sampled NCC log-HR agrees with the full-cohort coxph beta", {
  co <- make_ncc_cohort(beta_x = log(2.2), beta_z = log(1.5))
  cox <- survival::coxph(survival::Surv(t, d) ~ x + z, data = co)
  b_cohort <- unname(stats::coef(cox)["x"])
  se_cohort <- sqrt(stats::vcov(cox)["x", "x"])

  ncc <- withr::with_seed(5L, sample_ncc(co, time = "t", event = "d", m = 4L))
  fit <- matcha(
    ncc,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    confounders = ~z,
    estimator = "clogit"
  )
  res <- contrast(fit)
  b_ncc <- log(res$contrasts$estimate)
  se_ncc <- res$estimates$se
  # The NCC subsample targets the cohort Cox beta (OR = HR); a combined-SE band
  # is conservative because the two estimates are positively correlated.
  expect_lt(abs(b_ncc - b_cohort), 3.5 * sqrt(se_ncc^2 + se_cohort^2))
})

# --- relative efficiency m/(m+1) at the null -----------------------------

test_that("a sample_ncc NCC has (m+1)/m efficiency at the null", {
  # At beta = 0 the NCC partial-likelihood information per case is m/(m+1) of the
  # full-cohort information (Goldstein & Langholz 1992), so Var_ncc / Var_cohort
  # is (m+1)/m. Monte-Carlo pin (single seed, large cohort); generous band.
  m <- 2L
  co0 <- make_ncc_cohort(n = 8000L, beta_x = 0, beta_z = 0, seed = 99L)
  cox0 <- survival::coxph(survival::Surv(t, d) ~ x, data = co0)
  se_cohort <- sqrt(stats::vcov(cox0)["x", "x"])

  ncc0 <- withr::with_seed(5L, sample_ncc(co0, time = "t", event = "d", m = m))
  fit <- matcha(
    ncc0,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    estimator = "clogit"
  )
  se_ncc <- contrast(fit)$estimates$se
  var_ratio <- (se_ncc / se_cohort)^2
  expect_gt(var_ratio, 1.30)
  expect_lt(var_ratio, 1.75)
})

# --- additional matching on population strata ----------------------------

test_that("additional matching confines each control to the case's stratum", {
  co <- make_ncc_cohort(n = 1000L, beta_x = log(2))
  co$s <- factor(ifelse(co$z > 0, "hi", "lo"), levels = c("lo", "hi"))
  ncc <- withr::with_seed(
    7L,
    sample_ncc(
      co,
      time = "t",
      event = "d",
      m = 2L,
      match = ~s
    )
  )
  # Within every set the matching variable is constant (case and its controls
  # share the stratum).
  by_set <- tapply(as.character(ncc$s), ncc$set, function(v) length(unique(v)))
  expect_true(all(by_set == 1L))
  # Still analysis-ready end to end.
  fit <- matcha(
    ncc,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    estimator = "clogit"
  )
  expect_true(is.finite(contrast(fit)$contrasts$estimate))
})

test_that("matching crosses several columns into one stratum exactly", {
  # Two match columns must be crossed jointly, not merged: every set is constant
  # on BOTH simultaneously. Guards against a string-key collision conflating two
  # distinct (s1, s2) tuples into one stratum.
  co <- make_ncc_cohort(n = 1500L, beta_x = log(2))
  co$s1 <- factor(ifelse(co$z > 0, "hi", "lo"))
  co$s2 <- factor(ifelse(co$id %% 2L == 0L, "even", "odd"))
  ncc <- withr::with_seed(
    9L,
    sample_ncc(co, time = "t", event = "d", m = 2L, match = ~ s1 + s2)
  )
  by_set <- tapply(seq_len(nrow(ncc)), ncc$set, function(ix) {
    length(unique(paste(ncc$s1[ix], ncc$s2[ix]))) == 1L
  })
  expect_true(all(by_set))
})

test_that("a missing value in a match column is rejected", {
  # A missing matching value leaves the stratum undefined; it must not silently
  # merge with the literal string "NA" or match anyone, so it is rejected up
  # front rather than producing a quietly mis-stratified sample.
  co <- make_ncc_cohort(n = 200L)
  co$s <- factor(ifelse(co$z > 0, "hi", "lo"))
  co$s[1] <- NA
  expect_error(
    sample_ncc(co, time = "t", event = "d", m = 2L, match = ~s),
    class = "matchatr_bad_input"
  )
})

# --- delayed entry / left truncation -------------------------------------

test_that("delayed entry excludes subjects not yet under observation", {
  # Subject 4 enters follow-up at t = 5, AFTER the only case fails at t = 2, so
  # it is not at risk for that event and must never be sampled as a control --
  # even though its exit time (10) is the largest in the cohort.
  coh <- data.frame(
    id = 1:6,
    entry = c(0, 0, 0, 5, 0, 0),
    t = c(2, 8, 9, 10, 7, 6),
    d = c(1, 0, 0, 0, 0, 0),
    x = c(1, 0, 1, 0, 1, 0)
  )
  ncc <- withr::with_seed(
    1L,
    sample_ncc(
      coh,
      time = "t",
      event = "d",
      m = 4L,
      entry = "entry"
    )
  )
  expect_false(4L %in% ncc$id)
  # Without left truncation the same subject IS eligible (its exit time qualifies).
  ncc_noentry <- withr::with_seed(
    1L,
    sample_ncc(
      coh,
      time = "t",
      event = "d",
      m = 4L
    )
  )
  expect_true(4L %in% ncc_noentry$id)
})

# --- fewer than m eligible controls --> a smaller set, not an error -------

test_that("a late failure time with fewer than m controls yields a smaller set", {
  # Cases at t = 1 and t = 3 (max time t = 5 is a non-case, so no empty set). At
  # tc = 3 only two subjects remain at risk, fewer than m = 3.
  coh <- data.frame(
    id = 1:5,
    t = c(1, 2, 3, 5, 4),
    d = c(1, 0, 1, 0, 0),
    x = c(1, 0, 1, 0, 1)
  )
  ncc <- withr::with_seed(1L, sample_ncc(coh, time = "t", event = "d", m = 3L))
  set_sizes <- tabulate(ncc$set)
  # No set exceeds m + 1, and the late set is strictly smaller.
  expect_true(all(set_sizes <= 4L))
  expect_true(any(set_sizes < 4L))
})

# --- rejections ----------------------------------------------------------

test_that("a case with no eligible control is rejected (empty risk set)", {
  # The latest time is the case's, so its risk set has no other member.
  bad <- data.frame(id = 1:3, t = c(1, 2, 3), d = c(0, 0, 1), x = c(1, 0, 1))
  expect_error(
    sample_ncc(bad, time = "t", event = "d", m = 2L),
    class = "matchatr_empty_risk_set"
  )
})

test_that("a missing time / event column is rejected", {
  co <- make_ncc_cohort(n = 200L)
  expect_error(
    sample_ncc(co, time = "nope", event = "d", m = 2L),
    class = "matchatr_bad_design"
  )
  expect_error(
    sample_ncc(co, time = "t", event = "nope", m = 2L),
    class = "matchatr_bad_design"
  )
})

test_that("a non-0/1 or event-free event column is rejected", {
  co <- make_ncc_cohort(n = 200L)
  # Continuous "event" is not an indicator.
  co_cont <- co
  co_cont$d <- co_cont$z
  expect_error(
    sample_ncc(co_cont, time = "t", event = "d", m = 2L),
    class = "matchatr_bad_outcome"
  )
  # No events to anchor any risk set.
  co_none <- co
  co_none$d <- 0L
  expect_error(
    sample_ncc(co_none, time = "t", event = "d", m = 2L),
    class = "matchatr_bad_outcome"
  )
})

test_that("a non-whole / sub-1 / NULL m is rejected", {
  co <- make_ncc_cohort(n = 200L)
  for (bad_m in list(0L, 1.5, -1L, NA_integer_, Inf, NULL)) {
    expect_error(
      sample_ncc(co, time = "t", event = "d", m = bad_m),
      class = "matchatr_bad_ratio"
    )
  }
})

test_that("a bad match specification is rejected", {
  co <- make_ncc_cohort(n = 200L)
  # A string is not a formula.
  expect_error(
    sample_ncc(co, time = "t", event = "d", m = 2L, match = "z"),
    class = "matchatr_bad_input"
  )
  # A formula naming an absent column.
  expect_error(
    sample_ncc(co, time = "t", event = "d", m = 2L, match = ~nope),
    class = "matchatr_bad_design"
  )
})

test_that("an output-name collision is rejected", {
  co <- make_ncc_cohort(n = 200L)
  co$set <- 1L
  expect_error(
    sample_ncc(co, time = "t", event = "d", m = 2L),
    class = "matchatr_bad_input"
  )
})

test_that("a non-data.frame cohort and non-numeric time are rejected", {
  expect_error(
    sample_ncc(1:10, time = "t", event = "d", m = 2L),
    class = "matchatr_bad_input"
  )
  co <- make_ncc_cohort(n = 200L)
  co$t <- as.character(co$t)
  expect_error(
    sample_ncc(co, time = "t", event = "d", m = 2L),
    class = "matchatr_bad_input"
  )
})

# --- message snapshots ---------------------------------------------------

test_that("the empty-risk-set and collision errors read clearly", {
  bad <- data.frame(id = 1:3, t = c(1, 2, 3), d = c(0, 0, 1), x = c(1, 0, 1))
  expect_snapshot(
    sample_ncc(bad, time = "t", event = "d", m = 2L),
    error = TRUE
  )
  clash <- data.frame(id = 1:4, t = c(1, 2, 3, 4), d = c(1, 0, 1, 0), case = 0L)
  expect_snapshot(
    sample_ncc(clash, time = "t", event = "d", m = 1L),
    error = TRUE
  )
})
