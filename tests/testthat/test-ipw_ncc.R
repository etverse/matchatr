# IPW for nested case-control (PHASE_7 Chunk 1). Tests cover:
#   - sample_ncc(incl_prob = TRUE): structural checks on ipw_weight / .cohort_row
#   - Samuelsen KM weight correctness vs multipleNCC::KMprob oracle
#   - fit_ipw_cox() + contrast(): HR recovers full-cohort coxph (truth-based)
#   - multipleNCC::wpl() oracle: exact agreement on coefficients and SE
#   - Rejection paths: missing ipw_weight / .cohort_row columns
#
# Oracle: multipleNCC::wpl() with weight.method = "KM" and variance = "robust"
# computes the same Samuelsen weighted partial likelihood as matchatr's
# fit_ipw_cox(). Exact agreement is expected because both implement the same
# weighted score equations.

# ---- helpers ----------------------------------------------------------------

# A cohort with a known Cox log-HR for the binary exposure.
make_ipw_cohort <- function(
  n = 2000L,
  beta_x = log(2),
  beta_z = 0.4,
  base_rate = 0.08,
  tau = 5,
  seed = 77L
) {
  withr::with_seed(seed, {
    exposure <- stats::rbinom(n, 1L, 0.4)
    confounder <- stats::rnorm(n)
    rate <- base_rate * exp(beta_x * exposure + beta_z * confounder)
    tt <- stats::rexp(n, rate)
    d <- as.integer(tt <= tau)
    t_obs <- pmin(tt, tau)
    cohort <- data.frame(
      id = seq_len(n),
      t = t_obs,
      d = d,
      exposure = exposure,
      confounder = confounder
    )
    attr(cohort, "truth") <- c(beta_x = beta_x, beta_z = beta_z)
    cohort
  })
}

# samplestat vector for multipleNCC::wpl(): 0 = not sampled, 1 = sampled
# control, 2 = event (case). When a subject is both a case and appeared as a
# control in an earlier set, the case status (2) takes precedence.
build_samplestat <- function(cohort, ncc) {
  n <- nrow(cohort)
  sstat <- rep(0L, n)
  ctrl_rows <- unique(ncc$.cohort_row[ncc$case == 0L])
  sstat[ctrl_rows] <- 1L
  sstat[cohort$d == 1L] <- 2L # cases override controls
  sstat
}

# ---- structural tests for sample_ncc(incl_prob = TRUE) ---------------------

test_that("sample_ncc(incl_prob = TRUE) adds .cohort_row and ipw_weight", {
  cohort <- make_ipw_cohort()
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  expect_true(".cohort_row" %in% names(ncc))
  expect_true("ipw_weight" %in% names(ncc))
})

test_that("sample_ncc(incl_prob = TRUE): cases have weight 1", {
  cohort <- make_ipw_cohort()
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  case_weights <- ncc$ipw_weight[ncc$case == 1L]
  expect_equal(case_weights, rep(1, sum(ncc$case)))
})

test_that("sample_ncc(incl_prob = TRUE): control weights are >= 1", {
  cohort <- make_ipw_cohort()
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  ctrl_weights <- ncc$ipw_weight[ncc$case == 0L]
  expect_true(all(ctrl_weights >= 1 - .Machine$double.eps))
})

test_that("sample_ncc(incl_prob = TRUE): .cohort_row indexes original cohort", {
  cohort <- make_ipw_cohort()
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  expect_true(all(ncc$.cohort_row >= 1L & ncc$.cohort_row <= nrow(cohort)))
  # Cohort data matches: each NCC row's id equals cohort[.cohort_row]$id
  expect_equal(ncc$id, cohort$id[ncc$.cohort_row])
})

test_that("sample_ncc(incl_prob = FALSE) does not add weight columns", {
  cohort <- make_ipw_cohort()
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = FALSE)
  )
  expect_false(".cohort_row" %in% names(ncc))
  expect_false("ipw_weight" %in% names(ncc))
})

test_that("sample_ncc: clash on .cohort_row column aborts", {
  cohort <- make_ipw_cohort()
  cohort$.cohort_row <- 99L
  expect_error(
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE),
    class = "matchatr_bad_input"
  )
})

# ---- KM weight oracle: compare to multipleNCC::KMprob ----------------------

test_that("Samuelsen KM weights match multipleNCC::KMprob (inverse)", {
  skip_if_not_installed("multipleNCC")
  cohort <- make_ipw_cohort(n = 800L, seed = 13L)
  ncc <- withr::with_seed(
    5L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )

  # My weights (1/π_j) for sampled controls
  ctrl_mask <- ncc$case == 0L
  my_w <- ncc$ipw_weight[ctrl_mask]
  ctrl_rows <- ncc$.cohort_row[ctrl_mask]

  # KMprob returns π_j (probabilities); 1/π_j should equal my weights
  sstat <- build_samplestat(cohort, ncc)
  km_pi <- multipleNCC::KMprob(
    survtime = cohort$t,
    samplestat = sstat,
    m = 3
  )
  oracle_w <- 1 / km_pi[ctrl_rows]

  # Exact agreement is expected: both implement the same KM product formula.
  expect_equal(my_w, oracle_w, tolerance = 1e-8)
})

# ---- truth-based: IPW Cox recovers full-cohort Cox log-HR ------------------

test_that("IPW Cox HR recovers full-cohort Cox HR (truth-based)", {
  cohort <- make_ipw_cohort(n = 3000L, beta_x = log(2), seed = 41L)
  truth <- attr(cohort, "truth")

  # Full cohort Cox (the estimand the NCC IPW analysis targets)
  full_cox <- survival::coxph(
    survival::Surv(t, d) ~ exposure + confounder,
    data = cohort
  )
  target_log_hr <- unname(stats::coef(full_cox)["exposure"])

  ncc <- withr::with_seed(
    3L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    outcome = "d",
    exposure = "exposure",
    design = nested_cc(strata = "set", time = "t"),
    confounders = ~confounder,
    estimator = "ipw_cox"
  )
  res <- contrast(fit)
  log_hr <- res$estimates$estimate[1]
  se_hr <- res$estimates$se[1]

  # The IPW weighted partial likelihood is asymptotically unbiased for the
  # Cox hazard ratio. Use a 3.5-SE band (the truth-DGP convention).
  expect_lt(abs(log_hr - target_log_hr), 3.5 * se_hr)
})

# ---- multipleNCC::wpl oracle: exact HR and SE agreement --------------------

test_that("IPW Cox HR and SE match multipleNCC::wpl (KM weights, robust)", {
  skip_if_not_installed("multipleNCC")
  cohort <- make_ipw_cohort(n = 2000L, seed = 99L)

  ncc <- withr::with_seed(
    7L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    outcome = "d",
    exposure = "exposure",
    design = nested_cc(strata = "set", time = "t"),
    confounders = ~confounder,
    estimator = "ipw_cox"
  )
  res <- contrast(fit)
  log_hr <- res$estimates$estimate[1]
  se_hr <- res$estimates$se[1]

  # multipleNCC::wpl() is the canonical implementation of the Samuelsen
  # weighted partial likelihood. Exact agreement (to floating-point precision)
  # is expected because both optimise the same score equations.
  sstat <- build_samplestat(cohort, ncc)
  oracle <- multipleNCC::wpl(
    survival::Surv(t, d) ~ exposure + confounder,
    data = cohort,
    samplestat = sstat,
    m = 3L,
    weight.method = "KM",
    variance = "robust"
  )
  oracle_log_hr <- unname(oracle$coefficients[1])
  oracle_se <- sqrt(oracle$var[1, 1])

  expect_equal(log_hr, oracle_log_hr, tolerance = 1e-6)
  expect_equal(se_hr, oracle_se, tolerance = 1e-6)
})

# ---- contrast: type and ci_method rejection --------------------------------

test_that("contrast(type = 'or') is rejected for ipw_cox", {
  cohort <- make_ipw_cohort(n = 400L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    "d",
    "exposure",
    nested_cc(strata = "set", time = "t"),
    estimator = "ipw_cox"
  )
  expect_error(
    contrast(fit, type = "or"),
    class = "matchatr_unidentified_estimand"
  )
})

test_that("contrast(type = 'difference') is rejected for ipw_cox", {
  cohort <- make_ipw_cohort(n = 400L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    "d",
    "exposure",
    nested_cc(strata = "set", time = "t"),
    estimator = "ipw_cox"
  )
  expect_error(
    contrast(fit, type = "difference"),
    class = "matchatr_unidentified_estimand"
  )
})

test_that("ci_method = 'bootstrap' is rejected for ipw_cox", {
  cohort <- make_ipw_cohort(n = 400L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    "d",
    "exposure",
    nested_cc(strata = "set", time = "t"),
    estimator = "ipw_cox"
  )
  expect_error(
    contrast(fit, type = "hr", ci_method = "bootstrap"),
    class = "matchatr_unsupported_variance"
  )
})

# ---- rejection: missing ipw_weight / .cohort_row ---------------------------

test_that("matcha with ipw_cox and no ipw_weight column aborts", {
  cohort <- make_ipw_cohort(n = 400L, seed = 1L)
  ncc <- withr::with_seed(1L, sample_ncc(cohort, "t", "d", m = 2L)) # no incl_prob
  expect_error(
    matcha(
      ncc,
      "d",
      "exposure",
      nested_cc(strata = "set", time = "t"),
      estimator = "ipw_cox"
    ),
    class = "matchatr_missing_ipw_weights"
  )
})

test_that("matcha with ipw_cox and no .cohort_row column aborts", {
  cohort <- make_ipw_cohort(n = 400L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE)
  )
  ncc$.cohort_row <- NULL
  expect_error(
    matcha(
      ncc,
      "d",
      "exposure",
      nested_cc(strata = "set", time = "t"),
      estimator = "ipw_cox"
    ),
    class = "matchatr_missing_ipw_weights"
  )
})

test_that("estimator = 'ipw_cox' is rejected for non-nested-cc designs", {
  df <- make_cc_data()
  expect_error(
    matcha(df, "case", "x", matched_cc(strata = "set"), estimator = "ipw_cox"),
    class = "matchatr_bad_estimator"
  )
})

# ---- dispatch: ipw_cox engine wired correctly ------------------------------

test_that("matcha with ipw_cox returns correct engine", {
  cohort <- make_ipw_cohort(n = 400L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    "d",
    "exposure",
    nested_cc(strata = "set", time = "t"),
    estimator = "ipw_cox"
  )
  expect_identical(fit$engine, "ipw_cox")
  expect_s3_class(fit$model, "coxph")
})

test_that("default contrast type for ipw_cox is hr", {
  cohort <- make_ipw_cohort(n = 400L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    "d",
    "exposure",
    nested_cc(strata = "set", time = "t"),
    estimator = "ipw_cox"
  )
  res <- contrast(fit) # no type argument; defaults to "hr"
  expect_identical(res$type, "hr")
})

# ============================================================================
# PHASE_7 Chunk 2: working-model (GLM/GAM) weights and rejection paths
# ============================================================================

# ---- structural tests for compute_ncc_weights() ----------------------------

test_that("compute_ncc_weights GLM: cases have weight 1, controls >= 1", {
  cohort <- make_ipw_cohort(n = 800L, seed = 31L)
  ncc <- withr::with_seed(
    5L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  ncc_glm <- compute_ncc_weights(
    ncc,
    cohort = cohort,
    method = "glm",
    time = "t"
  )
  expect_equal(
    ncc_glm$ipw_weight[ncc_glm$case == 1L],
    rep(1, sum(ncc_glm$case))
  )
  expect_true(all(
    ncc_glm$ipw_weight[ncc_glm$case == 0L] >= 1 - .Machine$double.eps
  ))
})

test_that("compute_ncc_weights GLM: weights are finite", {
  cohort <- make_ipw_cohort(n = 800L, seed = 32L)
  ncc <- withr::with_seed(
    6L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  ncc_glm <- compute_ncc_weights(
    ncc,
    cohort = cohort,
    method = "glm",
    time = "t"
  )
  expect_true(all(is.finite(ncc_glm$ipw_weight)))
})

test_that("compute_ncc_weights GAM: cases have weight 1, controls >= 1", {
  skip_if_not_installed("mgcv")
  cohort <- make_ipw_cohort(n = 800L, seed = 33L)
  ncc <- withr::with_seed(
    7L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  ncc_gam <- compute_ncc_weights(
    ncc,
    cohort = cohort,
    method = "gam",
    selection_formula = ~ s(risk_time),
    time = "t"
  )
  expect_equal(
    ncc_gam$ipw_weight[ncc_gam$case == 1L],
    rep(1, sum(ncc_gam$case))
  )
  expect_true(all(
    ncc_gam$ipw_weight[ncc_gam$case == 0L] >= 1 - .Machine$double.eps
  ))
})

test_that("compute_ncc_weights: .cohort_row is preserved unchanged", {
  cohort <- make_ipw_cohort(n = 600L, seed = 34L)
  ncc <- withr::with_seed(
    8L,
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE)
  )
  ncc_glm <- compute_ncc_weights(
    ncc,
    cohort = cohort,
    method = "glm",
    time = "t"
  )
  expect_equal(ncc_glm$.cohort_row, ncc$.cohort_row)
})

test_that("compute_ncc_weights: custom selection_formula is respected", {
  cohort <- make_ipw_cohort(n = 800L, seed = 35L)
  ncc <- withr::with_seed(
    9L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  ncc_glm <- compute_ncc_weights(
    ncc,
    cohort = cohort,
    method = "glm",
    selection_formula = ~ risk_time + confounder,
    time = "t"
  )
  expect_true(all(is.finite(ncc_glm$ipw_weight)))
  expect_true(all(
    ncc_glm$ipw_weight[ncc_glm$case == 0L] >= 1 - .Machine$double.eps
  ))
})

# ---- truth-based: GLM-weighted IPW Cox recovers the cohort Cox HR ----------

test_that("GLM-weighted IPW Cox HR recovers full-cohort Cox HR (truth-based)", {
  cohort <- make_ipw_cohort(n = 3000L, beta_x = log(2), seed = 42L)
  full_cox <- survival::coxph(
    survival::Surv(t, d) ~ exposure + confounder,
    data = cohort
  )
  target_log_hr <- unname(stats::coef(full_cox)["exposure"])

  ncc <- withr::with_seed(
    11L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  ncc_glm <- compute_ncc_weights(
    ncc,
    cohort = cohort,
    method = "glm",
    time = "t"
  )
  fit <- matcha(
    ncc_glm,
    outcome = "d",
    exposure = "exposure",
    design = nested_cc(strata = "set", time = "t"),
    confounders = ~confounder,
    estimator = "ipw_cox"
  )
  res <- contrast(fit)
  log_hr <- res$estimates$estimate[1]
  se_hr <- res$estimates$se[1]

  expect_lt(abs(log_hr - target_log_hr), 3.5 * se_hr)
})

# ---- cross-implementation comparison with multipleNCC oracle ---------------

test_that("GLM-weighted HR agrees with multipleNCC::wpl(glm) within 2e-2", {
  skip_if_not_installed("multipleNCC")
  cohort <- make_ipw_cohort(n = 2000L, seed = 99L)

  ncc <- withr::with_seed(
    7L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  ncc_glm <- compute_ncc_weights(
    ncc,
    cohort = cohort,
    method = "glm",
    time = "t"
  )
  fit <- matcha(
    ncc_glm,
    outcome = "d",
    exposure = "exposure",
    design = nested_cc(strata = "set", time = "t"),
    confounders = ~confounder,
    estimator = "ipw_cox"
  )
  res <- contrast(fit)
  our_log_hr <- res$estimates$estimate[1]

  # Build samplestat vector for multipleNCC
  sstat <- build_samplestat(cohort, ncc)
  oracle <- multipleNCC::wpl(
    survival::Surv(t, d) ~ exposure + confounder,
    data = cohort,
    samplestat = sstat,
    m = 3L,
    weight.method = "glm",
    variance = "robust"
  )
  oracle_log_hr <- unname(oracle$coefficients[1])

  # Both implement the same working-model approach (logistic model for the
  # selection probability); small formula differences cause sub-2e-2 deviation.
  expect_equal(our_log_hr, oracle_log_hr, tolerance = 2e-2)
})

# ---- Python oracle: GLM weight values match statsmodels Logit exactly ------

test_that("GLM inclusion weights match Python statsmodels oracle (1e-6)", {
  skip_if(
    !file.exists(
      testthat::test_path("fixtures", "python", "glm_weights_results.csv")
    )
  )
  cohort <- read.csv(
    testthat::test_path("fixtures", "python", "glm_weights_cohort.csv")
  )
  ncc <- data.table::as.data.table(
    read.csv(testthat::test_path("fixtures", "python", "glm_weights_ncc.csv"))
  )
  py <- read.csv(
    testthat::test_path("fixtures", "python", "glm_weights_results.csv")
  )

  ncc_glm <- compute_ncc_weights(
    ncc,
    cohort = cohort,
    method = "glm",
    time = "t"
  )

  # Compare per-cohort-subject ipw_weight for each unique subject in the NCC.
  r_unique <- unique(as.data.frame(ncc_glm)[, c(".cohort_row", "ipw_weight")])
  merged <- merge(
    r_unique,
    py[, c(".cohort_row", "ipw_weight")],
    by = ".cohort_row",
    suffixes = c("_r", "_py")
  )
  diff <- abs(merged$ipw_weight_r - merged$ipw_weight_py)
  # Both sides apply the identical product formula on the same logistic fit;
  # numerical agreement should be well within double-precision rounding.
  expect_lt(max(diff), 1e-6)
})

# ---- rejection: Phase-1-missing (cohort = NULL or missing time column) -----

test_that("compute_ncc_weights aborts with matchatr_missing_phase1 when cohort = NULL", {
  cohort <- make_ipw_cohort(n = 400L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE)
  )
  expect_error(
    compute_ncc_weights(ncc, cohort = NULL, method = "glm", time = "t"),
    class = "matchatr_missing_phase1"
  )
})

test_that("compute_ncc_weights aborts with matchatr_missing_phase1 when time col missing", {
  cohort <- make_ipw_cohort(n = 400L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE)
  )
  bad_cohort <- cohort
  names(bad_cohort)[names(bad_cohort) == "t"] <- "time_renamed"
  expect_error(
    compute_ncc_weights(ncc, cohort = bad_cohort, method = "glm", time = "t"),
    class = "matchatr_missing_phase1"
  )
})

test_that("compute_ncc_weights aborts with matchatr_bad_input when .cohort_row is absent", {
  cohort <- make_ipw_cohort(n = 400L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE)
  )
  ncc$.cohort_row <- NULL
  expect_error(
    compute_ncc_weights(ncc, cohort = cohort, method = "glm", time = "t"),
    class = "matchatr_bad_input"
  )
})

test_that("compute_ncc_weights aborts with matchatr_bad_input when ncc lacks required cols", {
  cohort <- make_ipw_cohort(n = 400L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE)
  )
  ncc$set <- NULL # remove required column
  expect_error(
    compute_ncc_weights(ncc, cohort = cohort, method = "glm", time = "t"),
    class = "matchatr_bad_input"
  )
})

# ---- ci_method = "sandwich" is accepted (re-computes via sandwich::sandwich) --

test_that("ci_method = 'sandwich' is accepted and gives finite CIs", {
  cohort <- make_ipw_cohort(n = 600L, seed = 5L)
  ncc <- withr::with_seed(
    2L,
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    "d",
    "exposure",
    nested_cc(strata = "set", time = "t"),
    estimator = "ipw_cox"
  )
  res <- contrast(fit, type = "hr", ci_method = "sandwich")
  expect_true(is.finite(res$contrasts$ci_lower))
  expect_true(is.finite(res$contrasts$ci_upper))
  expect_true(
    res$contrasts$ci_lower > 0 &&
      res$contrasts$ci_upper > res$contrasts$ci_lower
  )
})
