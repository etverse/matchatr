# Case-cohort analysis: survival::cch pseudo-likelihood wrapper. fit_cch()
# builds the Surv(time, status) ~ x + covariates formula, subsets to
# cases + subcohort, and delegates to survival::cch. contrast_cch() reads the
# cch asymptotic variance (Self-Prentice / Lin-Ying / Borgan) and reports the
# exposure's hazard ratio. Oracles:
#   (1) survival::nwtco — the canonical survival::cch example; compare HR and
#       log-scale SE against a direct cch() call;
#   (2) make_case_cohort_data() DGP with known Cox log-HR for coverage;
#   (3) full-cohort coxph for the same DGP (agreement within sampling error).
# No Python (delicatessen) oracle: the Prentice / Borgan pseudo-likelihood is
# not an M-estimator and is not covered by delicatessen.

# --- oracle: survival::nwtco Prentice method ----------------------------------

test_that("fit_cch with Prentice matches survival::cch on nwtco", {
  skip_if_not_installed("survival")
  data(nwtco, package = "survival")

  # Direct survival::cch call on the pre-filtered data (cch requires only
  # cases + subcohort members; censored non-members must be excluded).
  nwtco2 <- nwtco[nwtco$in.subcohort | nwtco$rel == 1L, ]
  oracle <- survival::cch(
    survival::Surv(edrel, rel) ~ histol + stage + age,
    data = nwtco2,
    subcoh = ~in.subcohort,
    id = ~seqno,
    cohort.size = nrow(nwtco),
    method = "Prentice"
  )
  oracle_log_hr <- unname(coef(oracle)["histol"])
  oracle_se <- unname(sqrt(diag(vcov(oracle)))["histol"])

  fit <- matcha(
    nwtco,
    outcome = "rel",
    exposure = "histol",
    design = case_cohort(
      subcohort = "in.subcohort",
      time = "edrel",
      method = "Prentice",
      id = "seqno"
    ),
    confounders = ~ stage + age,
    estimator = "cch"
  )
  r <- contrast(fit, type = "hr")

  # contrast table reports on the HR scale; back out the log to compare.
  matchatr_log_hr <- log(r$contrasts$estimate)
  matchatr_se_hr <- r$contrasts$se # delta-method on HR scale
  matchatr_se_log <- matchatr_se_hr / r$contrasts$estimate # back to log scale

  expect_equal(matchatr_log_hr, oracle_log_hr, tolerance = 1e-6)
  expect_equal(matchatr_se_log, oracle_se, tolerance = 1e-6)
})

test_that("default contrast type for case_cohort is hr", {
  skip_if_not_installed("survival")
  data(nwtco, package = "survival")

  fit <- matcha(
    nwtco,
    outcome = "rel",
    exposure = "histol",
    design = case_cohort(
      subcohort = "in.subcohort",
      time = "edrel",
      method = "Prentice",
      id = "seqno"
    ),
    confounders = ~ stage + age,
    estimator = "cch"
  )
  r <- contrast(fit)
  expect_identical(r$type, "hr")
})

# --- oracle: nwtco SelfPrentice method ----------------------------------------

test_that("fit_cch with SelfPrentice matches survival::cch on nwtco", {
  skip_if_not_installed("survival")
  data(nwtco, package = "survival")

  nwtco2 <- nwtco[nwtco$in.subcohort | nwtco$rel == 1L, ]
  oracle <- survival::cch(
    survival::Surv(edrel, rel) ~ histol + stage + age,
    data = nwtco2,
    subcoh = ~in.subcohort,
    id = ~seqno,
    cohort.size = nrow(nwtco),
    method = "SelfPrentice"
  )
  oracle_log_hr <- unname(coef(oracle)["histol"])
  oracle_se <- unname(sqrt(diag(vcov(oracle)))["histol"])

  fit <- matcha(
    nwtco,
    outcome = "rel",
    exposure = "histol",
    design = case_cohort(
      subcohort = "in.subcohort",
      time = "edrel",
      method = "SelfPrentice",
      id = "seqno"
    ),
    confounders = ~ stage + age,
    estimator = "cch"
  )
  r <- contrast(fit, type = "hr")

  matchatr_log_hr <- log(r$contrasts$estimate)
  matchatr_se_log <- r$contrasts$se / r$contrasts$estimate

  expect_equal(matchatr_log_hr, oracle_log_hr, tolerance = 1e-6)
  expect_equal(matchatr_se_log, oracle_se, tolerance = 1e-6)
})

# --- oracle: nwtco LinYing method ---------------------------------------------

test_that("fit_cch with LinYing matches survival::cch on nwtco", {
  skip_if_not_installed("survival")
  data(nwtco, package = "survival")

  nwtco2 <- nwtco[nwtco$in.subcohort | nwtco$rel == 1L, ]
  oracle <- survival::cch(
    survival::Surv(edrel, rel) ~ histol + stage + age,
    data = nwtco2,
    subcoh = ~in.subcohort,
    id = ~seqno,
    cohort.size = nrow(nwtco),
    method = "LinYing"
  )
  oracle_log_hr <- unname(coef(oracle)["histol"])
  oracle_se <- unname(sqrt(diag(vcov(oracle)))["histol"])

  fit <- matcha(
    nwtco,
    outcome = "rel",
    exposure = "histol",
    design = case_cohort(
      subcohort = "in.subcohort",
      time = "edrel",
      method = "LinYing",
      id = "seqno"
    ),
    confounders = ~ stage + age,
    estimator = "cch"
  )
  r <- contrast(fit, type = "hr")

  matchatr_log_hr <- log(r$contrasts$estimate)
  matchatr_se_log <- r$contrasts$se / r$contrasts$estimate

  expect_equal(matchatr_log_hr, oracle_log_hr, tolerance = 1e-6)
  expect_equal(matchatr_se_log, oracle_se, tolerance = 1e-6)
})

# --- truth-based: DGP with known Cox log-HR -----------------------------------

test_that("fit_cch recovers known Cox log-HR within 3.5 SE (Prentice)", {
  cohort <- make_case_cohort_data()
  truth_beta_x <- unname(attr(cohort, "truth")["beta_x"])

  fit <- matcha(
    cohort,
    outcome = "d",
    exposure = "x",
    design = case_cohort(
      subcohort = "subcohort",
      time = "t",
      method = "Prentice",
      id = "id"
    ),
    confounders = ~z,
    estimator = "cch"
  )
  r <- contrast(fit, type = "hr")
  log_hr <- log(r$contrasts$estimate)
  se_log <- r$contrasts$se / r$contrasts$estimate

  expect_equal(log_hr, truth_beta_x, tolerance = 3.5 * se_log)
})

test_that("fit_cch recovers known Cox log-HR within 3.5 SE (LinYing)", {
  cohort <- make_case_cohort_data()
  truth_beta_x <- unname(attr(cohort, "truth")["beta_x"])

  fit <- matcha(
    cohort,
    outcome = "d",
    exposure = "x",
    design = case_cohort(
      subcohort = "subcohort",
      time = "t",
      method = "LinYing",
      id = "id"
    ),
    confounders = ~z,
    estimator = "cch"
  )
  r <- contrast(fit, type = "hr")
  log_hr <- log(r$contrasts$estimate)
  se_log <- r$contrasts$se / r$contrasts$estimate

  expect_equal(log_hr, truth_beta_x, tolerance = 3.5 * se_log)
})

# --- full-cohort coxph agreement ----------------------------------------------

test_that("fit_cch agrees with full-cohort coxph within combined SE", {
  skip_if_not_installed("survival")
  cohort <- make_case_cohort_data()

  fit_cc <- matcha(
    cohort,
    outcome = "d",
    exposure = "x",
    design = case_cohort(
      subcohort = "subcohort",
      time = "t",
      method = "Prentice",
      id = "id"
    ),
    confounders = ~z,
    estimator = "cch"
  )
  r_cc <- contrast(fit_cc, type = "hr")
  log_hr_cc <- log(r_cc$contrasts$estimate)
  se_cc <- r_cc$contrasts$se / r_cc$contrasts$estimate

  full <- survival::coxph(
    survival::Surv(t, d) ~ x + z,
    data = cohort
  )
  log_hr_full <- unname(coef(full)["x"])
  se_full <- unname(sqrt(diag(vcov(full)))["x"])

  # Both estimate the same Cox parameter; they should agree within the combined
  # SE (the case-cohort is noisier, so this is a loose oracle check).
  combined_se <- sqrt(se_cc^2 + se_full^2)
  expect_equal(log_hr_cc, log_hr_full, tolerance = 3.5 * combined_se)
})

# --- structural checks --------------------------------------------------------

test_that("fit_cch returns a cch object and matchatr_fit", {
  skip_if_not_installed("survival")
  data(nwtco, package = "survival")

  fit <- matcha(
    nwtco,
    outcome = "rel",
    exposure = "histol",
    design = case_cohort(
      subcohort = "in.subcohort",
      time = "edrel",
      method = "Prentice",
      id = "seqno"
    ),
    confounders = ~ stage + age,
    estimator = "cch"
  )
  expect_s3_class(fit, "matchatr_fit")
  expect_s3_class(fit$model, "cch")
})

test_that("contrast_cch returns a matchatr_result", {
  skip_if_not_installed("survival")
  data(nwtco, package = "survival")

  fit <- matcha(
    nwtco,
    outcome = "rel",
    exposure = "histol",
    design = case_cohort(
      subcohort = "in.subcohort",
      time = "edrel",
      method = "Prentice",
      id = "seqno"
    ),
    confounders = ~ stage + age,
    estimator = "cch"
  )
  r <- contrast(fit, type = "hr")
  expect_s3_class(r, "matchatr_result")
  expect_true(is.finite(r$contrasts$estimate))
  expect_true(r$contrasts$estimate > 0)
  expect_true(r$contrasts$ci_lower < r$contrasts$estimate)
  expect_true(r$contrasts$estimate < r$contrasts$ci_upper)
})

test_that("fit_cch does not mutate the input data", {
  cohort <- make_case_cohort_data()
  before <- names(cohort)
  matcha(
    cohort,
    outcome = "d",
    exposure = "x",
    design = case_cohort(
      subcohort = "subcohort",
      time = "t",
      method = "Prentice",
      id = "id"
    ),
    confounders = ~z,
    estimator = "cch"
  )
  expect_identical(names(cohort), before)
})

test_that("case_cohort design print includes method", {
  d <- case_cohort(subcohort = "in_sc", time = "t", method = "LinYing")
  out <- capture.output(print(d))
  expect_true(any(grepl("LinYing", out)))
})

# --- case_cohort() constructor rejections -------------------------------------

test_that("case_cohort rejects an invalid method", {
  expect_snapshot(
    error = TRUE,
    case_cohort(subcohort = "sc", time = "t", method = "BadMethod")
  )
})

test_that("case_cohort requires subcohort and time strings", {
  expect_error(
    case_cohort(subcohort = 1L, time = "t"),
    class = "matchatr_bad_input"
  )
  expect_error(
    case_cohort(subcohort = "sc", time = 1L),
    class = "matchatr_bad_input"
  )
})

test_that("case_cohort with id = non-string is rejected", {
  expect_error(
    case_cohort(subcohort = "sc", time = "t", id = 1L),
    class = "matchatr_bad_input"
  )
})

# --- contrast_cch() scale rejections ------------------------------------------

test_that("contrast_cch rejects type = or", {
  skip_if_not_installed("survival")
  data(nwtco, package = "survival")

  fit <- matcha(
    nwtco,
    outcome = "rel",
    exposure = "histol",
    design = case_cohort(
      subcohort = "in.subcohort",
      time = "edrel",
      method = "Prentice",
      id = "seqno"
    ),
    confounders = ~ stage + age,
    estimator = "cch"
  )
  expect_snapshot(
    error = TRUE,
    contrast(fit, type = "or")
  )
})

test_that("contrast_cch rejects type = difference", {
  skip_if_not_installed("survival")
  data(nwtco, package = "survival")

  fit <- matcha(
    nwtco,
    outcome = "rel",
    exposure = "histol",
    design = case_cohort(
      subcohort = "in.subcohort",
      time = "edrel",
      method = "Prentice",
      id = "seqno"
    ),
    confounders = ~ stage + age,
    estimator = "cch"
  )
  expect_error(
    contrast(fit, type = "difference"),
    class = "matchatr_unidentified_estimand"
  )
})

test_that("contrast_cch rejects ci_method = sandwich", {
  skip_if_not_installed("survival")
  data(nwtco, package = "survival")

  fit <- matcha(
    nwtco,
    outcome = "rel",
    exposure = "histol",
    design = case_cohort(
      subcohort = "in.subcohort",
      time = "edrel",
      method = "Prentice",
      id = "seqno"
    ),
    confounders = ~ stage + age,
    estimator = "cch"
  )
  expect_snapshot(
    error = TRUE,
    contrast(fit, ci_method = "sandwich")
  )
})

test_that("contrast_cch rejects ci_method = bootstrap", {
  skip_if_not_installed("survival")
  data(nwtco, package = "survival")

  fit <- matcha(
    nwtco,
    outcome = "rel",
    exposure = "histol",
    design = case_cohort(
      subcohort = "in.subcohort",
      time = "edrel",
      method = "Prentice",
      id = "seqno"
    ),
    confounders = ~ stage + age,
    estimator = "cch"
  )
  expect_error(
    contrast(fit, ci_method = "bootstrap"),
    class = "matchatr_unsupported_variance"
  )
})

# --- matcha() design-column validation ----------------------------------------

test_that("matcha rejects case_cohort with missing subcohort column", {
  cohort <- make_case_cohort_data()
  cohort$subcohort <- NULL

  expect_error(
    matcha(
      cohort,
      outcome = "d",
      exposure = "x",
      design = case_cohort(subcohort = "subcohort", time = "t"),
      estimator = "cch"
    ),
    class = "matchatr_bad_design"
  )
})

test_that("matcha rejects case_cohort with missing time column", {
  cohort <- make_case_cohort_data()
  cohort$t <- NULL

  expect_error(
    matcha(
      cohort,
      outcome = "d",
      exposure = "x",
      design = case_cohort(subcohort = "subcohort", time = "t"),
      estimator = "cch"
    ),
    class = "matchatr_bad_design"
  )
})

# --- oracle: nwtco I.Borgan and II.Borgan methods ------------------------------

test_that("fit_cch with I.Borgan matches survival::cch on nwtco", {
  skip_if_not_installed("survival")
  data(nwtco, package = "survival")

  # Borgan methods need per-stratum cohort sizes and a stratum formula.
  nwtco2 <- nwtco[nwtco$in.subcohort | nwtco$rel == 1L, ]
  strat_sizes <- table(nwtco$instit)
  oracle <- survival::cch(
    survival::Surv(edrel, rel) ~ histol + stage + age,
    data = nwtco2,
    subcoh = ~in.subcohort,
    id = ~seqno,
    cohort.size = strat_sizes,
    stratum = ~instit,
    method = "I.Borgan"
  )
  oracle_log_hr <- unname(coef(oracle)["histol"])
  oracle_se <- unname(sqrt(diag(vcov(oracle)))["histol"])

  fit <- matcha(
    nwtco,
    outcome = "rel",
    exposure = "histol",
    design = case_cohort(
      subcohort = "in.subcohort",
      time = "edrel",
      method = "I.Borgan",
      id = "seqno",
      stratum = "instit"
    ),
    confounders = ~ stage + age,
    estimator = "cch"
  )
  r <- contrast(fit, type = "hr")

  matchatr_log_hr <- log(r$contrasts$estimate)
  matchatr_se_log <- r$contrasts$se / r$contrasts$estimate

  expect_equal(matchatr_log_hr, oracle_log_hr, tolerance = 1e-6)
  expect_equal(matchatr_se_log, oracle_se, tolerance = 1e-6)
})

test_that("fit_cch with II.Borgan matches survival::cch on nwtco", {
  skip_if_not_installed("survival")
  data(nwtco, package = "survival")

  nwtco2 <- nwtco[nwtco$in.subcohort | nwtco$rel == 1L, ]
  strat_sizes <- table(nwtco$instit)
  oracle <- survival::cch(
    survival::Surv(edrel, rel) ~ histol + stage + age,
    data = nwtco2,
    subcoh = ~in.subcohort,
    id = ~seqno,
    cohort.size = strat_sizes,
    stratum = ~instit,
    method = "II.Borgan"
  )
  oracle_log_hr <- unname(coef(oracle)["histol"])
  oracle_se <- unname(sqrt(diag(vcov(oracle)))["histol"])

  fit <- matcha(
    nwtco,
    outcome = "rel",
    exposure = "histol",
    design = case_cohort(
      subcohort = "in.subcohort",
      time = "edrel",
      method = "II.Borgan",
      id = "seqno",
      stratum = "instit"
    ),
    confounders = ~ stage + age,
    estimator = "cch"
  )
  r <- contrast(fit, type = "hr")

  matchatr_log_hr <- log(r$contrasts$estimate)
  matchatr_se_log <- r$contrasts$se / r$contrasts$estimate

  expect_equal(matchatr_log_hr, oracle_log_hr, tolerance = 1e-6)
  expect_equal(matchatr_se_log, oracle_se, tolerance = 1e-6)
})

# --- truth-based: stratified case-cohort DGP with known Cox log-HR ------------

test_that("fit_cch I.Borgan recovers known log-HR within 3.5 SE", {
  cohort <- make_stratified_case_cohort_data()
  truth_beta_x <- unname(attr(cohort, "truth")["beta_x"])

  fit <- matcha(
    cohort,
    outcome = "d",
    exposure = "x",
    design = case_cohort(
      subcohort = "subcohort",
      time = "t",
      method = "I.Borgan",
      id = "id",
      stratum = "region"
    ),
    confounders = ~z,
    estimator = "cch"
  )
  r <- contrast(fit, type = "hr")
  log_hr <- log(r$contrasts$estimate)
  se_log <- r$contrasts$se / r$contrasts$estimate

  expect_equal(log_hr, truth_beta_x, tolerance = 3.5 * se_log)
})

test_that("fit_cch II.Borgan recovers known log-HR within 3.5 SE", {
  cohort <- make_stratified_case_cohort_data()
  truth_beta_x <- unname(attr(cohort, "truth")["beta_x"])

  fit <- matcha(
    cohort,
    outcome = "d",
    exposure = "x",
    design = case_cohort(
      subcohort = "subcohort",
      time = "t",
      method = "II.Borgan",
      id = "id",
      stratum = "region"
    ),
    confounders = ~z,
    estimator = "cch"
  )
  r <- contrast(fit, type = "hr")
  log_hr <- log(r$contrasts$estimate)
  se_log <- r$contrasts$se / r$contrasts$estimate

  expect_equal(log_hr, truth_beta_x, tolerance = 3.5 * se_log)
})

# --- Borgan rejections --------------------------------------------------------

test_that("Borgan I without stratum aborts with matchatr_bad_design", {
  skip_if_not_installed("survival")
  data(nwtco, package = "survival")

  expect_snapshot(
    error = TRUE,
    matcha(
      nwtco,
      outcome = "rel",
      exposure = "histol",
      design = case_cohort(
        subcohort = "in.subcohort",
        time = "edrel",
        method = "I.Borgan",
        id = "seqno"
      ),
      confounders = ~ stage + age,
      estimator = "cch"
    )
  )
})

test_that("Borgan II without stratum aborts with matchatr_bad_design", {
  cohort <- make_stratified_case_cohort_data()

  expect_error(
    matcha(
      cohort,
      outcome = "d",
      exposure = "x",
      design = case_cohort(
        subcohort = "subcohort",
        time = "t",
        method = "II.Borgan",
        id = "id"
      ),
      confounders = ~z,
      estimator = "cch"
    ),
    class = "matchatr_bad_design"
  )
})

test_that("case_cohort design print includes stratum when set", {
  d <- case_cohort(
    subcohort = "in_sc",
    time = "t",
    method = "I.Borgan",
    stratum = "region"
  )
  out <- capture.output(print(d))
  expect_true(any(grepl("region", out)))
  expect_true(any(grepl("I.Borgan", out)))
})

test_that("matcha rejects case_cohort with missing stratum column", {
  cohort <- make_stratified_case_cohort_data()
  cohort$region <- NULL

  expect_error(
    matcha(
      cohort,
      outcome = "d",
      exposure = "x",
      design = case_cohort(
        subcohort = "subcohort",
        time = "t",
        method = "I.Borgan",
        id = "id",
        stratum = "region"
      ),
      confounders = ~z,
      estimator = "cch"
    ),
    class = "matchatr_bad_design"
  )
})
