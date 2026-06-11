# The unmatched case-control logistic engine: matcha() fits stats::glm, and
# contrast(type = "or") / tidy() / summary() report the conditional odds ratio.
# Oracles: stats::glm itself (pass-through), the closed-form 2x2 OR with Woolf
# variance, the sandwich estimator, and a cohort DGP with a known log-OR.

# --- matcha now runs the engine -----------------------------------------

test_that("matcha fits a binomial glm for the logistic engine", {
  df <- make_cohort_cc()
  fit <- matcha(df, "case", "x", unmatched_cc(), confounders = ~age)
  expect_s3_class(fit$model, "glm")
  expect_identical(stats::family(fit$model)$family, "binomial")

  # Pass-through: matchatr's wrapper must reproduce a hand-fit glm exactly.
  oracle <- stats::glm(case ~ x + age, stats::binomial(), df)
  expect_equal(stats::coef(fit$model), stats::coef(oracle))
})

test_that("an engine with no wired estimator leaves model NULL", {
  df <- make_cc_data()
  # The CCW IPW / AIPW / TMLE engines are not yet wired; their fit carries
  # model = NULL (ccw_gformula is wired, exercised in test-ccw.R).
  fit <- matcha(
    df,
    "case",
    "x",
    unmatched_cc(prevalence = 0.05),
    estimator = "ccw_ipw"
  )
  expect_null(fit$model)
})

# --- truth-based: conditional OR recovers the cohort slope --------------

test_that("the conditional OR recovers the cohort log-OR (truth-based)", {
  df <- make_cohort_cc(beta_x = log(2), beta_age = 0.03)
  truth <- attr(df, "truth")
  fit <- matcha(df, "case", "x", unmatched_cc(), confounders = ~age)
  res <- contrast(fit, type = "or")

  # Case-control sampling preserves slopes: exp(beta_x) ~ 2 and the age slope
  # ~ 0.03 are recovered up to sampling error (large cohort, fixed seed).
  expect_equal(
    log(res$contrasts$estimate),
    unname(truth["beta_x"]),
    tolerance = 0.1
  )
  td <- tidy(fit)
  age_slope <- td$estimate[td$term == "age"]
  # Absolute check: the age slope (~0.003 SE here) recovers 0.03 to within 0.01.
  expect_lt(abs(age_slope - unname(truth["beta_age"])), 0.01)
})

# --- closed-form 2x2 oracle (OR + Woolf variance) -----------------------

test_that("the unadjusted OR matches the closed-form 2x2 with Woolf CI", {
  df <- make_2x2_cc(n11 = 60L, n10 = 40L, n01 = 30L, n00 = 70L)
  fit <- matcha(df, "case", "x", unmatched_cc())
  res <- contrast(fit, type = "or")

  or_closed <- (60 * 70) / (40 * 30) # 3.5
  expect_equal(res$contrasts$estimate, or_closed, tolerance = 1e-6)

  # A saturated logistic model reproduces the Woolf log-OR variance, so the
  # model-based Wald interval matches the closed-form Woolf interval (to glm's
  # IRLS convergence precision).
  woolf_se <- sqrt(1 / 60 + 1 / 40 + 1 / 30 + 1 / 70)
  z <- stats::qnorm(0.975)
  expect_equal(
    res$contrasts$ci_lower,
    exp(log(or_closed) - z * woolf_se),
    tolerance = 1e-5
  )
  expect_equal(
    res$contrasts$ci_upper,
    exp(log(or_closed) + z * woolf_se),
    tolerance = 1e-5
  )
})

# --- point/CI agreement with the glm oracle -----------------------------

test_that("adjusted OR and Wald CI match stats::glm exactly", {
  df <- make_cohort_cc()
  fit <- matcha(df, "case", "x", unmatched_cc(), confounders = ~age)
  res <- contrast(fit, type = "or")

  oracle <- stats::glm(case ~ x + age, stats::binomial(), df)
  b <- stats::coef(oracle)["x"]
  se <- sqrt(diag(stats::vcov(oracle)))["x"]
  z <- stats::qnorm(0.975)
  expect_equal(res$contrasts$estimate, unname(exp(b)), tolerance = 1e-10)
  expect_equal(res$contrasts$se, unname(exp(b) * se), tolerance = 1e-10)
  expect_equal(
    res$contrasts$ci_lower,
    unname(exp(b - z * se)),
    tolerance = 1e-10
  )
  expect_equal(
    res$contrasts$ci_upper,
    unname(exp(b + z * se)),
    tolerance = 1e-10
  )
})

test_that("the sandwich CI matches the robust glm variance", {
  df <- make_cohort_cc()
  fit <- matcha(df, "case", "x", unmatched_cc(), confounders = ~age)
  res <- contrast(fit, type = "or", ci_method = "sandwich")

  oracle <- stats::glm(case ~ x + age, stats::binomial(), df)
  se_robust <- sqrt(diag(sandwich::sandwich(oracle)))["x"]
  expect_equal(
    res$contrasts$se,
    unname(exp(stats::coef(oracle)["x"]) * se_robust),
    tolerance = 1e-10
  )
  # Sandwich and model-based SEs differ here (a real robustness check, not a
  # tautology).
  res_model <- contrast(fit, type = "or", ci_method = "model")
  expect_false(isTRUE(all.equal(res$contrasts$se, res_model$contrasts$se)))
})

test_that("conf_level widens / narrows the interval as expected", {
  df <- make_cohort_cc()
  fit <- matcha(df, "case", "x", unmatched_cc(), confounders = ~age)
  ci90 <- contrast(fit, type = "or", conf_level = 0.90)$contrasts
  ci99 <- contrast(fit, type = "or", conf_level = 0.99)$contrasts
  expect_gt(ci90$ci_lower, ci99$ci_lower)
  expect_lt(ci90$ci_upper, ci99$ci_upper)

  oracle <- stats::glm(case ~ x + age, stats::binomial(), df)
  b <- stats::coef(oracle)["x"]
  se <- sqrt(diag(stats::vcov(oracle)))["x"]
  expect_equal(
    ci90$ci_lower,
    unname(exp(b - stats::qnorm(0.95) * se)),
    tolerance = 1e-10
  )
})

# --- exposure encodings -------------------------------------------------

test_that("a two-level factor exposure gives the same OR as 0/1", {
  df <- make_cohort_cc()
  or_num <- contrast(
    matcha(df, "case", "x", unmatched_cc()),
    type = "or"
  )$contrasts$estimate

  dff <- df
  dff$x <- factor(
    ifelse(df$x == 1L, "exposed", "unexposed"),
    levels = c("unexposed", "exposed")
  )
  res_fac <- contrast(matcha(dff, "case", "x", unmatched_cc()), type = "or")
  expect_equal(nrow(res_fac$contrasts), 1L)
  expect_equal(res_fac$contrasts$estimate, or_num, tolerance = 1e-10)
})

test_that("a continuous exposure yields the per-unit OR", {
  df <- make_cohort_cc()
  # Use age as a continuous exposure (no confounders).
  fit <- matcha(df, "case", "age", unmatched_cc())
  res <- contrast(fit, type = "or")
  oracle <- stats::glm(case ~ age, stats::binomial(), df)
  expect_equal(
    res$contrasts$estimate,
    unname(exp(stats::coef(oracle)["age"])),
    tolerance = 1e-10
  )
})

test_that("a categorical (k>2) exposure yields one OR per non-reference level", {
  df <- make_categorical_cc()
  fit <- matcha(df, "case", "x", unmatched_cc(), confounders = ~age)
  res <- contrast(fit, type = "or")

  expect_equal(nrow(res$contrasts), 2L) # k = 3 levels -> 2 contrasts
  expect_identical(res$contrasts$comparison, c("xmed", "xhigh"))
  # Each OR is the level vs the recorded reference, matching glm exactly.
  expect_identical(res$reference, "low")
  oracle <- stats::glm(case ~ x + age, stats::binomial(), df)
  expect_equal(
    res$contrasts$estimate,
    unname(exp(stats::coef(oracle)[c("xmed", "xhigh")])),
    tolerance = 1e-10
  )
})

test_that("colliding factor-level coefficient names do not corrupt the OR", {
  # Exposure `ses` level `low` and confounder `se` level `slow` both produce the
  # glm coefficient name "seslow". The exposure OR / SE must be the `ses` term's,
  # selected by term position, not by the colliding coefficient name.
  withr::with_seed(2, {
    n <- 4000
    ses <- factor(sample(c("high", "low"), n, TRUE), levels = c("high", "low"))
    se <- factor(sample(c("fast", "slow"), n, TRUE), levels = c("fast", "slow"))
    case <- rbinom(
      n,
      1,
      plogis(-1 + 0.8 * (ses == "low") + 1.5 * (se == "slow"))
    )
  })
  df <- data.frame(case = case, ses = ses, se = se)
  fit <- matcha(df, "case", "ses", unmatched_cc(), confounders = ~se)
  res <- contrast(fit, type = "or")

  expect_equal(nrow(res$contrasts), 1L) # only the ses effect, not se
  oracle <- stats::glm(case ~ ses + se, stats::binomial(), df)
  ses_idx <- which(attr(stats::model.matrix(oracle), "assign") == 1L)
  expect_equal(
    res$contrasts$estimate,
    unname(exp(stats::coef(oracle)[ses_idx])),
    tolerance = 1e-10
  )
  # tidy() gives each colliding "seslow" its own (position-based) SE, not an NA
  # or a value borrowed from the first match.
  td <- tidy(fit)
  seslow_se <- td$std.error[td$term == "seslow"]
  expect_equal(length(seslow_se), 2L)
  expect_false(anyNA(seslow_se))
  expect_false(isTRUE(all.equal(seslow_se[1], seslow_se[2])))
})

test_that("the reference is the baseline used, not an unused declared level", {
  withr::with_seed(6, {
    n <- 1200
    lvl <- sample(c("med", "high"), n, TRUE)
    case <- rbinom(n, 1, plogis(-1 + 0.5 * (lvl == "high")))
  })
  # 'absent' is declared first but never occurs; glm's baseline is 'med'.
  df <- data.frame(
    case = case,
    x = factor(lvl, levels = c("absent", "med", "high"))
  )
  res <- contrast(matcha(df, "case", "x", unmatched_cc()), type = "or")
  expect_identical(res$reference, "med")
  expect_identical(res$contrasts$comparison, "xhigh")
})

test_that("an ordinal numeric exposure yields a single per-step trend OR", {
  df <- make_categorical_cc()
  # Integer scores 0/1/2 -> a single trend OR per one-level step.
  df$score <- as.integer(df$x) - 1L
  fit <- matcha(df, "case", "score", unmatched_cc(), confounders = ~age)
  res <- contrast(fit, type = "or")
  expect_equal(nrow(res$contrasts), 1L)
  oracle <- stats::glm(case ~ score + age, stats::binomial(), df)
  expect_equal(
    res$contrasts$estimate,
    unname(exp(stats::coef(oracle)["score"])),
    tolerance = 1e-10
  )
})

test_that("an ordered-factor exposure is rejected with guidance", {
  df <- make_categorical_cc()
  df$x <- factor(df$x, ordered = TRUE)
  # Polynomial contrasts are not per-level ORs; refuse before fitting (so no
  # wasted fit or misleading missing-data warning), pointing to the fix.
  expect_error(
    matcha(df, "case", "x", unmatched_cc()),
    class = "matchatr_bad_input"
  )
})

test_that("model_fn must be a function", {
  df <- make_categorical_cc()
  expect_error(
    matcha(df, "case", "x", unmatched_cc(), model_fn = "stats::glm"),
    class = "matchatr_bad_input"
  )
})

test_that("model_fn must accept a `family` argument", {
  df <- make_categorical_cc()
  no_family <- function(formula, data) {
    stats::glm(formula, family = stats::binomial(), data = data)
  }
  expect_error(
    matcha(df, "case", "x", unmatched_cc(), model_fn = no_family),
    class = "matchatr_bad_input"
  )
})

test_that("a non-binomial model_fn is rejected, not exponentiated", {
  df <- make_categorical_cc()
  # A fitter that ignores `family` and returns an OLS lm: its slope is not a
  # log odds ratio, so exp() of it must not be passed off as an OR.
  ols <- function(formula, family, data) stats::lm(formula, data = data)
  expect_error(
    matcha(df, "case", "x", unmatched_cc(), model_fn = ols),
    class = "matchatr_bad_model_fit"
  )
})

test_that("model_fn = gam reproduces glm when the confounder is linear", {
  skip_if_not_installed("mgcv")
  df <- make_categorical_cc()
  # gam with no smooth term is the same fit as glm: the exposure OR must match.
  res_glm <- contrast(
    matcha(df, "case", "x", unmatched_cc(), confounders = ~age),
    type = "or"
  )
  res_gam <- contrast(
    matcha(
      df,
      "case",
      "x",
      unmatched_cc(),
      confounders = ~age,
      model_fn = mgcv::gam
    ),
    type = "or"
  )
  expect_s3_class(res_gam, "matchatr_result")
  expect_equal(
    res_gam$contrasts$estimate,
    res_glm$contrasts$estimate,
    tolerance = 1e-6
  )
})

test_that("model_fn = gam with a smooth confounder runs under both CI methods", {
  skip_if_not_installed("mgcv")
  df <- make_categorical_cc()
  fit <- matcha(
    df,
    "case",
    "x",
    unmatched_cc(),
    confounders = ~ s(age),
    model_fn = mgcv::gam
  )
  expect_s3_class(fit$model, "gam")
  res_m <- contrast(fit, type = "or", ci_method = "model")
  res_s <- contrast(fit, type = "or", ci_method = "sandwich")
  expect_equal(nrow(res_m$contrasts), 2L)
  expect_true(all(is.finite(res_m$contrasts$estimate)))
  expect_true(all(res_m$contrasts$se > 0))
  expect_true(all(res_s$contrasts$se > 0))
})

test_that("tidy / summary on a gam fit exclude smooth-basis coefficients", {
  skip_if_not_installed("mgcv")
  df <- make_categorical_cc()
  fit <- matcha(
    df,
    "case",
    "x",
    unmatched_cc(),
    confounders = ~ s(age),
    model_fn = mgcv::gam
  )
  td <- tidy(fit)
  # Only parametric terms (intercept + the two exposure levels); no `s(age).k`.
  expect_setequal(td$term, c("(Intercept)", "xmed", "xhigh"))
  expect_false(any(grepl("^s\\(age\\)", td$term)))
  # Parametric SEs are finite (not a recycled / NA artefact).
  expect_true(all(is.finite(td$std.error)))
})

# Book value: the Ille-et-Vilaine esophageal-cancer case-control data (handbook
# Ch3). A categorical alcohol exposure adjusted for age and tobacco reproduces
# the canonical monotone dose-response, matching glm on the same expanded data.
test_that("the esoph alcohol odds ratios match glm (book-value oracle)", {
  rows <- expand_esoph()
  fit <- matcha(
    rows,
    "case",
    "alc",
    unmatched_cc(),
    confounders = ~ agegp + tobgp
  )
  res <- contrast(fit, type = "or")

  expect_identical(res$reference, "0-39g/day")
  expect_equal(nrow(res$contrasts), 3L)
  oracle <- stats::glm(case ~ alc + agegp + tobgp, stats::binomial(), rows)
  expect_equal(
    res$contrasts$estimate,
    unname(exp(stats::coef(oracle)[grep("^alc", names(stats::coef(oracle)))])),
    tolerance = 1e-8
  )
  # Monotone alcohol dose-response (the well-known esoph finding).
  expect_true(all(diff(res$contrasts$estimate) > 0))
  expect_gt(res$contrasts$estimate[1], 1) # even the lowest band is harmful
})

# --- result / table structure -------------------------------------------

test_that("the OR result carries log-scale estimates and OR-scale contrasts", {
  df <- make_cohort_cc()
  fit <- matcha(df, "case", "x", unmatched_cc(), confounders = ~age)
  res <- contrast(fit, type = "or")
  expect_s3_class(res, "matchatr_result")
  expect_identical(res$type, "or")
  expect_identical(res$estimand, "conditional OR")
  expect_identical(res$ci_method, "model")
  expect_identical(res$n, nrow(df))

  # estimates are on the log-odds scale, contrasts on the OR scale (exp).
  expect_s3_class(res$estimates, "data.table")
  expect_s3_class(res$contrasts, "data.table")
  expect_equal(res$contrasts$estimate, exp(res$estimates$estimate))
  # Only the exposure term is reported -- never the intercept.
  expect_false("(Intercept)" %in% res$contrasts$comparison)
  expect_identical(res$contrasts$comparison, "x")
})

test_that("tidy.matchatr_fit returns a broom-style coefficient table", {
  df <- make_cohort_cc()
  fit <- matcha(df, "case", "x", unmatched_cc(), confounders = ~age)

  td <- tidy(fit)
  expect_s3_class(td, "data.table")
  expect_named(
    td,
    c(
      "term",
      "estimate",
      "std.error",
      "statistic",
      "p.value",
      "conf.low",
      "conf.high"
    )
  )
  expect_true("(Intercept)" %in% td$term)

  # exponentiate flips the estimate and bounds to the OR scale; std.error stays
  # on the log-odds scale (broom convention).
  tde <- tidy(fit, exponentiate = TRUE)
  expect_equal(tde$estimate, exp(td$estimate))
  expect_equal(tde$std.error, td$std.error)
  expect_equal(tde$conf.low, exp(td$conf.low))

  # conf.int = FALSE drops the bounds.
  expect_false("conf.low" %in% names(tidy(fit, conf.int = FALSE)))

  # robust changes the SE.
  expect_false(isTRUE(all.equal(
    tidy(fit, robust = TRUE)$std.error,
    td$std.error
  )))
})

test_that("aliased terms get NA (not recycled) SEs, even under robust", {
  # An aliased term placed BEFORE an estimable one is the failing case: the
  # sandwich drops the aliased column, so a positional index shifted every SE
  # after it. SEs are now aligned by coefficient name.
  set.seed(11)
  n <- 300
  d <- data.frame(
    case = rep(c(1L, 0L), each = n / 2),
    x = stats::rbinom(n, 1, 0.4),
    z1 = stats::rnorm(n),
    z3 = stats::rnorm(n)
  )
  d$z2 <- d$z1 # exact collinearity -> z2 aliased and dropped by glm
  fit <- matcha(d, "case", "x", unmatched_cc(), confounders = ~ z1 + z2 + z3)
  oracle <- stats::glm(case ~ x + z1 + z2 + z3, stats::binomial(), d)
  get <- function(td, term, col) td[[col]][td$term == term]

  # No recycling warning under either variance source.
  expect_no_warning(tidy(fit))
  expect_no_warning(tidy(fit, robust = TRUE))
  td_m <- tidy(fit)
  td_r <- tidy(fit, robust = TRUE)

  # The aliased term carries an NA SE, not a borrowed one.
  expect_true(is.na(get(td_m, "z2", "std.error")))
  expect_true(is.na(get(td_r, "z2", "std.error")))

  # z3 comes after the aliased z2: its SE must be its own under both sources.
  expect_equal(
    get(td_m, "z3", "std.error"),
    unname(sqrt(diag(stats::vcov(oracle)))["z3"]),
    tolerance = 1e-10
  )
  expect_equal(
    get(td_r, "z3", "std.error"),
    unname(sqrt(diag(sandwich::sandwich(oracle)))["z3"]),
    tolerance = 1e-10
  )

  # contrast() indexes the exposure SE by name too: a 1x1 vcov named "x".
  res <- contrast(fit, type = "or", ci_method = "sandwich")
  se_x <- sqrt(diag(sandwich::sandwich(oracle)))["x"]
  expect_equal(
    res$contrasts$se,
    unname(exp(stats::coef(oracle)["x"]) * se_x),
    tolerance = 1e-10
  )
  expect_identical(rownames(res$vcov), "x")
})

test_that("tidy.matchatr_result tidies the contrasts", {
  df <- make_cohort_cc()
  res <- contrast(matcha(df, "case", "x", unmatched_cc()), type = "or")
  td <- tidy(res)
  expect_named(
    td,
    c("term", "estimate", "std.error", "type", "conf.low", "conf.high")
  )
  expect_identical(td$type, "or")
  expect_equal(td$estimate, res$contrasts$estimate)
})

test_that("summary prints the OR table and returns it invisibly", {
  df <- make_cohort_cc()
  fit <- matcha(df, "case", "x", unmatched_cc(), confounders = ~age)
  expect_output(summary(fit), "Conditional odds ratios")
  expect_output(summary(fit), "intercept is not an interpretable")
  out <- withVisible(summary(fit))
  expect_false(out$visible)
  expect_s3_class(out$value, "data.table")
})

# --- rejections ---------------------------------------------------------

test_that("RD / RR are rejected as unidentified from unmatched CC", {
  df <- make_cohort_cc()
  fit <- matcha(df, "case", "x", unmatched_cc(), confounders = ~age)
  expect_error(
    contrast(fit, type = "difference"),
    class = "matchatr_unidentified_estimand"
  )
  expect_error(
    contrast(fit, type = "ratio"),
    class = "matchatr_unidentified_estimand"
  )
})

test_that("rows with missing values are dropped and reported as the analysis n", {
  set.seed(42)
  n <- 300
  d <- data.frame(
    case = rep(c(1L, 0L), each = n / 2),
    x = stats::rbinom(n, 1, 0.4),
    age = stats::rnorm(n)
  )
  d$age[1:37] <- NA
  # The listwise deletion is surfaced at fit time.
  expect_warning(
    fit <- matcha(d, "case", "x", unmatched_cc(), confounders = ~age),
    class = "matchatr_dropped_rows"
  )
  res <- contrast(fit, type = "or")
  # The reported n is the complete-case count glm used, not the full sample.
  expect_equal(res$n, stats::nobs(fit$model))
  expect_equal(res$n, nrow(d) - 37L)
  expect_lt(res$n, nrow(d))
})

test_that("a constant (non-estimable) exposure is rejected", {
  # A constant exposure aliases to NA in glm; contrast() must abort, not return
  # an NA odds ratio.
  df_const <- data.frame(
    case = rep(c(1L, 0L), each = 100),
    x = rep(1L, 200),
    age = stats::rnorm(200)
  )
  fit_const <- matcha(df_const, "case", "x", unmatched_cc(), confounders = ~age)
  expect_error(
    contrast(fit_const, type = "or"),
    class = "matchatr_unestimable_exposure"
  )
})

test_that("bootstrap CI is rejected for the conditional OR", {
  df <- make_cohort_cc()
  fit <- matcha(df, "case", "x", unmatched_cc())
  expect_error(
    contrast(fit, type = "or", ci_method = "bootstrap"),
    class = "matchatr_unsupported_variance"
  )
})

test_that("logistic-contrast rejections read clearly", {
  df <- make_2x2_cc()
  fit <- matcha(df, "case", "x", unmatched_cc())
  expect_snapshot(contrast(fit, type = "difference"), error = TRUE)
  expect_snapshot(contrast(fit, type = "ratio"), error = TRUE)
  expect_snapshot(
    contrast(fit, type = "or", ci_method = "bootstrap"),
    error = TRUE
  )

  # Ordered-factor exposure: rejected at fit time; the guidance is user-facing.
  dord <- make_categorical_cc()
  dord$x <- factor(dord$x, ordered = TRUE)
  expect_snapshot(matcha(dord, "case", "x", unmatched_cc()), error = TRUE)
})
