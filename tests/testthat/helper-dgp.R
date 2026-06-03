# Deterministic sample builders for the design / dispatch / validation tests.
# Phase 1 ships no estimator, so these carry no estimation truth -- they exist
# to give the design layer realistic case-control / NCC / case-cohort data to
# validate and route. Later phases replace / extend this with a cohort DGP that
# has known causal truth and samples designs from it.

# A matched case-control / NCC sample: `n_sets` matched sets, each with exactly
# one case and `ratio` controls, so every set is informative for the
# conditional likelihood. Columns: case (0/1), x (binary exposure), age, smoke,
# set (matched-set id), t (a time scale for NCC / case-cohort).
make_cc_data <- function(n_sets = 30L, ratio = 2L, seed = 101L) {
  withr::with_seed(seed, {
    n_cases <- n_sets
    n_controls <- n_sets * ratio
    n <- n_cases + n_controls
    data.frame(
      case = c(rep(1L, n_cases), rep(0L, n_controls)),
      x = rbinom(n, 1L, 0.4),
      age = round(rnorm(n, 55, 9), 1),
      smoke = rbinom(n, 1L, 0.3),
      set = c(seq_len(n_sets), rep(seq_len(n_sets), each = ratio)),
      t = round(rexp(n, 0.1), 2),
      stringsAsFactors = FALSE
    )
  })
}

# A tiny matched sample whose set 2 holds only cases (no control), so the
# conditional likelihood drops it -- used to exercise the uninformative-stratum
# warning. Sets 1 and 3 are informative discordant pairs pointing in OPPOSITE
# directions (case exposed in set 1, unexposed in set 3), so the conditional
# MLE on the remaining sets is finite (beta -> 0) and clogit converges cleanly:
# the warning under test is the uninformative-stratum one, not an incidental
# non-convergence warning from a single separating pair.
make_uninformative_cc <- function() {
  data.frame(
    case = c(1L, 0L, 1L, 1L, 1L, 0L),
    x = c(1L, 0L, 1L, 0L, 0L, 1L),
    set = c(1L, 1L, 2L, 2L, 3L, 3L)
  )
}

# Truth-based DGP for the unmatched case-control conditional OR. A cohort is
# generated from a logistic risk model with KNOWN slopes, then an unmatched
# case-control sample is drawn by taking every case and a random `ratio`:1
# sample of controls. Case-control sampling shifts only the intercept (Prentice
# & Pyke, 1979), so the conditional OR for `x` recovers exp(beta_x) and the
# `age` slope recovers beta_age. Returned attributes carry the true parameters.
make_cohort_cc <- function(
  n = 1e5,
  alpha = -4,
  beta_x = log(2),
  beta_age = 0.03,
  ratio = 1L,
  seed = 7L
) {
  sample <- withr::with_seed(seed, {
    x <- rbinom(n, 1L, 0.3)
    age <- rnorm(n, 50, 10)
    # Center age so alpha sets the baseline log-odds at the mean age.
    lp <- alpha + beta_x * x + beta_age * (age - 50)
    case <- rbinom(n, 1L, plogis(lp))
    cohort <- data.frame(case = case, x = x, age = age)
    cases <- cohort[cohort$case == 1L, , drop = FALSE]
    controls_all <- cohort[cohort$case == 0L, , drop = FALSE]
    n_ctrl <- min(nrow(cases) * ratio, nrow(controls_all))
    controls <- controls_all[sample.int(nrow(controls_all), n_ctrl), ]
    rbind(cases, controls)
  })
  rownames(sample) <- NULL
  attr(sample, "truth") <- c(beta_x = beta_x, beta_age = beta_age)
  sample
}

# A case-control sample with a THREE-level categorical exposure (low/med/high)
# and a continuous confounder. The log-ORs versus the "low" reference are
# beta_med and beta_high, recovered (vs the glm oracle) per non-reference level.
make_categorical_cc <- function(
  n = 1500L,
  beta_med = 0.4,
  beta_high = 0.9,
  seed = 7L
) {
  withr::with_seed(seed, {
    lvl <- sample(c("low", "med", "high"), n, replace = TRUE)
    age <- rnorm(n, 55, 10)
    b <- c(low = 0, med = beta_med, high = beta_high)
    lp <- -1 + b[lvl] + 0.02 * (age - 55)
    case <- rbinom(n, 1L, plogis(lp))
    data.frame(
      case = case,
      x = factor(lvl, levels = c("low", "med", "high")),
      age = age,
      stringsAsFactors = FALSE
    )
  })
}

# Expand the grouped Ille-et-Vilaine esophageal-cancer case-control data
# (datasets::esoph: one row per agegp x alcgp x tobgp cell with case/control
# counts) into individual-level rows. The alcohol exposure is returned as an
# UNORDERED factor `alc` (esoph's alcgp is ordered, which would trigger
# polynomial contrasts); agegp / tobgp are kept as confounders.
expand_esoph <- function() {
  es <- datasets::esoph
  rows <- do.call(
    rbind,
    lapply(seq_len(nrow(es)), function(i) {
      data.frame(
        case = c(rep(1L, es$ncases[i]), rep(0L, es$ncontrols[i])),
        agegp = es$agegp[i],
        tobgp = es$tobgp[i],
        alcgp = es$alcgp[i]
      )
    })
  )
  rows$alc <- factor(as.character(rows$alcgp), levels = levels(es$alcgp))
  rows
}

# A stratified case-control sample with a binary exposure, for the
# Mantel-Haenszel oracle (cross-checked against stats::mantelhaen.test). Each
# stratum k has its own baseline log-odds; `x` carries a common log-OR of 0.7.
# Columns: case (0/1), x (0/1), agegrp (K-level factor), sex (2-level factor).
make_stratified_cc <- function(n_strata = 5L, beta_x = 0.7, seed = 10L) {
  withr::with_seed(seed, {
    do.call(
      rbind,
      lapply(seq_len(n_strata), function(k) {
        n <- sample(80:160, 1)
        x <- rbinom(n, 1L, 0.4)
        case <- rbinom(n, 1L, plogis(-1 + beta_x * x + 0.3 * (k - 3)))
        data.frame(
          case = case,
          x = x,
          agegrp = factor(k, levels = seq_len(n_strata)),
          sex = factor(sample(c("M", "F"), n, replace = TRUE)),
          stringsAsFactors = FALSE
        )
      })
    )
  })
}

# Truth-based DGP for the matched case-control conditional OR. Each of `n_sets`
# matched sets has `ratio` controls and exactly one case, generated from the
# CONDITIONAL likelihood itself: `ratio + 1` exposures are drawn at a set-level
# exposure prevalence (the matched-away nuisance), then the case is selected
# with probability proportional to exp(x * beta_x) (Breslow & Day 1980, the
# 1:M conditional-likelihood construction). The CMLE of `x` therefore recovers
# exp(beta_x) regardless of the set-level baseline. A binary nuisance covariate
# `z` (non-matching confounder) is carried for adjustment tests. Sets whose
# exposure is constant are uninformative and dropped by `clogit` (harmless).
# Columns: case (0/1), x (0/1 exposure), z (0/1 covariate), set (matched-set id).
make_matched_cc <- function(
  n_sets = 250L,
  ratio = 3L,
  beta_x = log(2.5),
  seed = 11L
) {
  withr::with_seed(seed, {
    m <- ratio + 1L
    parts <- lapply(seq_len(n_sets), function(i) {
      # Set-level exposure prevalence: the matching variable's effect, removed
      # by conditioning on the set, so it never biases the CMLE.
      p_i <- stats::plogis(stats::rnorm(1, 0, 1))
      x <- stats::rbinom(m, 1L, p_i)
      z <- stats::rbinom(m, 1L, 0.5)
      # Conditional-likelihood case selection: exactly one case per set, drawn
      # with weight exp(x * beta_x). This is the matched-CC analogue of the
      # risk-set sampling that makes clogit the design-faithful estimator.
      w <- exp(x * beta_x)
      case_idx <- sample.int(m, 1L, prob = w)
      case <- integer(m)
      case[case_idx] <- 1L
      data.frame(case = case, x = x, z = z, set = i)
    })
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    attr(out, "truth") <- c(beta_x = beta_x)
    out
  })
}

# A deterministic 2x2 case-control table as a data frame, for the closed-form
# odds-ratio / Woolf-variance oracle. With these cell counts the OR is exactly
# (n11 * n00) / (n10 * n01) = (60 * 70) / (40 * 30) = 3.5, and a saturated
# logistic model reproduces both the OR and the Woolf log-OR variance
# (1/n11 + 1/n10 + 1/n01 + 1/n00) exactly.
make_2x2_cc <- function(n11 = 60L, n10 = 40L, n01 = 30L, n00 = 70L) {
  data.frame(
    case = c(rep(1L, n11 + n10), rep(0L, n01 + n00)),
    x = c(rep(1L, n11), rep(0L, n10), rep(1L, n01), rep(0L, n00))
  )
}
