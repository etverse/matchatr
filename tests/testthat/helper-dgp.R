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

# Truth-based DGP for the polytomous (multinomial) subtype odds ratios. A cohort
# is drawn from a baseline-category multinomial model with KNOWN exposure log-ORs
# for each non-reference outcome group: caseA and caseB each contrast against the
# "control" reference. The exposure `x` has conditional log-OR `beta_a` in the
# caseA equation and `beta_b` in the caseB equation; a continuous confounder
# `age` carries its own per-equation slope. Separate group sampling would offset
# only the intercepts (the multinomial analogue of Prentice & Pyke 1979), so the
# slopes are recovered from the cohort itself, which is what the multinom oracle
# checks. Columns: g (factor: control/caseA/caseB), x (0/1), age. The true
# exposure log-ORs are returned in the "truth" attribute.
make_polytomous_cc <- function(
  n = 4000L,
  beta_a = log(2.5),
  beta_b = log(0.5),
  seed = 23L
) {
  sample <- withr::with_seed(seed, {
    x <- rbinom(n, 1L, 0.4)
    age <- rnorm(n, 50, 10)
    # Linear predictors for the two non-reference equations (control = baseline,
    # linear predictor fixed at 0). Distinct age slopes keep the equations from
    # collapsing onto each other.
    lp_a <- -0.4 + beta_a * x + 0.010 * (age - 50)
    lp_b <- -0.7 + beta_b * x + 0.025 * (age - 50)
    denom <- 1 + exp(lp_a) + exp(lp_b)
    p_ctrl <- 1 / denom
    p_a <- exp(lp_a) / denom
    u <- runif(n)
    g <- ifelse(
      u < p_ctrl,
      "control",
      ifelse(u < p_ctrl + p_a, "caseA", "caseB")
    )
    data.frame(
      g = factor(g, levels = c("control", "caseA", "caseB")),
      x = x,
      age = age,
      stringsAsFactors = FALSE
    )
  })
  rownames(sample) <- NULL
  attr(sample, "truth") <- c(caseA.x = beta_a, caseB.x = beta_b)
  sample
}

# A deterministic 3-group case-control table with exact cell counts, so the
# saturated multinomial coefficients AND their variances have a closed form
# (used by the polytomous subtype-OR and the homogeneity-test oracles). Rows:
# control (reference), caseA, caseB; columns: x = 0 / 1. For each non-reference
# group k the saturated multinomial log OR is the 2x2 Woolf value and its
# variance is the Woolf sum 1/n_{k1} + 1/n_{k0} + 1/n_{ref1} + 1/n_{ref0}; the
# caseA and caseB equations share the reference cells, so their log-OR
# covariance is 1/n_{ref1} + 1/n_{ref0}.
make_3group_table <- function(
  ctrl1 = 80L,
  ctrl0 = 120L,
  a1 = 60L,
  a0 = 40L,
  b1 = 30L,
  b0 = 70L
) {
  g <- c(
    rep("control", ctrl1 + ctrl0),
    rep("caseA", a1 + a0),
    rep("caseB", b1 + b0)
  )
  x <- c(
    rep(1L, ctrl1),
    rep(0L, ctrl0),
    rep(1L, a1),
    rep(0L, a0),
    rep(1L, b1),
    rep(0L, b0)
  )
  data.frame(
    g = factor(g, levels = c("control", "caseA", "caseB")),
    x = x,
    stringsAsFactors = FALSE
  )
}

# A deterministic 4-group case-control table with exact cell counts, extending
# `make_3group_table()` to three non-reference subtypes (caseA, caseB, caseC) so
# the homogeneity test runs with df = M - 1 = 2 and a closed-form oracle. The
# saturated multinomial reproduces the closed-form Woolf log-ORs and a 3x3
# covariance whose diagonal is the per-subtype Woolf sum and whose every
# off-diagonal is the shared reference contribution 1/n_ref1 + 1/n_ref0. The
# default counts give caseA / caseB a common OR of 2.25 and caseC an OR of 0.5,
# so the homogeneity test rejects. Rows: control (reference), caseA, caseB,
# caseC; columns: x = 0 / 1.
make_4group_table <- function(
  ctrl1 = 80L,
  ctrl0 = 120L,
  a1 = 60L,
  a0 = 40L,
  b1 = 45L,
  b0 = 30L,
  c1 = 20L,
  c0 = 60L
) {
  g <- c(
    rep("control", ctrl1 + ctrl0),
    rep("caseA", a1 + a0),
    rep("caseB", b1 + b0),
    rep("caseC", c1 + c0)
  )
  x <- c(
    rep(1L, ctrl1),
    rep(0L, ctrl0),
    rep(1L, a1),
    rep(0L, a0),
    rep(1L, b1),
    rep(0L, b0),
    rep(1L, c1),
    rep(0L, c0)
  )
  data.frame(
    g = factor(g, levels = c("control", "caseA", "caseB", "caseC")),
    x = x,
    stringsAsFactors = FALSE
  )
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
# exp(beta_x) regardless of the set-level baseline. A non-matching covariate `z`
# carries its own known conditional log-OR `beta_z` (default 0 = pure noise) so
# adjustment tests can check recovery, not just forwarding. Sets whose exposure
# is constant are uninformative and dropped by `clogit` (harmless).
# Columns: case (0/1), x (0/1 exposure), z (0/1 covariate), set (matched-set id).
make_matched_cc <- function(
  n_sets = 250L,
  ratio = 3L,
  beta_x = log(2.5),
  beta_z = 0,
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
      # with weight exp(x * beta_x + z * beta_z). The matched-CC analogue of the
      # risk-set sampling that makes clogit design-faithful; gives x (and z, when
      # beta_z != 0) a known conditional log-OR.
      w <- exp(x * beta_x + z * beta_z)
      case_idx <- sample.int(m, 1L, prob = w)
      case <- integer(m)
      case[case_idx] <- 1L
      data.frame(case = case, x = x, z = z, set = i)
    })
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    attr(out, "truth") <- c(beta_x = beta_x, beta_z = beta_z)
    out
  })
}

# Truth-based DGP for effect modification in matched case-control data. Each
# matched set has `ratio` controls and one case, drawn from the CONDITIONAL
# likelihood with a per-member selection weight exp(x * beta_{m}), so the
# exposure's conditional log-OR is `betas[level]` within each modifier level
# `m`. The conditional logistic model `case ~ x * m + strata(set)` therefore
# recovers beta_x = betas[ref] and beta_x + beta_{x:level} = betas[level].
#
# `within_set = FALSE` makes the modifier CONSTANT within each set (the modifier
# is itself a matching variable): the sets split into disjoint groups by level
# and, for 1:1 matching, each level's conditional likelihood reduces to McNemar
# on that level's discordant pairs -- the independent point + variance oracle.
# `within_set = TRUE` lets the modifier vary within a set, so its main effect is
# estimable (not aliased) and the more general interaction path is exercised.
# Columns: case (0/1), x (0/1 exposure), m (factor modifier), set (matched-set
# id). The true per-level log-ORs are returned in the "truth" attribute.
make_matched_cc_em <- function(
  n_sets = 300L,
  ratio = 1L,
  betas = c(a = log(2), b = log(5)),
  within_set = FALSE,
  seed = 41L
) {
  levs <- names(betas)
  m_size <- ratio + 1L
  out <- withr::with_seed(seed, {
    parts <- lapply(seq_len(n_sets), function(i) {
      mod <- if (within_set) {
        sample(levs, m_size, replace = TRUE)
      } else {
        rep(sample(levs, 1L), m_size)
      }
      # Set-level exposure prevalence (the matched-away nuisance).
      p_i <- stats::plogis(stats::rnorm(1, 0, 1))
      x <- stats::rbinom(m_size, 1L, p_i)
      # Conditional-likelihood case selection with a level-specific exposure
      # weight, giving x a known conditional log-OR within each modifier level.
      w <- exp(x * betas[mod])
      case_idx <- sample.int(m_size, 1L, prob = w)
      case <- integer(m_size)
      case[case_idx] <- 1L
      data.frame(case = case, x = x, m = mod, set = i, stringsAsFactors = FALSE)
    })
    do.call(rbind, parts)
  })
  rownames(out) <- NULL
  out$m <- factor(out$m, levels = levs)
  attr(out, "truth") <- betas
  out
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

# Truth-based DGP for the nested case-control (NCC) hazard ratio. A cohort is
# generated from a proportional-hazards model with a constant baseline hazard
# (exponential survival): the rate for subject i is base_rate * exp(beta_x * x +
# beta_z * z), so coxph(Surv(t, d) ~ x + z) on the FULL cohort recovers beta_x /
# beta_z. Continuous exponential times make ties negligible, so the sampled risk
# sets are clean. Administrative censoring at `tau`. Columns: id, t (observed
# time), d (event indicator), x (binary exposure), z (continuous confounder).
# The true Cox log hazard ratios are returned in the "truth" attribute. An NCC
# sample is drawn from this cohort with sample_ncc_riskset(); the conditional
# partial likelihood on that sample recovers the SAME beta (OR = HR exactly under
# risk-set sampling), which is the design-faithful oracle.
make_ncc_cohort <- function(
  n = 2500L,
  beta_x = log(2.2),
  beta_z = log(1.5),
  base_rate = 0.08,
  tau = 4,
  seed = 51L
) {
  cohort <- withr::with_seed(seed, {
    x <- stats::rbinom(n, 1L, 0.4)
    z <- stats::rnorm(n)
    # PH hazard: a constant baseline times exp(linear predictor). Inverting the
    # exponential survival gives T ~ Exp(rate = base_rate * exp(lp)).
    rate <- base_rate * exp(beta_x * x + beta_z * z)
    tt <- stats::rexp(n, rate)
    d <- as.integer(tt <= tau)
    t_obs <- pmin(tt, tau)
    data.frame(id = seq_len(n), t = t_obs, d = d, x = x, z = z)
  })
  attr(cohort, "truth") <- c(beta_x = beta_x, beta_z = beta_z)
  cohort
}

# A richer cohort for IPW NCC oracle-agreement tests: a continuous exposure `xc`,
# a continuous confounder `z1`, and a three-level factor confounder `z2`.
# Exponential rates keep hazards positive; the censoring cap at the upper
# quartile puts ~25% of subjects at the same exit time (heavy ties). The additive
# / AFT model fit to it need not be correctly specified — these tests check that
# matchatr and the external oracle agree on the fit for a complex covariate set,
# not that they recover a known parameter.
make_complex_ncc_cohort <- function(n = 3000L, seed = 101L) {
  withr::with_seed(seed, {
    xc <- stats::rnorm(n)
    z1 <- stats::rnorm(n)
    z2 <- factor(
      sample(c("a", "b", "c"), n, replace = TRUE),
      levels = c("a", "b", "c")
    )
    rate <- 0.05 * exp(0.3 * xc + 0.2 * z1 + 0.4 * (as.integer(z2) - 1L))
    tt <- stats::rexp(n, rate)
    tau <- stats::quantile(tt, 0.75)
    data.frame(
      id = seq_len(n),
      t = pmin(tt, tau),
      d = as.integer(tt <= tau),
      xc = xc,
      z1 = z1,
      z2 = z2
    )
  })
}

# Risk-set (incidence-density) NCC sampler used as a deterministic test fixture.
# Delegates to the exported sample_ncc() (the production risk-set sampler) under a
# fixed seed, so the analysis-path tests exercise the real generator rather than a
# parallel copy. For each case (d == 1) at its failure time, sample_ncc() draws m
# controls without replacement from those at risk then (t >= t_case, no left
# truncation); a subject sampled as a control may itself fail later in the cohort
# (the classical NCC structure), so `case` is the per-set indicator, not the
# cohort's `d`. Returns a data.frame in the historical column order: the cohort
# columns, then case (per-set 0/1), set (stratum id), risk_time (the failure time).
sample_ncc_riskset <- function(cohort, m = 2L, seed = 71L) {
  out <- withr::with_seed(
    seed,
    sample_ncc(cohort, time = "t", event = "d", m = m)
  )
  out <- as.data.frame(out)
  added <- c("case", "set", "risk_time")
  out <- out[, c(setdiff(names(out), added), added), drop = FALSE]
  rownames(out) <- NULL
  out
}

# Truth-based DGP for the case-cohort hazard ratio. A cohort is generated from
# a proportional-hazards model with a constant baseline hazard (exponential
# survival), then a random subcohort of size `sub_frac * n` is sampled at
# baseline. Every case is included (as in the standard case-cohort design).
# Columns: id, t (observed time), d (event indicator), x (binary exposure),
# z (continuous confounder), subcohort (0/1). The true Cox log hazard ratios
# are returned in the "truth" attribute.
make_case_cohort_data <- function(
  n = 3000L,
  beta_x = log(2),
  beta_z = log(1.5),
  base_rate = 0.06,
  tau = 5,
  sub_frac = 0.25,
  seed = 77L
) {
  cohort <- withr::with_seed(seed, {
    x <- stats::rbinom(n, 1L, 0.4)
    z <- stats::rnorm(n)
    rate <- base_rate * exp(beta_x * x + beta_z * z)
    tt <- stats::rexp(n, rate)
    d <- as.integer(tt <= tau)
    t_obs <- pmin(tt, tau)
    # Subcohort sampled once at baseline, independently of x and z.
    n_sub <- round(n * sub_frac)
    sc_idx <- sample.int(n, n_sub)
    subcohort <- integer(n)
    subcohort[sc_idx] <- 1L
    data.frame(
      id = seq_len(n),
      t = t_obs,
      d = d,
      x = x,
      z = z,
      subcohort = subcohort,
      stringsAsFactors = FALSE
    )
  })
  attr(cohort, "truth") <- c(beta_x = beta_x, beta_z = beta_z)
  cohort
}

# Truth-based DGP for the stratified case-cohort hazard ratio. Like
# make_case_cohort_data() but the subcohort is sampled *within* two strata
# (e.g. region A / B) at different fractions, which is the design Borgan I/II
# are built for. Columns: id, t, d, x, z, region (factor 2-level), subcohort
# (0/1, sampled within region). True Cox log-HRs are in the "truth" attribute.
make_stratified_case_cohort_data <- function(
  n = 3000L,
  beta_x = log(2),
  beta_z = log(1.5),
  base_rate = 0.06,
  tau = 5,
  sub_frac_a = 0.20,
  sub_frac_b = 0.35,
  seed = 88L
) {
  cohort <- withr::with_seed(seed, {
    region <- factor(sample(c("A", "B"), n, replace = TRUE))
    x <- stats::rbinom(n, 1L, 0.4)
    z <- stats::rnorm(n)
    rate <- base_rate * exp(beta_x * x + beta_z * z)
    tt <- stats::rexp(n, rate)
    d <- as.integer(tt <= tau)
    t_obs <- pmin(tt, tau)
    # Stratified subcohort: different sampling fractions per region.
    n_a <- sum(region == "A")
    n_b <- sum(region == "B")
    sc_a <- sample.int(n_a, round(n_a * sub_frac_a))
    sc_b <- sample.int(n_b, round(n_b * sub_frac_b))
    subcohort <- integer(n)
    idx_a <- which(region == "A")
    idx_b <- which(region == "B")
    subcohort[idx_a[sc_a]] <- 1L
    subcohort[idx_b[sc_b]] <- 1L
    data.frame(
      id = seq_len(n),
      t = t_obs,
      d = d,
      x = x,
      z = z,
      region = region,
      subcohort = subcohort,
      stringsAsFactors = FALSE
    )
  })
  attr(cohort, "truth") <- c(beta_x = beta_x, beta_z = beta_z)
  cohort
}

# Truth-based DGP for case-control-weighted MARGINAL causal effects. A cohort is
# drawn from a logistic outcome model with a binary exposure `x` confounded by a
# continuous `w` (x depends on w), so the conditional and marginal effects differ
# (non-collapsibility) and confounding must be adjusted. The g-formula truth is
# computed analytically on the FULL cohort using the true coefficients:
# m1 = E_w[expit(alpha + beta_x + gamma_w w)], m0 = E_w[expit(alpha + gamma_w w)],
# giving the marginal risk difference m1 - m0, risk ratio m1/m0, and marginal odds
# ratio (m1/(1-m1)) / (m0/(1-m0)). The CONDITIONAL odds ratio is exp(beta_x), which
# differs from the marginal OR. An unmatched case-control sample (every case + a
# `ratio`:1 random control sample) is returned; case-control sampling shifts only
# the intercept (Prentice & Pyke 1979), so case-control weighting back to the
# source prevalence q0 = P(Y=1) recovers the marginal truth. The "truth" attribute
# carries c(rd, rr, mor, cond_or); the "q0" attribute the source prevalence.
make_cohort_ccw <- function(
  n = 2e5,
  alpha = -3,
  beta_x = log(2.5),
  gamma_w = 1.0,
  ratio = 5L,
  seed = 13L
) {
  withr::with_seed(seed, {
    w <- stats::rnorm(n)
    # Exposure confounded by w, so adjustment is required and the marginal effect
    # is not the conditional one.
    x <- stats::rbinom(n, 1L, stats::plogis(0.2 * w))
    # Counterfactual risks under treat-all / treat-none, for the g-formula truth.
    p1 <- stats::plogis(alpha + beta_x + gamma_w * w)
    p0 <- stats::plogis(alpha + gamma_w * w)
    y <- stats::rbinom(n, 1L, stats::plogis(alpha + beta_x * x + gamma_w * w))
    m1 <- mean(p1)
    m0 <- mean(p0)
    truth <- c(
      rd = m1 - m0,
      rr = m1 / m0,
      mor = (m1 / (1 - m1)) / (m0 / (1 - m0)),
      cond_or = exp(beta_x)
    )
    q0 <- mean(y)
    coh <- data.frame(case = y, x = x, w = w)
    cases <- coh[coh$case == 1L, , drop = FALSE]
    controls_all <- coh[coh$case == 0L, , drop = FALSE]
    n_ctrl <- min(nrow(cases) * ratio, nrow(controls_all))
    controls <- controls_all[
      sample.int(nrow(controls_all), n_ctrl),
      ,
      drop = FALSE
    ]
    samp <- rbind(cases, controls)
    rownames(samp) <- NULL
    attr(samp, "truth") <- truth
    attr(samp, "q0") <- q0
    samp
  })
}

# Truth-based DGP for the DOUBLE-ROBUSTNESS of CCW-AIPW. Two cohorts, each with a
# single confounder `w` and a constant conditional treatment effect (0.7 on the
# logit), drawn so that exactly ONE of the two `~ w`-linear working models is
# correctly specified — the misspecification is a missing quadratic term:
#
#   kind = "out_wrong":  propensity is linear in w (a ~ expit(0.9 w)), so a `~ w`
#     propensity model is CORRECT; the outcome is nonlinear in w (a w^2 term), so
#     a `~ w` outcome model is WRONG. CCW-IPW (correct propensity) and CCW-AIPW
#     (doubly robust) recover the truth; CCW-g-formula (wrong outcome) is biased.
#   kind = "prop_wrong": the outcome is linear in w, so a `~ w` outcome model is
#     CORRECT; the propensity is nonlinear in w (a w^2 term), so a `~ w`
#     propensity model is WRONG. CCW-g-formula (correct outcome) and CCW-AIPW
#     recover the truth; CCW-IPW (wrong propensity) is biased.
#
# Both pass `confounders = ~ w` to matcha(), so the misspecification is a
# functional-form one (the working models cannot see w^2). The marginal
# risk-difference truth is the g-formula on the full cohort under the TRUE
# outcome model, m1 - m0 with m_a = E_w[expit(lp_a)], which does not depend on
# the propensity. An unmatched case-control sample (every case + a `ratio`:1
# control sample) is returned with the "truth" (marginal RD) and "q0" attributes.
make_dr_cohort_ccw <- function(
  kind = c("out_wrong", "prop_wrong"),
  n = 1.5e5,
  ratio = 5L,
  seed = 20L
) {
  kind <- match.arg(kind)
  withr::with_seed(seed, {
    w <- stats::rnorm(n)
    b0 <- -2.3
    if (kind == "out_wrong") {
      # Strong selection on w (so the treated / control w-distributions differ
      # sharply) plus a strong w^2 outcome term: the `~ w` outcome model's failure
      # to capture the curvature then biases the standardized risk difference
      # clearly (a marginal RD is otherwise fairly robust to mild outcome
      # misspecification because the curvature averages out).
      a <- stats::rbinom(n, 1L, stats::plogis(2.0 * w)) # linear propensity: correct
      lp <- function(av) b0 + 0.7 * av + 2.0 * w - 2.5 * w^2 # nonlinear outcome: wrong
    } else {
      a <- stats::rbinom(n, 1L, stats::plogis(-0.6 + 0.7 * w + 1.4 * w^2)) # nonlinear propensity: wrong
      lp <- function(av) b0 + 0.7 * av + 1.0 * w # linear outcome: correct
    }
    y <- stats::rbinom(n, 1L, stats::plogis(lp(a)))
    truth_rd <- mean(stats::plogis(lp(1)) - stats::plogis(lp(0)))
    q0 <- mean(y)
    coh <- data.frame(case = y, x = a, w = w)
    cases <- coh[coh$case == 1L, , drop = FALSE]
    controls_all <- coh[coh$case == 0L, , drop = FALSE]
    n_ctrl <- min(nrow(cases) * ratio, nrow(controls_all))
    controls <- controls_all[
      sample.int(nrow(controls_all), n_ctrl),
      ,
      drop = FALSE
    ]
    samp <- rbind(cases, controls)
    rownames(samp) <- NULL
    attr(samp, "truth") <- truth_rd
    attr(samp, "q0") <- q0
    samp
  })
}

# Counter-matched NCC sampler used as a deterministic test fixture.
# Delegates to the exported sample_ncc_counter_matched() under a fixed seed.
# When cohort is NULL, builds from make_ncc_cohort() and attaches z_bin = x as
# the binary surrogate (maximally correlated with the true exposure, so counter-
# matching concentrates all study resources in the exposure-surrogate boundary).
# Returns a data.frame (like sample_ncc_riskset) for column-order stability.
sample_ncc_counter_matched_fixture <- function(
  cohort = NULL,
  m = 1L,
  seed = 83L
) {
  if (is.null(cohort)) {
    cohort <- make_ncc_cohort()
    cohort$z_bin <- cohort$x
  }
  out <- withr::with_seed(
    seed,
    sample_ncc_counter_matched(
      cohort,
      time = "t",
      event = "d",
      surrogate = "z_bin",
      m = m
    )
  )
  out <- as.data.frame(out)
  added <- c("case", "set", "risk_time", "log_w")
  out <- out[, c(setdiff(names(out), added), added), drop = FALSE]
  rownames(out) <- NULL
  out
}
