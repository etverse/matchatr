---
paths:
  - "R/case_cohort.R"
  - "R/ipw_cox.R"
  - "R/weighted_cox.R"
  - "R/aft_ncc.R"
  - "R/additive_ncc.R"
  - "R/lin_ying.R"
  - "R/aalen_cumulative.R"
  - "R/excess_risk.R"
  - "R/absolute_risk*.R"
  - "R/risk_set_sampling.R"
  - "R/multi_endpoint.R"
  - "R/causal_survival_sampled.R"
  - "R/contrast_surv_sampled.R"
  - "R/variance_surv_sampled.R"
  - "tests/testthat/test-case_cohort.R"
  - "tests/testthat/test-ipw_ncc.R"
  - "tests/testthat/test-absolute_risk*.R"
  - "tests/testthat/test-causal_survival_sampled.R"
  - "tests/testthat/test-weighted_cox.R"
  - "tests/testthat/helper-ncc-oracle.R"
---

# Nested case-control, case-cohort and design-weighted survival

These are deliberate. Do not report them as bugs without a numerical reproducer.

- The case-cohort pseudo-likelihood is not a likelihood: Prentice and Borgan estimators reuse controls across failure times, so the score terms are dependent, the information matrix gives no valid SE and likelihood-ratio statistics are not χ². The variance is the Self-Prentice or Borgan asymptotic one or a robust sandwich, from `survival::cch`.
- For a stratified case-cohort design (Borgan II) the robust sandwich overstates the variance substantially (Borgan et al. 2018, section 16.4.4); the plug-in asymptotic estimator is preferred, so a gap between them is expected.
- `fit_ipw_cox()` uses `ties = "breslow"`, not coxph's Efron default, so the coefficients agree with the Breslow baseline hazard used for absolute risk (`ipw_breslow_ncc()`), checked to machine precision against `survival::survfit` under heavy ties. Under incidence-density sampling, failure times are usually distinct and the two coincide.
- Absolute-risk linear predictors for `ipw_cox` are built from the fitted model's terms (`ar_lp_from_newdata()`: `model.matrix(delete.response(terms(model)), model.frame(..., xlev = model$xlevels))`), so a data-dependent basis such as `poly()`, `ns()` or `scale()` is reproduced from the fit's `predvars`; recomputing it from `newdata` gives a different basis under the same coefficient names. The `cch` path uses the original formula, because `survival::cch` has no `predvars` that match its coefficients; keep the two paths separate.
- Design-weighted marginal survival (`surv_gcomp`) g-computes on matchatr's own weighted Cox (`fit_cch()`, `ipw_breslow_cch()`, with a cohort-resampling bootstrap) and does not call survatr. A full-cohort truth showed that survatr's pooled-logistic route biases the contrast under sampling: constant per-subject weights cannot express the time-varying case-cohort risk-set weighting, which attenuates β by 20 to 30 Monte Carlo SEs. The reasoning is in the roxygen comment of `fit_surv_gcomp()` in `R/causal_survival_sampled.R`.
- Known bug: the model SE of absolute risk under nested case-control IPW is too large. In 150 simulated cohorts it was 0.065 against an empirical SD of 0.039. `assemble_absolute_risk()` computes `sqrt(x'Vx + var_log_Λ₀)` without the β̂–Λ̂₀ covariance or the second-phase sampling term, on top of a robust Cox variance that is conservative for Samuelsen weights. No test compares an SE with Monte Carlo spread yet.
- Before reporting a numerical problem, compare nested case-control and case-cohort results with `multipleNCC` and `survival::cch`. `survival`, `multipleNCC`, `survey` and `Epi` are oracles or delegated engines.
