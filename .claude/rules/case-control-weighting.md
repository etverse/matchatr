---
paths:
  - "R/ccw*.R"
  - "R/tmle_ccw.R"
  - "R/variance_ccw.R"
  - "R/weights_cc.R"
  - "R/summary.R"
  - "R/tidy.R"
  - "tests/testthat/test-ccw.R"
  - "tests/testthat/test-tmle_ccw.R"
  - "tests/testthat/test-variance_ccw.R"
  - "tests/testthat/helper-tmle-oracle.R"
---

# Case-control weighting

These are deliberate. Do not report them as bugs without a numerical reproducer.

- The g-formula, IPW and AIPW case-control-weighting estimators call causatr with the q₀ weights as `weights` (`fit_ccw()` maps `ccw_gformula`, `ccw_ipw`, `ccw_aipw` to `causat(estimator = "gcomp", "ipw", "aipw")`), and all three need `confounders`. causatr has no targeted learning, so the targeting step of `ccw_tmle` is matchatr's own code.
- `fit_ccw()` fits the outcome model with `family = "quasibinomial"`: the mean model and sandwich are those of the binomial, without the warning about non-integer successes that fractional weights cause. The `family` argument governs the outcome or marginal-mean model of all three estimators; causatr detects the propensity family itself.
- A case-control-weighting fit's `model` is a `causatr_fit`, which has no `coef()` or `vcov()`. `tidy.matchatr_fit` and `summary.matchatr_fit` branch on `inherits(model, "causatr_fit")` and report the marginal contrast (the risk difference by default) instead of a coefficient table.
- Its contrast records `ci_method = "sandwich"` whether `"model"` or `"sandwich"` was asked for, because a marginal g-formula contrast has no separate information-matrix variance.
- `ci_method = "bootstrap"` (`ccw_bootstrap_ci()`) resamples cases and controls separately, so n₁ and n₀ and therefore the known-q₀ weights stay fixed across replicates; it refits per replicate and reports a percentile interval around the analytic point estimate. Resampling the whole sample would mix the strata.
- An estimated prevalence (`unmatched_cc(prevalence = q0, prevalence_n = N)`) widens the interval and leaves the estimate unchanged. The analytic path adds (∂ψ/∂q₀)² · q₀(1 − q₀)/N, with ∂ψ/∂q₀ by a central finite difference on the reported scale (`ccw_estimated_q0_term()`, `ccw_apply_estimated_q0()`); the bootstrap draws q₀* ~ Binomial(N, q₀)/N per replicate. `ccw_boot_point()` must set `fit$design$prevalence_n <- NULL` before computing the contrast, because both variance branches call back into it; without that line they recurse without end.
- Case-control weighting runs on `unmatched_cc()` and `matched_cc()`, never on `nested_cc()`, whose risk-set controls are not a case-control sample. `matcha()` rejects that pairing with `matchatr_bad_estimator`, pointing to `ipw_cox`, before checking the prevalence. On a matched design the matching variable is a baseline covariate that must be in `confounders`, and the matched sets are ignored (Rose and van der Laan 2009); the only test of case-control weighting on a matched design checks the point estimate of `ccw_gformula`, `ccw_aipw` and `ccw_tmle` against a known risk difference (`test-ccw.R`), with no `ccw_ipw`, SE, coverage or oracle check, and the method has not been checked against the paper.
- All four estimators drop incomplete rows once, in `ccw_prepare()`, with a `matchatr_dropped_rows` warning, and compute `cc_weights()` on the complete cases so the weighted case fraction still equals q₀. Computing weights on the full sample and dropping afterwards breaks the mapping. A non-converging targeting step warns `matchatr_tmle_convergence` and falls back to the initial fit.
- Before reporting a variance problem, check whether q₀ is known or estimated (the estimated case has the extra term), compare sandwich with bootstrap, and compare the marginal contrast with causatr run on the explicitly reweighted pseudo-cohort.
- The double-robustness test checks `expect_lt(abs(est − truth), 0.01)`. The risk difference is about 0.05, where `expect_equal()` switches between relative and absolute tolerance unpredictably; the data-generating process puts consistent estimators within 0.003 and the misspecified one beyond 0.03.
