---
paths:
  - "R/clogit.R"
  - "R/mcnemar.R"
  - "R/unconditional.R"
  - "R/mantel_haenszel.R"
  - "R/polytomous.R"
  - "R/effect_modification.R"
  - "R/homogeneity.R"
  - "R/contrast.R"
  - "tests/testthat/test-clogit.R"
  - "tests/testthat/test-mcnemar.R"
  - "tests/testthat/test-nested_cc.R"
  - "tests/testthat/test-polytomous.R"
  - "tests/testthat/test-unconditional.R"
---

# Matched and nested case-control: conditional likelihood

These are deliberate. Do not report them as bugs without a numerical reproducer.

- Matched case-control uses the conditional likelihood (`survival::clogit`), never unconditional logistic regression on matched-set indicators, which biases the odds ratio (towards OR² for 1:1 matching; Pike et al. 1980, Breslow and Day 1980).
- The classical nested case-control analysis is the same conditional likelihood with each sampled risk set (a case and its controls) as a stratum (Thomas 1977); do not drop the strata.
- Under risk-set (incidence-density) sampling the odds ratio equals the hazard ratio with no rare-disease assumption (Miettinen 1976; Prentice and Breslow 1978).
- One scale per design: `clogit` gives `exp(β)`, which a matched design reports as an odds ratio (`type = "or"`) and a nested design as a hazard ratio (`type = "hr"`). `default_contrast_type()` follows the design, and `reject_offdesign_conditional_scale()` rejects the other label with `matchatr_unidentified_estimand`: the number is the same, the estimand is not. The design's `time` column records how controls were sampled; the conditional likelihood reads the risk set from `strata` only.
- `estimator = "mcnemar"` is the 1:1 binary-exposure case only: OR = n₁₀/n₀₁ with Var(log OR) = 1/n₁₀ + 1/n₀₁ over discordant pairs. A set with more than one case or control raises `matchatr_not_one_to_one` and points to `clogit`; an empty side (n₁₀ = 0 or n₀₁ = 0) raises `matchatr_unestimable_exposure`. Exposures are paired by `order()` on the stratum factor, which is correct under any row order and set ids.
- The McNemar closed form is the only oracle for the matched conditional variance that is independent of `clogit`; comparing a wrapper's SE with `clogit`'s own `vcov()` only checks forwarding. A recovery test against a truth uses a band scaled by the SE.
- `nnet::multinom` does not set a rank-deficient coefficient to `NA`: it splits it across collinear columns and returns an attenuated odds ratio. `reject_collinear_exposure()` (`R/polytomous.R`) checks the design-matrix rank and aborts with `matchatr_unestimable_exposure` only when dropping the exposure columns does not lower the rank by their number. Collinearity among confounders alone leaves the exposure identified and is not rejected.
