---
paths:
  - "R/matcha.R"
  - "R/dispatch.R"
  - "R/cc_design.R"
  - "R/checks*.R"
  - "R/resolve.R"
  - "R/constructors.R"
  - "R/weights_*.R"
  - "R/coef_extract.R"
  - "tests/testthat/test-matcha.R"
  - "tests/testthat/test-cc_design.R"
  - "tests/testthat/test-dispatch.R"
  - "tests/testthat/test-rejections.R"
  - "tests/testthat/test-weights_cc.R"
---

# Designs, weights and input checks

These are deliberate. Do not report them as bugs without a numerical reproducer.

- matchatr owns the sampling-design and weighting layer. Estimation goes to causatr (g-computation, IPW, AIPW and their variance), `survival` (`clogit`, `coxph`, `cch`), `nnet` (polytomous) and `survey` (two-phase); do not reimplement them here.
- `estimator =` chooses the analysis and `design =` the sampling structure. The design object carries strata, time, the prevalence q₀ and inclusion weights.
- Case-control weights and design weights are different objects with different variance. Case-control weights (Rose and van der Laan) map a case-control sample back to the source population through the marginal prevalence q₀; design weights (Samuelsen, Borgan) are inverse probabilities of being sampled as a control or subcohort member. Never share a code path that assumes one variance form for both.
- Weights enter a fit as the `weights` argument, carried by the design object, never as a column of the data.
- Rejections with classed errors: a non-binary case-control outcome (`matchatr_bad_outcome`, the same check for logical, two-level factor and 0/1 encodings, and requiring both a case and a control); a case-control-weighting estimator without a prevalence (`matchatr_missing_prevalence`); an outcome or exposure that is also a confounder or a design column (`matchatr_bad_input`). Confounders may overlap design columns, since adjusting for a frequency-matched variable is valid.
- A stratum with no cases or no controls is uninformative; it warns `matchatr_uninformative_stratum`.
- Never index an SE vector by position against `coef()`: `stats::vcov()` keeps aliased coefficients as `NA` rows and `sandwich::sandwich()` drops them. Reduce to the estimable coefficients and match by name (`estimable_vcov()`); a non-estimable exposure raises `matchatr_unestimable_exposure`.
- Odds-ratio intervals are Wald on the log scale, exponentiated. The `se` in `contrasts` is the delta-method `OR · SE(log OR)` and does not reconstruct the interval; the log-scale SE is in `estimates`.
- Known bug: `two_phase()` is exported, but `matcha()` on it returns a fit with `model = NULL`, and the user finds out only at `contrast()`.
- For a `clogit` or `coxph` fit, `model$n` is the number of rows used and `nobs()` the number of events; the analysis size and `n_dropped = nrow(data) − model$n` read `model$n`. A stratum `clogit` drops as uninformative is still counted in `model$n`, so it raises no `matchatr_dropped_rows` warning. For a `survival::cch` fit `model$n` double-counts subcohort cases, so no drop count can be computed from it; `NA` handling there is left to `cch`'s `na.action`.
