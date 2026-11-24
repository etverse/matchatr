---
paths:
  - "R/**"
  - "tests/testthat/**"
---

# Combination axes

A combination audit runs over these axes. The design-and-estimator pairs come from `dispatch_table()` in `R/dispatch.R`, plus the case-control-weighting estimators that `matcha()` routes to causatr; whether a combination is supported is decided by the rejection paths in `R/checks*.R`, `R/resolve.R` and the engines, and recorded per combination in `FEATURE_COVERAGE_MATRIX.md`. Re-read the table when an axis may have changed.

| Design (`design =`) | Estimators (`estimator =`) |
|---|---|
| `unmatched_cc()` | `logistic`, `mh`, `polytomous`, `ccw_gformula`, `ccw_ipw`, `ccw_aipw`, `ccw_tmle` |
| `matched_cc()` | `clogit`, `mcnemar`, `ccw_gformula`, `ccw_ipw`, `ccw_aipw`, `ccw_tmle` |
| `nested_cc()` | `clogit`, `ipw_cox`, `ipw_aft`, `ipw_aalen` |
| `case_cohort()` | `cch`, `surv_gcomp` |
| `two_phase()` | `survey` |
| `counter_matched()` | `weighted_cox` |

| Axis | Values | Where |
|---|---|---|
| contrast scale | `difference`, `ratio`, `or`, `hr`, `af`, `excess`, `rmst` | `contrast(type =)` |
| variance | `model`, `sandwich`, `bootstrap` | `contrast(ci_method =)` |
| weights | none, case-control weights from a known or estimated prevalence (`prevalence`, `prevalence_n`), design or inclusion-probability weights | the design constructors in `R/cc_design.R`, `R/weights_*.R` |
| outcome | binary (case-control), time to event (nested case-control, case-cohort), several case groups (`polytomous`) | `matcha()` |
