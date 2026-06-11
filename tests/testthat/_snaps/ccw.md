# ccw_gformula rejects a non-binary exposure

    Code
      matcha(cc, outcome = "case", exposure = "xc", design = unmatched_cc(prevalence = attr(
        cc, "q0")), confounders = ~w, estimator = "ccw_gformula")
    Condition
      Error in `fit_ccw()`:
      ! The CCW g-formula estimator requires a binary exposure; `xc` is not binary (logical, two-level factor, or numeric 0/1).
      i For a categorical (k>2) or continuous exposure use `a conditional estimator (e.g. estimator = "logistic")`.

# ccw_gformula requires confounders to standardize over

    Code
      matcha(cc, outcome = "case", exposure = "x", design = unmatched_cc(prevalence = attr(
        cc, "q0")), estimator = "ccw_gformula")
    Condition
      Error in `fit_ccw()`:
      ! The CCW g-formula requires `confounders` to standardize over.
      i Supply an adjustment set, e.g. `confounders = ~ age + smoke`, on `matcha()`.

# ccw_gformula rejects bootstrap variance and off-scale contrasts

    Code
      contrast(fit, type = "difference", ci_method = "bootstrap")
    Condition
      Error in `contrast()`:
      ! `ci_method = "bootstrap"` is not available for the CCW g-formula.
      i Use `ci_method = "model"` or `ci_method = "sandwich"` (causatr's influence-function variance on the weighted fit).

---

    Code
      contrast(fit, type = "hr")
    Condition
      Error in `contrast()`:
      ! The CCW g-formula reports a marginal effect, not `type = "hr"`.
      i Use `type = "difference"` (risk difference), `"ratio"` (risk ratio), or `"or"` (marginal odds ratio).

