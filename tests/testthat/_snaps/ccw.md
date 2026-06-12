# ccw_gformula rejects a non-binary exposure

    Code
      matcha(cc, outcome = "case", exposure = "xc", design = unmatched_cc(prevalence = attr(
        cc, "q0")), confounders = ~w, estimator = "ccw_gformula")
    Condition
      Error in `ccw_prepare()`:
      ! The CCW g-formula estimator requires a binary exposure; `xc` is not binary (logical, two-level factor, or numeric 0/1).
      i For a categorical (k>2) or continuous exposure use `a conditional estimator (e.g. estimator = "logistic")`.

# ccw_gformula requires confounders to standardize over

    Code
      matcha(cc, outcome = "case", exposure = "x", design = unmatched_cc(prevalence = attr(
        cc, "q0")), estimator = "ccw_gformula")
    Condition
      Error in `ccw_prepare()`:
      ! The case-control-weighted estimators require `confounders` for the adjustment model(s).
      i Supply an adjustment set, e.g. `confounders = ~ age + smoke`, on `matcha()`.

# ccw_gformula rejects an off-scale contrast

    Code
      contrast(fit, type = "hr")
    Condition
      Error in `contrast()`:
      ! A case-control-weighted estimator reports a marginal effect, not `type = "hr"`.
      i Use `type = "difference"` (risk difference), `"ratio"` (risk ratio), or `"or"` (marginal odds ratio).

