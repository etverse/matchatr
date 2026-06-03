# logistic-contrast rejections read clearly

    Code
      contrast(fit, type = "difference")
    Condition
      Error in `contrast()`:
      ! The risk difference is not identified from an unmatched case-control sample without the source-population prevalence q0.
      i Report the conditional odds ratio with `type = "or"`.
      i For a marginal risk difference / ratio, supply `prevalence =` on the design and use a case-control-weighted estimator (e.g. `estimator = "ccw_gformula"`).

---

    Code
      contrast(fit, type = "ratio")
    Condition
      Error in `contrast()`:
      ! The risk ratio is not identified from an unmatched case-control sample without the source-population prevalence q0.
      i Report the conditional odds ratio with `type = "or"`.
      i For a marginal risk difference / ratio, supply `prevalence =` on the design and use a case-control-weighted estimator (e.g. `estimator = "ccw_gformula"`).

---

    Code
      contrast(fit, type = "or", ci_method = "bootstrap")
    Condition
      Error in `contrast()`:
      ! Bootstrap confidence intervals are not provided for the conditional odds ratio.
      i Use `ci_method = "model"` (Wald) or `ci_method = "sandwich"` (robust).

---

    Code
      matcha(dord, "case", "x", unmatched_cc())
    Condition
      Error in `matcha()`:
      ! Exposure `x` is an ordered factor, which is fit with polynomial contrasts (not per-level odds ratios).
      i Pass a numeric score for a trend OR, or an unordered factor for per-level ORs.

