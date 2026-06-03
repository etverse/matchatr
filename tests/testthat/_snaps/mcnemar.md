# McNemar rejection messages read clearly

    Code
      matcha(df3, "case", "x", matched_cc(strata = "set"), estimator = "mcnemar")
    Condition
      Error in `matcha()`:
      ! The McNemar estimator requires 1:1 matched pairs, but 80 matched set(s) have more than one case or more than one control.
      i Use `estimator = "clogit"` for M:1, variable-ratio, or richer matching.

---

    Code
      contrast(fit, type = "or", ci_method = "bootstrap")
    Condition
      Error in `contrast()`:
      ! `ci_method = "bootstrap"` is not available for the McNemar estimator.
      i It reports the McNemar interval (Var(log OR) = 1/n10 + 1/n01); use `ci_method = "model"` (the default).

