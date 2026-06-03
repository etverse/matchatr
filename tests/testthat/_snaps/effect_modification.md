# effect-modification rejection messages read clearly

    Code
      matcha(df, "case", "x", matched_cc(strata = "set"), effect_modifier = "mnum",
      estimator = "clogit")
    Condition
      Error in `matcha()`:
      ! `effect_modifier` `mnum` must be categorical (logical, character, or factor).
      i Bin a continuous modifier or wrap it in `factor()` to report one odds ratio per level.

