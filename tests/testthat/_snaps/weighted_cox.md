# the empty-opposite-stratum error reads clearly

    Code
      sample_ncc_counter_matched(bad, time = "t", event = "d", surrogate = "z_bin")
    Condition
      Error in `sample_ncc_counter_matched()`:
      ! 1 case(s) had no eligible control at their failure time, so a risk set could not be formed.
      i Affected failure time(s): 3.
      i Check the `time` origin/scale, the `entry` column, or whether `match` strata are too fine to leave any at-risk control.

# the missing weights column error reads clearly

    Code
      matcha(ncc, "case", "x", counter_matched(strata = "set", time = "risk_time"),
      estimator = "weighted_cox")
    Condition
      Error in `fit_weighted_cox()`:
      ! A counter-matched design requires a `weights` column.
      i Supply the name of the log-weight column via `counter_matched(weights = "log_w")`. `sample_ncc_counter_matched()` appends it automatically.

# the OR-from-counter-matched error reads clearly

    Code
      contrast(fit, type = "or")
    Condition
      Error in `contrast()`:
      ! A counter-matched design is reported on the hazard-ratio scale.
      i The counter-matched partial likelihood identifies the hazard ratio (Langholz & Borgan 1995). Use `type = "hr"` (the default).

