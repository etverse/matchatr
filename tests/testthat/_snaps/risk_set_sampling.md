# the empty-risk-set and collision errors read clearly

    Code
      sample_ncc(bad, time = "t", event = "d", m = 2L)
    Condition
      Error in `sample_ncc()`:
      ! 1 case(s) had no eligible control at their failure time, so a risk set could not be formed.
      i Affected failure time(s): 3.
      i Check the `time` origin/scale, the `entry` column, or whether `match` strata are too fine to leave any at-risk control.

---

    Code
      sample_ncc(clash, time = "t", event = "d", m = 1L)
    Condition
      Error in `sample_ncc()`:
      ! `cohort` already has column(s) `case`, which `sample_ncc()` appends to the output.
      i Rename them so the sampled set id / case indicator / risk time are unambiguous.

