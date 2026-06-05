# design-constructor rejections read clearly

    Code
      unmatched_cc(prevalence = 1.5)
    Condition
      Error in `unmatched_cc()`:
      ! `prevalence` (q0) must be a single number strictly between 0 and 1.

---

    Code
      matched_cc(strata = "set", ratio = 1.5)
    Condition
      Error in `matched_cc()`:
      ! `ratio` must be a single whole number >= 1 (controls per case).

---

    Code
      matched_cc(strata = character(0))
    Condition
      Error in `matched_cc()`:
      ! `strata` must be a non-empty character vector of column names.

---

    Code
      nested_cc(strata = "set", time = 5)
    Condition
      Error in `nested_cc()`:
      ! `time` must be a single non-empty character string.

# matcha rejections read clearly

    Code
      matcha(df, "case", "x", unmatched_cc(), estimator = "bogus")
    Condition
      Error in `matcha()`:
      ! Estimator `bogus` is not available for design `unmatched_cc`.
      i Supported estimators for this design: "logistic", "mh", "polytomous", "ccw_gformula", "ccw_ipw", "ccw_aipw", "ccw_tmle".

---

    Code
      matcha(df, "case", "x", unmatched_cc(), estimator = "ccw_ipw")
    Condition
      Error in `matcha()`:
      ! Estimator `ccw_ipw` needs the source-population prevalence q0 to reweight the sample.
      i Supply it on the design, e.g. `unmatched_cc(prevalence = 0.02)`.

---

    Code
      matcha(df, "case", "missing_col", unmatched_cc())
    Condition
      Error in `matcha()`:
      ! Column(s) for `exposure` not found in `data`: `missing_col`.

---

    Code
      matcha(df, "case", "x", matched_cc(strata = "no_such_set"))
    Condition
      Error in `matcha()`:
      ! Column(s) for `design` not found in `data`: `no_such_set`.

---

    Code
      matcha(df, "case", "case", unmatched_cc())
    Condition
      Error in `matcha()`:
      ! `outcome` and `exposure` must be different columns.

---

    Code
      matcha(bad_y, "case", "x", unmatched_cc())
    Condition
      Error in `matcha()`:
      ! Outcome `case` must be a binary case indicator (logical, two-level factor, or numeric 0/1).
      i Multiple case / control groups are handled by a polytomous estimator.

---

    Code
      matcha(df, "case", "x", "not a design")
    Condition
      Error in `matcha()`:
      ! `design` must be a `matchatr_design` object.
      i Build one with e.g. `unmatched_cc()`, `matched_cc()`, or `nested_cc()`.

# the uninformative-stratum warning reads clearly

    Code
      invisible(matcha(bad, "case", "x", matched_cc(strata = "set")))
    Condition
      Warning in `matcha()`:
      1 stratum/strata have no cases or no controls and carry no information for the conditional likelihood.
      i `survival::clogit` drops these sets; check the matched/risk-set sampling.

