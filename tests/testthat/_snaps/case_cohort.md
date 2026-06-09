# case_cohort rejects an invalid method

    Code
      case_cohort(subcohort = "sc", time = "t", method = "BadMethod")
    Condition
      Error in `case_cohort()`:
      ! `method` must be one of: "Prentice", "SelfPrentice", "LinYing", "I.Borgan", "II.Borgan".
      i Default is `"Prentice"`.

# contrast_cch rejects type = or

    Code
      contrast(fit, type = "or")
    Condition
      Error in `contrast()`:
      ! A case-cohort design is reported on the hazard-ratio scale.
      i The Prentice / Self-Prentice / Lin-Ying pseudo-likelihood identifies the hazard ratio. Use `type = "hr"` (the default).

# contrast_cch rejects ci_method = sandwich

    Code
      contrast(fit, ci_method = "sandwich")
    Condition
      Error in `contrast()`:
      ! `ci_method = "sandwich"` is not available for the case-cohort estimator.
      i It reports the pseudo-likelihood asymptotic variance; use `ci_method = "model"` (the default).

# Borgan I without stratum aborts with matchatr_bad_design

    Code
      matcha(nwtco, outcome = "rel", exposure = "histol", design = case_cohort(
        subcohort = "in.subcohort", time = "edrel", method = "I.Borgan", id = "seqno"),
      confounders = ~ stage + age, estimator = "cch")
    Condition
      Error in `fit_cch()`:
      ! `method = "I.Borgan"` requires a subcohort stratification column.
      i Supply `stratum = "<column>"` in `case_cohort()` to specify which column defines the subcohort sampling strata.

