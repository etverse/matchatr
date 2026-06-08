# print.matchatr_design renders each design

    Code
      print(unmatched_cc(prevalence = 0.02))
    Output
      <matchatr_design>
       Type:       Unmatched case-control
       Prevalence: 0.02 (q0)
       Weights:    case-control (q0)

---

    Code
      print(matched_cc(strata = c("age_grp", "sex"), ratio = 2))
    Output
      <matchatr_design>
       Type:       Matched case-control
       Strata:     age_grp, sex
       Ratio:      2:1
       Weights:    none

---

    Code
      print(nested_cc(strata = "set", time = "t", ratio = 3))
    Output
      <matchatr_design>
       Type:       Nested case-control
       Strata:     set
       Time:       t
       Ratio:      3:1
       Weights:    inclusion-probability

---

    Code
      print(case_cohort(subcohort = "in_subcohort", time = "t"))
    Output
      <matchatr_design>
       Type:       Case-cohort
       Time:       t
       Subcohort:  in_subcohort
       Method:     Prentice
       Weights:    inclusion-probability

---

    Code
      print(two_phase(phase1 = "stratum", phase2 = "in_phase2"))
    Output
      <matchatr_design>
       Type:       Two-phase
       Phase 1:    stratum
       Phase 2:    in_phase2
       Weights:    design (two-phase)

---

    Code
      print(counter_matched(strata = "surrogate", time = "t"))
    Output
      <matchatr_design>
       Type:       Counter-matched
       Strata:     surrogate
       Time:       t
       Weights:    counter-matching

# print.matchatr_fit renders the resolved analysis

    Code
      print(fit_adj)
    Output
      <matchatr_fit>
       Design:     Unmatched case-control
       Estimator:  logistic  (engine: glm_logistic)
       Outcome:    case
       Exposure:   x
       Confounders: ~age + smoke
       N:          60  (cases: 20, controls: 40)

---

    Code
      print(fit_clogit)
    Output
      <matchatr_fit>
       Design:     Matched case-control
       Estimator:  clogit  (engine: clogit)
       Outcome:    case
       Exposure:   x
       Confounders: none
       N:          60  (cases: 20, controls: 40)

---

    Code
      print(fit_ccw)
    Output
      <matchatr_fit>
       Design:     Unmatched case-control
       Estimator:  ccw_gformula  (engine: ccw_gformula)
       Outcome:    case
       Exposure:   x
       Confounders: ~age
       N:          60  (cases: 20, controls: 40)

