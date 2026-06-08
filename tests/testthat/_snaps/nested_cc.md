# NCC contrast rejection and print messages read clearly

    Code
      contrast(fit, type = "or")
    Condition
      Error in `contrast()`:
      ! A nested case-control design is reported on the hazard-ratio scale.
      i Risk-set (incidence-density) sampling identifies the hazard ratio (OR = HR exactly; Prentice & Breslow 1978). Use `type = "hr"` (the default).

---

    Code
      contrast(matched, type = "hr")
    Condition
      Error in `contrast()`:
      ! A matched case-control design does not identify a hazard ratio.
      i The hazard ratio needs risk-set (incidence-density) sampling -- a nested case-control design (`nested_cc()`). Report the conditional odds ratio with `type = "or"` (the default).

---

    Code
      print(contrast(fit))
    Output
      <matchatr_result>
       Estimator:  clogit  (engine: clogit)
       Estimand:   hazard ratio
       Contrast:   Hazard ratio
       CI method:  model
       N:          693
      
      Contrasts:
         comparison estimate        se ci_lower ci_upper
             <char>    <num>     <num>    <num>    <num>
      1:          x 1.882519 0.3047074 1.370764 2.585331

