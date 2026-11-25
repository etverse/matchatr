# matchatr (development version)

No version has been released yet. What the package does today:

* `matcha()` fits a case-control analysis and `contrast()` reports the effect, for unmatched, matched, nested case-control, case-cohort and counter-matched designs declared with `unmatched_cc()`, `matched_cc()`, `nested_cc()`, `case_cohort()` and `counter_matched()`. `tidy()` and `summary()` give the result as a table.
* Unmatched and matched case-control data give odds ratios: logistic regression, Mantel-Haenszel, conditional logistic regression (`clogit`), McNemar for 1:1 pairs, stratum-specific effects through `effect_modifier`, and polytomous odds ratios for several case groups with `test_homogeneity()`.
* Nested case-control data give the risk-set hazard ratio. `sample_ncc()` draws a risk-set sample from a cohort, with Samuelsen inclusion weights on request. `ipw_cox`, `ipw_aft` and `ipw_aalen` fit weighted Cox, accelerated failure time and additive-hazards models, `absolute_risk()` and `excess_risk()` report absolute and excess risk, and `reuse_ncc_endpoint()` reuses one control set for a second endpoint.
* Case-cohort data give Prentice, Self-Prentice, Lin-Ying and Borgan hazard ratios through `survival::cch`, with absolute risk; `surv_gcomp` reports marginal risk differences, risk ratios and restricted mean survival time differences by g-computation on the design's own weighted Cox model.
* `ccw_gformula`, `ccw_ipw`, `ccw_aipw` and `ccw_tmle` report marginal risk differences, risk ratios and odds ratios from an unmatched or matched case-control sample with a known or estimated prevalence, using case-control weights; intervals come from the sandwich or influence-function variance or from a within-stratum bootstrap.
* `two_phase()` declares a two-phase design, but no estimator is wired to it yet.
