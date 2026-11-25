# matchatr roadmap

matchatr analyses sampled designs: unmatched and matched case-control, nested case-control, case-cohort and two-phase designs, with the classical estimators each design supports and marginal causal effects through case-control and inclusion weights. Its methods follow the *Handbook of Statistical Methods for Case-Control Studies* (Borgan, Breslow, Chatterjee, Gail, Scott and Wild 2018). The [etverse design](https://etverse.github.io/design.html) fixes its scope and interface. Each item's status is in the [issues](https://github.com/etverse/matchatr/issues), not in this file.

## 0.1.0

The first release comes once every etverse estimator has been checked against independent software, the bugs those checks find are fixed, and the package follows the design. Its work is the [0.1.0 milestone](https://github.com/etverse/matchatr/milestone/1):

- Correct the standard error of absolute risk under nested case-control IPW ([#27](https://github.com/etverse/matchatr/issues/27)).
- Fix `two_phase()`, which returns a fit with no model ([#28](https://github.com/etverse/matchatr/issues/28)).
- Check case-control weighting on matched data, including whether the matching variable is adjusted for ([#29](https://github.com/etverse/matchatr/issues/29)).

## Causal survival in sampled designs

- Survival g-computation for nested case-control ([#30](https://github.com/etverse/matchatr/issues/30)), with RMST-ratio and quantile scales ([#31](https://github.com/etverse/matchatr/issues/31)).
- Categorical ([#33](https://github.com/etverse/matchatr/issues/33)) and continuous or shifted exposures ([#34](https://github.com/etverse/matchatr/issues/34)).
- An analytic sandwich ([#35](https://github.com/etverse/matchatr/issues/35)) and a doubly robust estimator ([#36](https://github.com/etverse/matchatr/issues/36)).
- Competing risks and cumulative incidence ([#45](https://github.com/etverse/matchatr/issues/45)).

## Two-phase designs and weights

- The two-phase survey engine ([#32](https://github.com/etverse/matchatr/issues/32)) and calibrated weights ([#37](https://github.com/etverse/matchatr/issues/37)).
- Chen weights and IPW for counter-matched nested case-control ([#46](https://github.com/etverse/matchatr/issues/46)).

## Missing data

- Multiple imputation for a covariate missing by design ([#38](https://github.com/etverse/matchatr/issues/38)), for missing confounders and outcomes under case-control weighting ([#39](https://github.com/etverse/matchatr/issues/39)), and congenial with conditional logistic regression ([#40](https://github.com/etverse/matchatr/issues/40)). These pool through causatr's developer interface.

## Other estimators

- Semiparametric maximum likelihood, once an independent oracle exists ([#41](https://github.com/etverse/matchatr/issues/41)).
- Small-sample and separation methods ([#42](https://github.com/etverse/matchatr/issues/42)), alternative risk models ([#43](https://github.com/etverse/matchatr/issues/43)) and a value-against-value odds ratio for a smooth exposure ([#47](https://github.com/etverse/matchatr/issues/47)).
- Secondary-outcome analysis ([#44](https://github.com/etverse/matchatr/issues/44)).

## Open scope questions

- Time-varying covariates and exposures in sampled designs ([#48](https://github.com/etverse/matchatr/issues/48)).

## Not planned here

- Self-controlled case series and power calculations, which the design puts outside matchatr.
- Longitudinal response-selective sampling, which would belong with dynamatr.
- Cohort survival estimands: survatr.
