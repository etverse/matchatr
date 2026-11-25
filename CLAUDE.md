# matchatr

matchatr analyses sampled designs: unmatched and matched case-control, nested case-control, case-cohort, two-phase and counter-matched designs. It provides the classical estimators each design supports (conditional logistic regression, weighted Cox, case-cohort pseudo-likelihood) and marginal causal effects through case-control weighting (g-formula, IPW, AIPW, TMLE) and design-based inclusion weights.

It is part of the etverse, whose design (https://etverse.github.io/design.html) fixes what matchatr owns and the interface every package shares. matchatr does not follow that design yet; until it is brought into line, this file and the package documentation describe what the code does today.

## API

```r
fit <- matcha(data, outcome = "case", exposure = "x",
              design = unmatched_cc(prevalence = 0.02),
              confounders = ~ age + smoke, estimator = "ccw_gformula")
contrast(fit, type = "difference", ci_method = "sandwich")

fit <- matcha(ncc, outcome = "d", exposure = "x",
              design = nested_cc(strata = "set", time = "t"), estimator = "ipw_cox")
contrast(fit)   # hazard ratio
```

The valid design and estimator pairs are in `dispatch_table()` (`R/dispatch.R`). `sample_ncc()` draws a nested case-control sample from a cohort. S3 classes: `matchatr_fit` (from `matcha()`), `matchatr_result` (from `contrast()`) and `matchatr_design` (from the design constructors); all errors carry `matchatr_error`.

## Code layout

Files in `R/` are named by what they hold; `tests/testthat/test-foo.R` tests `R/foo.R`.

- Design and API: `cc_design.R`, `matcha.R`, `dispatch.R`, `contrast.R`, `constructors.R`, `checks.R`, `checks_design.R`, `resolve.R`
- Sampling and weights: `risk_set_sampling.R`, `weights_design.R`, `weights_cc.R`, `multi_endpoint.R`
- Classical estimators: `unconditional.R`, `mantel_haenszel.R`, `clogit.R`, `mcnemar.R`, `polytomous.R`, `effect_modification.R`, `homogeneity.R`, `weighted_cox.R`, `ipw_cox.R`, `aft_ncc.R`, `additive_ncc.R`, `lin_ying.R`, `aalen_cumulative.R`, `excess_risk.R`, `case_cohort.R`, `absolute_risk*.R`, `coef_extract.R`
- Case-control weighting: `ccw_prepare.R`, `ccw.R`, `tmle_ccw.R`, `ccw_tmle_contrast.R`, `variance_ccw.R`
- Design-weighted causal survival: `causal_survival_sampled.R`, `contrast_surv_sampled.R`, `variance_surv_sampled.R`
- S3 methods: `print.R`, `summary.R`, `tidy.R`; package doc in `matchatr-package.R`

Invariants that look like bugs and are not, and the axes a combination audit runs over, are in `.claude/rules/`; each file loads when you open the code it covers.

## Commands

The `Makefile` has `make document`, `make test`, `make check`, `make format` (Air), `make lint` and `make site` (altdoc). While iterating, run one file:

```bash
Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-clogit.R")'
```

The full suite takes about three minutes.

## Planned work

Planned work and known bugs are GitHub issues (`gh issue list`, `gh issue view <n>`), grouped in `ROADMAP.md`; the reason behind shipped code is in the roxygen of the function it explains. The methods follow the *Handbook of Statistical Methods for Case-Control Studies* (Borgan, Breslow, Chatterjee, Gail, Scott and Wild 2018). `FEATURE_COVERAGE_MATRIX.md` records which combinations are tested.

## Dependencies

Imports: boot, causatr, data.table, generics, nnet, numDeriv, rlang, sandwich, survival. matchatr calls causatr's exported `causat()`, `contrast()` and `static()` for case-control weighting. survatr is in Suggests and not called: design-weighted survival uses matchatr's own weighted Cox. Oracles in Suggests include Epi, multipleNCC, timereg, tmle, survey and WeightIt.
