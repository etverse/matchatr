---
paths:
  - "tests/testthat/**"
---

# Tests

- Run a file with `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-<name>.R")'`.
- Truth tests simulate a cohort with a known effect and then sample a case-control, nested case-control or case-cohort design from it; the data-generating processes are in `helper-dgp.R`.
- Python oracles for the classical estimators (`statsmodels`) are committed data and result files under `fixtures/python/`, compared in `test-python-oracle.R`; `delicatessen` is reserved for the causal and sandwich estimands.
- R oracles: `survival::clogit` and `survival::cch` (classical), `Epi::ccwc` (risk-set sampling), `multipleNCC` (nested case-control IPW), `timereg::aalen` (additive hazards), `tmle` (the targeting step), and causatr on the explicitly reweighted pseudo-cohort (case-control weighting).
