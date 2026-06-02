
<!-- README.md is generated from README.Rmd. Please edit that file -->

# matchatr

<!-- badges: start -->

[![R-CMD-check](https://github.com/etverse/matchatr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/etverse/matchatr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**matchatr** provides causal inference for **(matched) case-control,
nested case-control (NCC), and case-cohort** study designs. It pairs
design-faithful classical estimators with marginal causal effects, and
integrates with the [etverse](https://github.com/etverse) ecosystem —
delegating estimation to [`causatr`](https://github.com/etverse/causatr)
(g-computation / IPW / AIPW with sandwich and bootstrap variance) and
[`survatr`](https://github.com/etverse/survatr) (causal survival on
person-period data).

> **Status: early development.** The package scaffold and the full
> `PHASE_*.md` design roadmap are in place; estimator implementations
> are being built phase by phase. The API below is the planned interface
> (fixed in `PHASE_1`).

## What it does

Two orthogonal axes: a **`design`** object encodes the sampling
structure (strata, matching ratio, time scale, prevalence, inclusion
weights); an **`estimator`** chooses the analysis.

| Design | Classical estimand | Causal (marginal) estimand |
|----|----|----|
| Unmatched case-control | conditional OR, Mantel-Haenszel | RD / RR / marginal OR (case-control weighting) |
| Matched case-control | conditional OR (conditional logistic) | RD / RR via standardization |
| Nested case-control | risk-set HR; Samuelsen IPW Cox | marginal survival contrasts (design-weighted) |
| Case-cohort | Prentice / Self-Prentice / Borgan HR | absolute risk, RD(t), RMST |

Marginal causal effects use **case-control weighting** (the Rose & van
der Laan g-formula / IPW / AIPW / TMLE family) and **design-based
inclusion weighting** (Samuelsen, Borgan): the weights are passed as
observation weights into the etverse engines, so they compose directly
with existing estimators.

## Installation

You can install the development version of matchatr from
[GitHub](https://github.com/etverse/matchatr) with:

``` r
# install.packages("pak")
pak::pak("etverse/matchatr")
```

## Example

``` r
library(matchatr)

# Matched case-control -> conditional odds ratio
fit <- matcha(
  data,
  outcome = "case", exposure = "x",
  design = matched_cc(strata = "set"),
  confounders = ~ age + smoke, estimator = "clogit"
)

# Marginal causal risk difference from an unmatched case-control sample
fit <- matcha(
  data,
  outcome = "case", exposure = "x",
  design = unmatched_cc(prevalence = 0.02),   # source-population prevalence q0
  confounders = ~ age + smoke, estimator = "ccw_gformula"
)
contrast(fit, type = "difference", ci_method = "sandwich")
```

## Roadmap

The design is documented in `PHASE_1`–`PHASE_20` at the repository root,
mapping the *Handbook of Statistical Methods for Case-Control Studies*
(Borgan et al., 2018) to an implementation plan. See `CLAUDE.md` for the
phase index and `FEATURE_COVERAGE_MATRIX.md` for what is implemented and
tested.

## Part of the etverse

matchatr is one package in the [etverse](https://github.com/etverse)
family for causal inference and methodological triangulation, alongside
`causatr` (causal effect estimation) and `survatr` (causal survival
analysis).
