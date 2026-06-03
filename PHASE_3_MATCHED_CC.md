# Phase 3 — Matched Case-Control: Conditional Logistic Regression

> **Status: DESIGN**
> Book chapters: 4 (Matched Case-Control Studies).

## Scope

**In:** individually and frequency-matched case-control analysis via conditional
maximum likelihood (CMLE) — `survival::clogit`. McNemar closed form for 1:1 binary
exposure; M:1 and variable matching ratios; adjustment for non-matching covariates and
residual (imperfect-match) differences; effect modification across matching strata.

**Out:** marginal effects from matched data (needs external population distribution of
matching variables — Phase 9 standardization), risk-set/incidence-density matching for
time-to-event (that is NCC — Phase 5).

## Key design decisions

- **Conditional likelihood, never unconditional MLE on matched-set indicators.** For
  1:1 matching the unconditional MLE → OR² in large samples (Pike et al. 1980; Breslow
  & Day 1980). matchatr fits via `survival::clogit` (= Cox partial likelihood with each
  matched set as a stratum). The conditional likelihood for a stratum with one case and
  M controls is ∏ exp(x_case·β) / Σ_{j} exp(x_j·β); for 1:1 it reduces to
  expit{(x_case − x_control)·β}. This is the central correctness invariant
  (`hard-rules.md`).
- **Matching-variable coefficients are not estimable** under CMLE (conditioned away).
  Document that only the exposure/adjustment ORs are reported; matching variables are
  controlled implicitly.
- **Covariate adjustment vs matched-set adjustment.** When matching is on a few
  observable confounders, ordinary covariate adjustment in logistic regression is also
  valid and may be offered (Pearce 2016); but the safe default for matched designs is
  CMLE. The design object's `matched_cc()` selects `estimator = "clogit"` by default.
- **Strata with 0 cases or 0 controls are uninformative** — `clogit` drops them; emit a
  classed warning `matchatr_uninformative_stratum` with the count dropped.
- **Bias from matching on exposure-related variables cannot be removed** without
  external data — documented as a design caveat, not a runtime check.

## API design

```r
fit <- matcha(infert, outcome = "case", exposure = "induced",
              design = matched_cc(strata = "stratum"),
              confounders = ~ spontaneous, estimator = "clogit")
tidy(fit, exponentiate = TRUE)   # conditional OR + CI

# 1:1 McNemar closed form
matcha(pairs, outcome = "case", exposure = "x",
       design = matched_cc(strata = "pair", ratio = 1), estimator = "mcnemar")
```

## Support matrix

| Matching | Estimator | Estimand | Contrast | Variance | Status |
|---|---|---|---|---|---|
| 1:1, binary x | mcnemar | cond. OR | OR | McNemar (1/m10+1/m01) | needs-test |
| 1:1 / M:1 | clogit | cond. OR | OR | partial-lik info | needs-test |
| M:1 + covariates | clogit | cond. OR | OR | partial-lik info | needs-test |
| variable ratio | clogit | cond. OR | OR | partial-lik info | needs-test |
| with `x:modifier` | clogit | stratum-specific OR | OR | partial-lik info | needs-test |
| clogit | — | RD/RR/marginal | — | — | ⛔ `matchatr_unidentified_estimand` |

## Implementation plan

- `R/clogit.R` — `fit_clogit()` (wraps `survival::clogit`, building the
  `case ~ x + cov + strata(set)` formula), McNemar closed form `fit_mcnemar()`,
  stratum-validity check, uninformative-stratum warning.
- S3: `tidy`/`summary`/`print` (conditional OR table). Reuse Phase 2's OR contrast
  assembly on the log scale.

## Variance / inference notes

CMLE: the inverse partial-likelihood information matrix (what `clogit` returns). Wald
CIs on the log-OR scale, exponentiated. McNemar: Var(log OR) = 1/m10 + 1/m01. No
sandwich needed for the conditional fit; cluster-robust is available via `clogit`'s
`cluster()` if controls are reused (relevant when bridging to NCC reuse, Phase 7).

## Oracle testing strategy

- `survival::clogit` is the engine; assert our wrapper reproduces it exactly on the
  `infert` dataset (handbook §4.4 induced-abortion example, CMLE matched-set ORs:
  OR(1+) ≈ 4.0, OR(2+) ≈ 16.8).
- Truth-based: simulate a cohort with known conditional OR, draw 1:M matched sets, check
  CMLE recovers the OR and that unconditional MLE shows the OR² bias for 1:1 (a
  demonstration test pinning the invariant).
- McNemar vs `stats::mcnemar.test` / `exact2x2`.

## Chunk plan

1. `fit_clogit()` + clogit oracle on `infert`. **(done)** — `R/clogit.R`
   wraps `survival::clogit` (`outcome ~ exposure + confounders + strata(set)`),
   reuses the shared `conditional_or_result()` OR assembly, and reports the
   conditional OR with the partial-likelihood Wald interval. Validated against
   `survival::clogit` (exact pass-through), the handbook §4.4 induced-abortion
   ORs, and a matched-set DGP built from the conditional likelihood.
2. McNemar closed form + truth-based OR²-bias demonstration.
3. Effect modification across strata + variable-ratio handling.

## Deferred items

Marginal effects from matched data via standardization with external matching-variable
distribution (Phase 9). Counter-matching (Phase 5). Conditional exact / Firth for
sparse strata (Phase 15).
