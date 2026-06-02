# Phase 15 — Small-Sample Methods: Firth, Exact, and Penalized Likelihood

> **Status: DESIGN**
> Book chapters: 8 (Small Sample Methods).

## Scope

**In:** estimation of log ORs under sparse data / separation, where ordinary (conditional)
MLE is biased or nonexistent — Firth penalized likelihood (unmatched and matched/
conditional), exact and Monte-Carlo / MCMC conditional logistic, and log-F penalties
(data augmentation). Both unmatched (Phase 2) and matched (Phase 3) settings.

**Out:** the standard large-sample fits (Phases 2–3); Bayesian full-model inference.

## Key design decisions

- **Firth penalized likelihood** (Jeffreys-prior penalty) is the recommended default for
  separation — `logistf::logistf()` (unmatched) and `logistf::clogistf()` / `coxphf`
  (matched/conditional). Matched Firth conditional is the newest/most-recommended option.
- **For n ≤ ~100, log-F(2,2) penalty via data augmentation** can outperform Firth (lowest
  bias/MSE; Greenland-Mansournia). Offer as an option.
- **Exact conditional / mid-p** for the smallest tables (network enumeration / MCMC via
  `elrm`); exact p-values are conservative (discreteness) — provide the mid-p correction.
- **Use penalized-likelihood-ratio (not Wald) CIs** for penalized fits.
- These integrate as a `small_sample =` modifier on Phase 2/3 fits rather than a separate
  estimator, since they are alternative *fitting* methods for the same model.

## API design

```r
matcha(data, outcome = "case", exposure = "x", design = matched_cc(strata = "set"),
       estimator = "clogit", small_sample = "firth")          # logistf::clogistf
matcha(data, outcome = "case", exposure = "x", design = unmatched_cc(),
       estimator = "logistic", small_sample = "logF")         # data augmentation
matcha(..., small_sample = "exact")                            # elrm / mid-p
```

## Support matrix

| Setting | Method | Engine | Inference | Status |
|---|---|---|---|---|
| unmatched, separation | firth | logistf::logistf | penalized-LR CI | needs-test |
| matched, sparse strata | firth (conditional) | logistf::clogistf / coxphf | penalized-LR | needs-test |
| n ≤ 100 | log-F(2,2) | data augmentation + glm | LR | needs-test |
| tiny tables | exact / mid-p | elrm | exact / mid-p | smoke |

## Implementation plan

- `R/small_sample.R` — dispatch `small_sample` ∈ {firth, logF, exact} to the right
  engine; build the augmented dataset for log-F; route penalized-LR CIs.

## Variance / inference notes

Penalized fits: use penalized-likelihood-ratio CIs (Wald is unreliable under separation).
Exact: conservative; offer mid-p. Document that point estimates are bias-reduced, not
unbiased.

## Oracle testing strategy

- `logistf` (Firth, unmatched + conditional), `elrm` (MCMC exact) as engines + oracles.
- Truth-based / separation: construct a separated dataset where MLE diverges; confirm
  Firth gives finite, sensibly-shrunk estimates and the LR CI has correct coverage.
- log-F via augmentation vs `logistf` agreement on small n.

## Chunk plan

1. Firth (unmatched + conditional) via `logistf` + separation test + penalized-LR CI.
2. log-F(2,2) data augmentation.
3. Exact / mid-p via `elrm` (smoke).

## Deferred items

Saddlepoint approximations; full Bayesian models; exact methods for large strata
(intractable).
