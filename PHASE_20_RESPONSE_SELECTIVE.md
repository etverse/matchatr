# Phase 20 — Response-Selective (Outcome-Dependent) Sampling from Longitudinal Cohorts

> **Status: DESIGN (extension / lower priority)**
> Book chapters: 15 (Response Selective Study Designs).

## Scope

**In:** outcome-dependent sampling from a longitudinal cohort — Phase-1 follows all
subjects, Phase-2 samples on outcome summaries (intercept, slope, or both); valid
estimation of fixed-effect coefficients in linear mixed-effects / GEE models via the
ascertainment-corrected likelihood (ACL), weighted likelihood, full likelihood, or
multiple imputation. Quantitative and binary longitudinal outcomes (the latter via
sequential offsetted regression, SOR).

**Out:** the cross-sectional / survival case-control designs (Phases 1–18). This is a
longitudinal extension that overlaps with `survatr`'s person-period world.

## Key design decisions

- **Known sampling fractions π(q)** on outcome summaries q (by design). Four routes
  (Ch15): (1) **ACL** — score equations with ascertainment correction integrating over
  strata; (2) **weighted likelihood** — inverse-π weighting (inefficient); (3) **full
  likelihood** — use sampled + unsampled subjects (information from the unsampled on the
  population mixture); (4) **multiple imputation** of the expensive covariate (nearly
  full-cohort efficiency even for inefficient designs).
- **Binary longitudinal outcomes**: sequential offsetted regression (SOR) for GEE with a
  two-stage auxiliary sampling-ratio model (Rathouz & Schildcrout) — genuinely new
  machinery.
- **Reuse where possible**: weighted GEE via `survey::svydesign` + `geepack::geeglm`;
  MI via `mice`/`smcfcs`; ACL/SOR/full-likelihood are new. Consider housing the
  person-period plumbing on top of `survatr`.
- **Sampling level matters**: observation-level vs subject-level sampling affects the
  GEE working-independence assumption — enforce/validate.

## API design

```r
matcha(longdata, outcome = "y", exposure = "x",
       design = response_selective(id = "id", time = "t",
                                   sample_on = "slope", fractions = pi_spec),
       estimator = "acl")          # or "weighted" / "full" / "mi" / "sor" (binary)
```

## Support matrix

| Outcome | Route | Engine | Variance | Status |
|---|---|---|---|---|
| quantitative | ACL | custom | info matrix | smoke |
| quantitative | weighted | survey + lme | sandwich | needs-test |
| quantitative | full likelihood | custom | info matrix | smoke |
| quantitative/binary | MI | mice/smcfcs | Rubin | needs-test |
| binary longitudinal | SOR (GEE) | custom + geepack | sandwich | smoke |

## Implementation plan

- `R/response_selective.R` — `response_selective()` design; weighted GEE route
  (`survey` + `geepack`); MI route (`mice`/`smcfcs`); ACL / full-likelihood / SOR
  custom fitters (optional, smoke first).

## Variance / inference notes

Sandwich for weighted/SOR; information matrix for ACL/full; Rubin for MI. Full likelihood
can tighten confidence regions vs weighting; MI nearly matches full-cohort efficiency even
for inefficient designs (CAMP example, Ch15). Watch the observation- vs subject-level
sampling assumption for GEE.

## Oracle testing strategy

- `survey` + `geepack` (weighted), `mice`/`smcfcs` (MI) as engines + oracles. Truth-based:
  simulate a longitudinal cohort, sample on slope, confirm the naive analysis is biased
  and ACL/weighted/full/MI recover the fixed-effect truth.

## Chunk plan

1. Weighted GEE route + MI route + biased-naive-vs-corrected truth test.
2. ACL / full-likelihood (quantitative) — smoke.
3. SOR for binary longitudinal — smoke.

## Deferred items

Time-varying-effect sampling subtleties; integration with survatr's longitudinal causal
machinery; this is an optional extension off the core case-control path.
