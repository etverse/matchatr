# Phase 10 — Design-Weighted Causal Survival for Sampled Cohorts

> **Status: DESIGN**
> Book chapters: 16, 17, 19 (sampling + absolute risk); implements Track 3 of
> `PHASE_8_CAUSAL_STRATEGY` via `survatr`.

## Scope

**In:** marginal causal survival estimands (absolute risk F_x(t), risk difference,
risk ratio, RMST difference) under NCC and case-cohort sampling, by feeding the
design / inclusion-probability weights (Phase 7 Samuelsen, Phase 6 Borgan) into
`survatr`'s weighted person-period causal-survival engine.

**Out:** classical HR estimation (Phases 5–7), the non-survival CCW family (Phase 9),
weight calibration (Phase 12).

## Key design decisions

- **Delegate causal survival to `survatr`.** survatr already does pooled-logistic /
  ICE g-computation, IPW, and AIPW on person-period data with sandwich/bootstrap variance
  and risk-difference / risk-ratio / RMST contrasts. matchatr's job is to (a) convert the
  sampled design to person-period form and (b) supply the inclusion-probability weights
  as survatr observation weights.
- **Design weights, not case-control weights.** This track uses Samuelsen/Borgan
  inverse-inclusion-probability weights (Phase 6/7), which reweight the sampled cohort to
  the full cohort. (Distinct from Phase 9's q₀ weights — `hard-rules.md`.)
- **Variance is the harder part.** survatr's sandwich treats weights as fixed; the
  sampling variation of estimated inclusion weights needs the Samuelsen/Borgan
  correction or a design-aware bootstrap (resample the cohort, re-sample the design,
  refit). Default to bootstrap for marginal contrasts under sampling; offer the
  fixed-weight sandwich as a fast (slightly anticonservative) option with a classed
  warning, mirroring survatr's existing unbalanced-panel warning pattern.
- **Reject** combinations survatr cannot weight (e.g. designs without a usable
  person-period mapping) with a classed error.

## API design

```r
# Case-cohort -> marginal absolute risk + risk difference
fit <- matcha(cohort_pp, outcome = "event", exposure = "x",
              design = case_cohort(subcohort = "sub", time = "t"),
              confounders = ~ age, estimator = "surv_gcomp")   # -> survatr weighted
contrast(fit, type = "difference", times = c(2, 5))            # marginal RD(t)
rmst_difference(fit, horizon = 5)

# NCC with Samuelsen KM weights -> marginal survival contrast
matcha(ncc_pp, outcome = "event", exposure = "dose",
       design = nested_cc(strata = "set", time = "t", weights = ncc_weights("km")),
       estimator = "surv_ipw")
```

## Support matrix

| Design | Weight | Estimator (survatr) | Estimand | Variance | Status |
|---|---|---|---|---|---|
| case-cohort | Borgan IPW | surv_gcomp | F_x(t), RD, RR, RMST | bootstrap | needs-test |
| NCC | Samuelsen KM | surv_ipw | RD, RR | bootstrap | needs-test |
| case-cohort | Borgan | surv_aipw | RD (DR) | bootstrap | needs-test |
| any | fixed-weight sandwich | — | RD | sandwich + warning | smoke |
| design w/o PP mapping | — | — | — | ⛔ `matchatr_no_person_period` |

## Implementation plan

- `R/causal_survival_sampled.R` — `fit_surv_sampled()`: convert design → person-period
  (reuse `survatr::to_person_period`-style prep), attach inclusion weights, call
  `survatr::surv_fit()` with the chosen estimator + weights, wrap the result; route
  `contrast()`/`rmst_difference()` to survatr.
- `R/variance_samuelsen.R` (shared with Phase 7) — sampling-variance correction;
  design-aware bootstrap refitter.

## Variance / inference notes

Marginal contrasts under sampling: design-aware bootstrap is the reference (resample the
full cohort, redraw the NCC/case-cohort sample, recompute weights, refit via survatr).
The fixed-weight survatr sandwich ignores weight-estimation variance → anticonservative;
gate it behind a classed warning. For known/fixed inclusion probabilities the sandwich is
appropriate.

## Oracle testing strategy

- **survatr on the full cohort** is the truth: simulate a cohort with known marginal
  RD(t)/RMST, run survatr on the full cohort for truth, then draw an NCC/case-cohort
  sample and confirm the design-weighted survatr estimate recovers it.
- Cross-check classical HR (Phase 6/7 via `cch`/`multipleNCC`) is consistent with the
  marginal survival curves where the model is correct.
- Bootstrap vs fixed-weight sandwich: confirm the sandwich is anticonservative and the
  bootstrap attains nominal coverage (Monte-Carlo).

## Chunk plan

1. Person-period conversion + `fit_surv_sampled()` (surv_gcomp) + full-cohort truth
   oracle + no-PP rejection.
2. NCC Samuelsen-weighted surv_ipw + design-aware bootstrap.
3. surv_aipw (DR) + fixed-weight sandwich option with warning.

## Deferred items

Weight calibration for efficiency (Phase 12), competing risks under sampling (survatr
supports cause-specific hazards — a later extension), transportability.
