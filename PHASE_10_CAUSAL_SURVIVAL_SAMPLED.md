# Phase 10 — Design-Weighted Causal Survival for Sampled Cohorts

> **Status: DESIGN** (revised after the 2026-06-11 survatr-reuse audit — see the
> **survatr prerequisites** section: survatr's `surv_ipw` rejects external weights
> and there is no `surv_aipw` yet, so the support matrix and chunk plan changed).
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

- **Delegate causal survival to `survatr`.** survatr does pooled-logistic discrete-time
  hazard g-computation (`surv_gcomp`), propensity IPW (`surv_ipw`), and time-varying ICE
  (`surv_ice`) on person-period data, with sandwich/bootstrap variance and
  survival / risk / risk-difference / risk-ratio / RMST / RMTL / quantile / CIF contrasts.
  matchatr's job is to (a) convert the sampled design to person-period form and (b) supply
  the inclusion-probability weights as survatr observation weights.
- **Only `surv_gcomp` accepts external observation weights** (verified: `gcomp_survival.R`
  broadcasts `weights` to the at-risk rows and switches to quasibinomial). `surv_ipw`,
  `surv_ice`, and competing-risks **reject** external weights (their weights are the
  fitted propensity / IPCW weights, a different object from design weights). So the
  design-weighted path for **every** sampled design routes through `surv_gcomp + design
  weights`, NOT `surv_ipw` — the inclusion weights are observation weights into the
  standardized hazard model, not a substitute for a propensity model. (Composing design
  weights with survatr's stabilized IPW would need a survatr change — see prerequisites.)
- **survatr requires rectangular person-period data** (it rejects a ragged panel). A
  sampled NCC / case-cohort is not rectangular, so `fit_surv_sampled()` must pad the
  retained subjects to a common time grid before calling `surv_fit()`.
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

# NCC with Samuelsen KM weights -> marginal survival contrast.
# The Samuelsen inclusion weights are observation weights into surv_gcomp (the
# only survatr estimator that accepts external weights), NOT a propensity model,
# so this is surv_gcomp + design weights, not surv_ipw.
matcha(ncc_pp, outcome = "event", exposure = "dose",
       design = nested_cc(strata = "set", time = "t", weights = ncc_weights("km")),
       estimator = "surv_gcomp")
```

## Support matrix

| Design | Weight | Estimator (survatr) | Estimand | Variance | Status |
|---|---|---|---|---|---|
| case-cohort | Borgan IPW | surv_gcomp | F_x(t), RD, RR, RMST | bootstrap | needs-test |
| NCC | Samuelsen KM | surv_gcomp (design weights) | F_x(t), RD, RR, RMST | bootstrap | needs-test |
| case-cohort / NCC | Borgan / Samuelsen | surv_aipw (DR) | RD (DR) | bootstrap | **blocked on survatr `surv_aipw`** |
| any | fixed-weight sandwich | surv_gcomp | RD | sandwich + warning | smoke |
| design w/o PP mapping | — | — | — | ⛔ `matchatr_no_person_period` |

Both sampled designs now route through **`surv_gcomp` + design weights** (the only
survatr estimator that accepts external weights); the earlier `surv_ipw` row for NCC was
removed because survatr's `surv_ipw` rejects external weights. The DR (`surv_aipw`) row
is blocked on a survatr addition (see prerequisites).

## survatr prerequisites (from the 2026-06-11 reuse audit)

Two survatr capabilities this phase assumes do not exist yet. Resolve before the
corresponding chunk:

1. **`surv_gcomp` external weights — present (no change needed).** Confirmed survatr's
   `surv_fit(estimator = "gcomp", weights = )` accepts external observation weights, so the
   primary design-weighted path (Chunks 1–2) needs no survatr change.
2. **`surv_aipw` (doubly-robust survival) — MISSING; must be added to survatr.** survatr
   ships only `gcomp` / `ipw` / `ice`; there is no AIPW/DR survival estimator, so the DR
   row (Chunk 3) is blocked until survatr gains one. A DR treatment-specific survival /
   RMST estimator is **theoretically sound and well-established**: the locally efficient
   augmented-IPW (AIPWCC) estimator is consistent for the treatment-specific survival
   distribution if **either** the outcome (survival/hazard) model **or** the propensity +
   censoring models are correct. Primary references (verified):
   - Robins & Rotnitzky (1992) — the AIPW / locally efficient estimating-equation theory.
   - Hubbard, van der Laan & Robins (2000), *Statistical Models in Epidemiology, the
     Environment, and Clinical Trials* — locally efficient survival estimation.
   - **Zhang & Schaubel (2012)**, "Contrasting treatment-specific survival using
     double-robust estimators", *Statistics in Medicine* 31(30): 4255–4268 — DR for
     treatment-specific survival **and RMST**, consistent if either a logistic treatment
     model or a Cox death-hazard model is correct. The closest fit to survatr's
     pooled-logistic-hazard + propensity framework, and the recommended template.
   - **Bai, Tsiatis & O'Brien (2013)**, "Doubly-robust Estimators of Treatment-specific
     Survival Distributions in Observational Studies with Stratified Sampling",
     *Biometrics* 69(4) — locally efficient AIPWCC for treatment-specific survival
     **under stratified sampling**, directly relevant to the NCC / case-cohort sampled
     designs this phase targets.

   **Action:** file a survatr feature request for `estimator = "surv_aipw"` (a DR survival
   estimator accepting external design weights), templated on Zhang & Schaubel (2012);
   matchatr's Chunk 3 then delegates to it exactly as Chunks 1–2 delegate to `surv_gcomp`.
   Until then, Chunk 3 is deferred (the singly-robust `surv_gcomp` path still ships in
   Chunks 1–2).
3. **`surv_ipw` external weights — by design rejected; not pursued.** survatr's `surv_ipw`
   is a propensity-weighted MSM, so it rejects external design weights. Composing design
   weights with stabilized IPW is a possible future survatr extension but is NOT needed:
   the design-weighted path uses `surv_gcomp` (decision above), so no survatr change is
   requested here.

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

1. Person-period conversion (+ rectangular padding) + `fit_surv_sampled()` (surv_gcomp,
   case-cohort + Borgan weights) + full-cohort truth oracle + no-PP rejection.
2. NCC Samuelsen-weighted **surv_gcomp** (design weights, not surv_ipw) + design-aware
   bootstrap + fixed-weight sandwich option with a classed warning.
3. surv_aipw (DR) — **blocked on a survatr `surv_aipw` estimator** (see prerequisites);
   deferred until survatr adds it (templated on Zhang & Schaubel 2012). The singly-robust
   surv_gcomp path from Chunks 1–2 ships regardless.

## Deferred items

Weight calibration for efficiency (Phase 12), competing risks under sampling (survatr
supports cause-specific hazards — a later extension), transportability.
