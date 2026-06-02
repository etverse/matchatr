# Phase 19 — Self-Controlled Case Series (SCCS)

> **Status: DESIGN (extension / lower priority)**
> Book chapters: 22 (The Self-Controlled Case Series Method).

## Scope

**In:** the SCCS design — a self-matched method using only cases (r ≥ 1 events),
conditioning on the observation period and event count so that fixed multiplicative
confounders cancel; estimation of the log relative incidence for transient exposures
(e.g. vaccine safety). Parametric (age groups + risk periods via conditional Poisson /
clogit) and semiparametric (cumulative age function) baselines.

**Out:** the case-control / cohort-sampling designs (Phases 1–18). SCCS is a distinct
module that shares only the conditional-likelihood machinery.

## Key design decisions

- **SCCS likelihood** = ∏ λ(t_ij | exposure) / [∫ λ(s) ds]^{r_i}, fit as a **conditional
  Poisson** (Poisson GLM with offsets) or via `survival::clogit` — only cases are needed,
  no controls. Fixed confounders eliminated by the self-control conditioning.
- **Proportional incidence** λ_i(t) = φ_i θ(t) exp(exposure·β): parametric age via step
  functions / risk periods; semiparametric via a cumulative age function jumping at event
  times (high compute cost). Spline / fractional-polynomial age as middle-ground.
- **Key assumptions to enforce/validate**: events must not alter the observation period
  (event-dependent observation) or future exposure (event-dependent exposure); for the
  basic model events are rare and non-recurrent (or recurrences independent). Provide
  classed warnings/checks for these.
- **Largely reuse** the conditional Poisson / clogit + `mgcv` spline machinery; only the
  semiparametric jump-at-events baseline is genuinely new (and optional).
- **Separate design constructor** `sccs(id, event_times, obs_period, exposure_history)`,
  distinct from the case-control designs.

## API design

```r
fit <- matcha(cases, outcome = "event_times", exposure = "vaccine_window",
              design = sccs(id = "id", obs = c("start", "end"), age_groups = breaks),
              estimator = "sccs")
contrast(fit)    # relative incidence (RR) + CI
```

## Support matrix

| Baseline | Engine | Variance | Status |
|---|---|---|---|
| parametric age + risk periods | conditional Poisson / clogit | Poisson info | needs-test |
| spline / fractional-poly age | mgcv | info | smoke |
| semiparametric (jumps at events) | custom | Fisher info | smoke |
| event-dependent obs/exposure | — | (warn / classed check) | needs-test |

## Implementation plan

- `R/sccs.R` — `sccs()` design constructor; expand cases to risk-period × age-group
  intervals; fit conditional Poisson (`glm`) / `clogit`; relative-incidence contrast;
  assumption checks. Optional semiparametric baseline.

## Variance / inference notes

Poisson standard errors (parametric); Fisher information (semiparametric). Asymptotic
relative efficiency to a cohort is near 1 for common exposures / long risk periods, lower
for rare exposures / short periods (document). Sample-size formula from Ch22 §22.4.

## Oracle testing strategy

- `SCCS` CRAN package (if available) and conditional Poisson via `glm`/`clogit` as
  oracles. Truth-based: simulate a case series with a known relative incidence and risk
  window; confirm recovery and that fixed confounders genuinely cancel.

## Chunk plan

1. `sccs()` design + interval expansion + conditional-Poisson fit + oracle + assumption
   checks.
2. Spline / fractional-polynomial age baselines.
3. Semiparametric jump-at-events baseline (smoke).

## Deferred items

Event-dependent exposure/observation advanced handling; recurrent-event extensions;
calendar-time effects. This is an optional module — not on the core case-control path.
