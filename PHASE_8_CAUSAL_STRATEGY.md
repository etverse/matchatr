# Phase 8 — Causal Estimation Strategy for Case-Control Designs

> **Status: DELIVERED (strategy / decision doc — ships no estimator code).**
> The strategy, the q0 / weight contract, and the build-vs-reuse map below are in
> place: `unmatched_cc(prevalence = q0)` + `check_prevalence()`, the `ccw_*`
> dispatch family valid on any design, the `matchatr_missing_prevalence` guard,
> and `causatr` / `survatr` as Imports. Track 2 (case-control weighting) is
> implemented in Phase 9; Track 3 (design-weighted survival) in Phase 10.
> Book chapters: 6 (Causal Inference framing), 12 (two-phase view); methods literature:
> Rose & van der Laan (2008, 2009, 2011, 2014) case-control-weighted estimation.

## Purpose

This is the doc the package owner flagged as the open question: *which causal estimation
strategy is best for case-control-type designs?* It lays out the tracks, states what
each identifies, what it costs, and what it reuses vs needs new code, and recommends a
default. It does not ship code — it sets the contract Phases 9–10 implement.

## The core problem

Classical case-control analysis (Phases 2–7) targets a **conditional** odds ratio or
hazard ratio. That is a valid associational/conditional causal parameter, but it is:
non-collapsible (changes with covariate adjustment even absent confounding), not a
marginal effect, and not on the risk-difference / risk-ratio scale that population
decisions need. The case-control sample is also a **biased sample** of the source
population (cases oversampled). Marginal causal estimands (ATE as RD, RR, marginal OR;
absolute risk) require correcting that sampling bias.

## The four tracks

### Track 1 — Conditional model (handbook-native) — Phases 2–7
Conditional logistic / weighted Cox → conditional OR / HR. **Identifies** the
within-stratum effect. **Cost** none beyond classical fitting. **Reuse:** `survival`.
Always available; the baseline. Not a marginal effect.

### Track 2 — Case-control weighting (Rose & van der Laan) — Phase 9 *(recommended headline)*
Reweight the CC sample back to the source population using the **marginal prevalence
q₀ = P(Y=1)**. With case weights ∝ q₀ and control weights ∝ (1−q₀), the weighted
empirical distribution mimics the cohort, so any cohort estimator applied to the
weighted sample targets the **marginal** estimand. The full family:
- **CCW-g-formula** — fit an outcome model on the weighted sample, standardize → marginal
  RD/RR/OR. *Reuses `causatr` g-comp* with the q₀ weights.
- **CCW-IPW** — fit a propensity model on the weighted sample, Horvitz-Thompson/Hájek →
  marginal effect. *Reuses `causatr` IPW.*
- **CCW-AIPW** — doubly-robust augmented IPW on the weighted sample. *Reuses `causatr`
  AIPW.* This is the etverse's doubly-robust estimator.
- **CCW-TMLE** — targeted maximum likelihood: initial outcome fit + a targeting
  (fluctuation) step solving the efficient influence equation on the weighted sample.
  **NEW CODE** — causatr has no targeted-learning machinery, so the fluctuation step is
  matchatr's own (its own sub-phase in Phase 9).

**Identifies** marginal RD/RR/OR (and works for independent, matched, and nested CC).
**Cost** requires q₀ (known, or estimated from the full cohort). **Variance:** influence
function with an extra term when q₀ is estimated; or bootstrap.

### Track 3 — Design / inclusion weighting (Samuelsen / Borgan) — Phases 6, 7, 10
Inverse probability of *selection* weights (subcohort / risk-set sampling). Reweights a
sampled cohort to the full cohort. **Identifies** cohort-level regression parameters and
absolute risk under the sampling design. **Distinct from Track 2** (selection bias vs
case-oversampling bias; different variance). **Reuse:** `multipleNCC`, `survival::cch`,
and `survatr` for person-period causal survival (Phase 10).

### Track 4 — Two-phase / survey + full MLE / MI (efficiency) — Phases 11–14
Treat the design as a two-phase sample; use survey calibration (Track-2/3 weights tuned
to Phase-1 totals) for efficiency, or semiparametric MLE / multiple imputation for the
missing-by-design covariate. Layers on top of Tracks 2–3.

## Recommendation

- **Default causal path = Track 2 (case-control weighting).** It is the cleanest mapping
  onto the etverse: the q₀ weights are *observation weights*, so CCW-g-formula / IPW /
  AIPW are immediate reuses of `causatr`, and the package gains marginal RD/RR/OR for
  every CC design with one weighting concept.
- **Ship CCW-g-formula and CCW-AIPW first** (g-formula = simplest correct marginal
  estimator; AIPW = doubly-robust, already in causatr). Add **CCW-TMLE** as a dedicated
  sub-phase since it is the only genuinely new engine.
- **Track 3 owns the time-to-event sampled designs** (NCC/case-cohort absolute risk and
  survival contrasts) via `survatr` (Phase 10).
- **Track 4 is efficiency**, added once Tracks 2–3 are solid.

## Rejected / deferred alternatives

- *Reading marginal effects off a conditional logistic fit* — rejected: non-collapsible,
  and no q₀ means no baseline risk. Marginal effects require Track 2.
- *Building a general targeted-learning engine in matchatr* — deferred: only the CCW
  fluctuation step is needed now; a full TL framework is out of scope (and arguably an
  etverse-wide concern, not matchatr's).
- *Conflating q₀ weights with Samuelsen design weights* — rejected by design (distinct
  objects, distinct variance).

## The q₀ (prevalence) interface (contract for Phase 9)

- `unmatched_cc(prevalence = q0)` / `nested_cc(...)` carry q₀.
- q₀ may be **known** (literature/registry) → treated as fixed in the IF; or **estimated
  from the full cohort** → adds a variance term. The fit records which.
- Missing q₀ for any `ccw_*` estimator → classed error `matchatr_missing_prevalence`.

## What this phase delivers

A written strategy + the q₀/weight contracts + the build-vs-reuse map. No estimator code.
Phases 9 and 10 implement Tracks 2 and 3 respectively.
