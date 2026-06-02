# Phase 5 — Nested Case-Control: Risk-Set Sampling and Conditional Partial Likelihood

> **Status: DESIGN**
> Book chapters: 16 (Cohort Sampling Overview), 18 (NCC Counting-Process Approach).

## Scope

**In:** the classical nested case-control analysis within a cohort — risk-set
(incidence-density) sampling of m controls per case, and estimation of Cox hazard
ratios via the conditional partial likelihood (each sampled risk set as a stratum).
Simple random control sampling and counter-matching (stratified risk-set sampling);
additional matching on population strata. A control-sampling helper to *generate* an
NCC dataset from a cohort (for simulation/teaching).

**Out:** breaking the matching / reusing controls via IPW (Phase 7), absolute-risk and
non-PH models (Phase 7), marginal causal survival under sampling (Phase 10).

## Key design decisions

- **NCC partial likelihood = `clogit`/`coxph` with the sampled risk set as a stratum.**
  L_ncc(β) = ∏_j exp(β·x_{i_j}) / Σ_{k ∈ R̃(t_j)} exp(β·x_k), identical in form to the
  matched-CC conditional likelihood (Phase 3). Fit via `survival::clogit` on
  case-control sets, or `survival::coxph` with strata. Same engine, different design
  semantics. (Invariant in `hard-rules.md`.)
- **OR = HR exactly** under risk-set matching — no rare-disease assumption (Prentice &
  Breslow 1978). Contrasts report HR directly; do not attach a rare-disease caveat.
- **Counter-matching uses log-sampling-weights as offsets.** L_cm(β) carries weights
  w_k(t) = n_{s(k)}(t)/m_{s(k)}; implemented by entering log-weights as a Cox offset
  (Langholz & Borgan 1995). The `counter_matched()` design carries the stratum sizes.
- **Time scale must be fixed before control selection** for the partial-likelihood
  analysis (unlike case-cohort). Document; the alternative time scale needs the IPW
  reformulation (Phase 7).
- **Reject** non-time-to-event use of `nested_cc` and risk sets with no eligible
  controls (`matchatr_empty_risk_set`).

## API design

```r
# Analyse an existing NCC dataset (case + m matched controls per set)
fit <- matcha(ncc, outcome = "case", exposure = "dose",
              design = nested_cc(strata = "setid", time = "agexit"),
              confounders = ~ ageRx, estimator = "clogit")
contrast(fit)                       # HR + CI

# Generate an NCC sample from a cohort (helper)
ncc <- sample_ncc(cohort, time = "t", event = "d", m = 2,
                  match = ~ entry_stratum)   # wraps Epi::ccwc
```

## Support matrix

| Sampling | Estimator | Estimand | Variance | Status |
|---|---|---|---|---|
| simple m:1 risk-set | clogit | HR | partial-lik info | needs-test |
| m:1 + extra matching (pop. strata) | clogit (stratified) | HR | partial-lik info | needs-test |
| counter-matched | coxph + offset | HR | weighted partial-lik | needs-test |
| nested_cc, non-survival outcome | — | — | — | ⛔ `matchatr_bad_outcome` |

## Implementation plan

- `R/clogit.R` (shared with Phase 3) — NCC formula builder routing the sampled risk set
  to `strata()`.
- `R/risk_set_sampling.R` — `sample_ncc()` (wrap `Epi::ccwc`; fall back to a base
  implementation), counter-matching sampler, eligibility/empty-risk-set checks.
- `R/weighted_cox.R` — counter-matching offset handling (the weighted partial
  likelihood; full weight machinery in Phase 7).

## Variance / inference notes

Conditional partial-likelihood information matrix (standard Cox/clogit output) — valid
because controls are NOT reused across sets in the classical NCC analysis (each set is
an independent stratum). Counter-matching: weighted partial-likelihood information with
the offset. Robust/cluster variance only becomes necessary when controls are reused
(Phase 7).

## Oracle testing strategy

- `survival::clogit` / `survival::coxph` on the same NCC data (we wrap them; assert
  exact agreement).
- Truth-based: simulate a cohort with known Cox β, draw an NCC sample with `sample_ncc`,
  confirm the partial-likelihood estimate recovers β with efficiency ≈ m/(m+1) at β=0.
- Counter-matching: cross-check against `Epi` / `multipleNCC` examples (radiation/breast-
  cancer, handbook §19.3.1).
- Compare full-cohort `coxph` β to NCC β on the same simulated cohort (should agree
  within sampling error).

## Chunk plan

1. NCC `clogit` analysis path + survival-outcome rejection + oracle.
2. `sample_ncc()` control-sampling helper (+ `Epi::ccwc` wrap) + efficiency test.
3. Counter-matching offset path + oracle.

## Deferred items

Control reuse / IPW (Phase 7), absolute risk and non-PH models (Phase 7), marginal
causal contrasts under NCC sampling (Phase 10).
