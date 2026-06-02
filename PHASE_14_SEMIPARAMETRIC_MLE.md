# Phase 14 — Semiparametric Maximum Likelihood (NPMLE) for Case-Cohort / NCC

> **Status: DESIGN**
> Book chapters: 21 (MLE for Case-Cohort and Nested CC).

## Scope

**In:** the semiparametric-efficient full-likelihood estimator (NPMLE) that uses the full
cohort information — marginalizing over the covariate distribution via an EM algorithm —
for proportional-hazards and semiparametric transformation models (Zeng & Lin 2006/2014).
Discrete Phase-1 (finite atoms) and continuous Phase-1 (kernel-smoothed) cases.

**Out:** the simpler IPW / pseudo-likelihood (Phases 6–7) and MI (Phase 13) alternatives,
which are usually preferred in practice; transformation-model zoo beyond PH/PO.

## Key design decisions

- **NPMLE attains the semiparametric efficiency bound** — the most efficient estimator,
  but the heaviest to compute and the least available in CRAN. Position it as the
  "efficiency ceiling" reference, not the default.
- **EM algorithm**: E-step computes conditional expectations of the missing X for
  unsampled units; M-step maximizes the complete-data semiparametric likelihood for
  (β, Λ₀, covariate-distribution atoms). Discrete Z → finite atoms p_{sl}; continuous Z →
  kernel local likelihood with bandwidth a_N = N^{-1/(3+d_z)}.
- **GENUINELY NEW CODE** — there is no off-the-shelf CRAN package (only Zeng & Lin's own
  code; `addhazard` covers the additive-hazards special case). This is the most
  implementation-heavy phase; gate it behind clear "advanced/optional" framing and make
  it the lowest-priority core phase.
- **Variance** via the observed-data information matrix (Louis 1982) for discrete Z; a
  profile-likelihood numerical estimate (with perturbation constant h_N) for continuous Z.

## API design

```r
fit <- matcha(cohort, outcome = "case", exposure = "x",
              design = case_cohort(subcohort = "sub", time = "t"),
              confounders = ~ z, estimator = "npmle",
              transform = "ph")     # or "po" / Box-Cox
contrast(fit)    # efficient HR + CI
```

## Support matrix

| Phase-1 Z | Model | Method | Variance | Status |
|---|---|---|---|---|
| discrete | PH | EM-NPMLE | observed info (Louis) | needs-test |
| discrete | transformation (PO/Box-Cox) | EM | observed info | smoke |
| continuous | PH | kernel EM | profile likelihood | smoke |

## Implementation plan

- `R/npmle.R` — the EM driver: E-step (conditional X expectations / atom weights),
  M-step (semiparametric likelihood maximization for β, Λ₀, atoms), convergence control;
  transformation-model link.
- `R/variance_npmle.R` — Louis observed-information (discrete) and profile-likelihood
  (continuous) variance.

## Variance / inference notes

Asymptotic normality with covariance attaining the efficiency bound. Discrete Z: Louis
formula on the observed-data information. Continuous Z: numerical profile likelihood,
bandwidth-sensitive. Document the computational cost and the bandwidth/perturbation
choices as the main practical caveats.

## Oracle testing strategy

- No CRAN oracle for the general NPMLE. Validate by: (a) limiting-case agreement — when
  the subcohort = full cohort, NPMLE = full-cohort `coxph`; (b) efficiency comparison —
  NPMLE SE ≤ IPW (Phase 7) SE on the same DGP; (c) `addhazard` for the additive-hazards
  special case; (d) Monte-Carlo coverage of the profile-likelihood CI. Flag as
  "smoke/limiting-case only — no general oracle" in the coverage matrix.

## Chunk plan

1. Discrete-Z EM-NPMLE for PH + Louis variance + full-cohort limiting-case oracle.
2. Transformation models (PO/Box-Cox).
3. Continuous-Z kernel EM + profile-likelihood variance (smoke).

## Deferred items

High-dimensional Z reduction; internal time-dependent covariates; broader transformation
families. This is the lowest-priority core phase — IPW/MI (Phases 7/13) cover most needs.
