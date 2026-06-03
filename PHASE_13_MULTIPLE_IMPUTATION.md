# Phase 13 — Multiple Imputation for the Missing-by-Design Covariate

> **Status: DESIGN**
> Book chapters: 20 (Multiple Imputation for Sampled Cohort Data).

## Scope

**In:** multiple imputation of the expensive covariate X (missing by design for unsampled
cohort members) followed by full-cohort Cox regression and Rubin pooling — the
White-Royston approximate method and the Bartlett rejection-sampling (`smcfcs`) method;
auxiliary variables; competing risks / left-truncation / stratified Cox extensions.

**In (added 2026-06-02):** design-aware MI for *sporadic* missing covariates in
matched/unmatched CC — congenial with **conditional** logistic regression for matched CC
(Seaman & Keogh 2015), reusing the `mice`/`smcfcs` engine but with matchatr-owned setup
(see the literature review below). Unmatched CC can lean on `causatr::causat_mice`; matched
CC cannot (it is not clogit-congenial).

**Out:** calibration-via-MI (Phase 12).

## Key design decisions

- **MAR by design**: X is missing only for unsampled units; the imputation model is the
  conditional X | (Z, δ, Λ̂₀(t)) implied by the Cox likelihood. Two methods:
  - **Approximate (White & Royston 2009)** — regress X on (Z, δ, Nelson-Aalen Λ̂₀(t)),
    draw K parameter vectors, impute. Fast; minor downward bias for large effects. Via
    `mice`.
  - **Rejection sampling (Bartlett et al. 2015)** — sample X from a proposal, accept via
    the Cox-likelihood rule; unbiased, handles nonlinear terms. Via `smcfcs`.
- **Include auxiliary Phase-1 predictors of X** (e.g. an FFQ for diet) for large
  efficiency gains.
- **Delegate imputation to `mice` / `smcfcs`** and pooling to Rubin's rules; matchatr
  supplies the Nelson-Aalen term and the design-aware imputation setup, then fits the
  full-cohort Cox per imputation. **Reuse causatr's pooling** (`pool_rubin`) where it
  applies, for consistency across the etverse.
- **Do not impute the outcome / event** — use Y, δ as predictors in the imputation model
  (mirrors causatr's MI invariant).

## API design

```r
fit <- matcha(cohort, outcome = "case", exposure = "x",   # x measured on sample only
              design = case_cohort(subcohort = "sub", time = "t"),
              confounders = ~ z, estimator = "mi_cox",
              impute = mi_spec(method = "smcfcs", m = 20, aux = ~ ffq))
contrast(fit)    # pooled HR (Rubin)
```

## Support matrix

| Method | Design | Engine | Pooling | Status |
|---|---|---|---|---|
| approximate (W&R) | NCC / case-cohort | mice + survival::coxph | Rubin | needs-test |
| rejection (Bartlett) | NCC / case-cohort | smcfcs | Rubin | needs-test |
| + auxiliary | both | mice/smcfcs | Rubin | needs-test |
| competing risks / truncation | both | smcfcs | Rubin | smoke |
| impute the outcome | — | — | ⛔ `matchatr_impute_outcome` |

## Implementation plan

- `R/mi_cox.R` — `mi_spec()`; build the Nelson-Aalen term; run `mice`/`smcfcs`; loop
  full-cohort `coxph` over completed datasets; pool via Rubin (reuse `causatr::pool_rubin`
  if exported, else local). Reject outcome imputation.

## Variance / inference notes

Rubin's rules: total variance W̄ + (1+1/K)B̄. The approximate method has a small downward
bias for large effects (document); rejection sampling is unbiased. CIs use the
Barnard-Rubin df.

## Oracle testing strategy

- `mice` / `smcfcs` engines + oracle. Truth-based: simulate a cohort with known Cox β,
  set X missing for unsampled units, impute + pool, confirm β recovery and that an
  informative auxiliary improves efficiency (handbook §20.6.2).
- Compare to the IPW (Phase 7) and full-cohort estimates on the same DGP — MI should be
  at least as efficient when the imputation model is good.

## Chunk plan

1. Approximate MI (`mice` + Nelson-Aalen) + Rubin pooling + truth oracle + outcome-impute
   rejection.
2. Rejection sampling via `smcfcs` + auxiliary variables.
3. Competing-risks / left-truncation / stratified extensions (smoke).

## Deferred items

MI-for-calibration (Phase 12), full semiparametric MLE alternative (Phase 14).

> **Revision (2026-06-02):** the earlier plan to defer *all* non-survival CC
> confounder MI to `causatr::causat_mice` is **methodologically insufficient for
> matched CC** — see the literature review below. `causat_mice` is a generic
> wrapper and does not enforce congeniality with **conditional** logistic
> regression or the matched-set imputation structure. matchatr must own the
> design-aware MI setup for matched/unmatched CC; only the generic engine and
> Rubin pooling are reused.

## Literature review: missing-data methods across these designs (2026-06-02)

Two **distinct** missing-data problems arise here; they need different machinery.

**(1) Missing *by design* — the expensive covariate X unmeasured for unsampled
cohort members (NCC, case-cohort).** This is the section above. The covariate is
MAR given the sampling, and the full cohort supplies the imputation information.
- Keogh & White (2013, *Statistics in Medicine*) — impute X for unsampled members
  using full-cohort data; the substudy-only and intermediate variants.
- Keogh, Seaman, Bartlett & Wood (2018, *Biometrics* 74(4):1438) — unifies
  missing-by-design **and** missing-by-chance for NCC/case-cohort; adapts the
  White-Royston approximate method (MI-approx) and the Bartlett SMC method
  (MI-SMC), plus the Seaman-Keogh "MI matched set" approach for NCC without
  full-cohort info. The intermediate approach is more robust to imputation-model
  misspecification than the full-cohort approach.
- Marti & Chavance (2011, *Statistics in Medicine*) — MI for case-cohort HRs.

**(2) Missing *by chance* — sporadic missing confounders/exposure in (matched or
unmatched) case-control.** This is the gap the current doc under-specified.
- **The analysis model dictates the imputation model (congeniality;** Meng 1994,
  *Statistical Science* 9(4):538–558; Bartlett et al. 2015). For **matched** CC the
  substantive model is **conditional** logistic regression, so the imputation
  model must be compatible with *it*, not with a marginal logistic fit.
- **Seaman & Keogh (2015, *Biometrics* 71(4):1150–1159)** is the key reference for
  matched CC. Two approaches: (a) **per-individual** — impute conditioning on the
  matching variables *and* disease status, then restore the matching for the
  conditional-logistic analysis; (b) **per-matched-set** — impute the whole set
  jointly with a set-level random effect, avoiding a parametric model for the
  matching variables (useful when they are unavailable/awkward). FCS with a
  **restricted general-location model is compatible** with conditional logistic;
  naive **normal / latent-normal joint imputation is NOT** compatible. FCS is
  available in standard tooling (`mice` in R) — no bespoke Bayesian software
  needed.
- **Always include the outcome (case status) as a predictor in the imputation
  model** (Sterne et al. 2009, *BMJ* 338:b2393; Moons et al. 2006). Omitting it
  biases the imputed associations toward the null. This is *why* MI is valid for
  case-control data at all: the logistic odds ratio is symmetric in the roles of
  exposure and outcome, so imputing the missing covariate **given** case status is
  congenial with the prospective OR analysis.
- **Substantive-Model-Compatible FCS** (Bartlett, Seaman, White & Carpenter 2015,
  *Statistical Methods in Medical Research* 24(4):462–487; R package `smcfcs`) is
  the method of choice when the substantive model has **nonlinear terms,
  interactions, or is a Cox model**. `smcfcs` supports linear/logistic/Poisson/Cox
  substantive models out of the box; **conditional logistic is not a built-in
  smcfcs family**, so matched CC needs either the Seaman-Keogh FCS construction or
  the per-matched-set route — confirming this cannot be a plain `causat_mice` call.
- **Survival substantive model:** include the event indicator δ and the
  **Nelson-Aalen cumulative hazard** Λ̂(t) (not raw time) as imputation predictors
  (White & Royston 2009, *Statistics in Medicine* 28:1982–1998).

**Pitfalls / what NOT to do.** The **missing-indicator method is biased** for
case-control associations and is not an MI substitute. Complete-case analysis is
unbiased only when missingness is independent of the outcome given the covariates
in the model (often violated) and is always inefficient. Compatibility is
necessary but **not sufficient** — a compatible-but-misspecified imputation model
still biases estimates (recent simulation work), so the imputation model needs the
right functional form and informative auxiliaries, not just structural compatibility.

**Implications for matchatr.** (i) Keep the missing-by-design survival MI above for
NCC/case-cohort (Phase 13 core). (ii) Add design-aware MI for matched/unmatched CC:
unmatched CC can lean on `causat_mice` / `smcfcs` with case status + a logistic
substantive model, but **matched CC requires a conditional-logistic-congenial
setup** (Seaman-Keogh per-individual or per-matched-set), which matchatr must
construct itself — generic `causat_mice` will not. (iii) In every case the
imputation model must include the outcome, and matchatr should reject imputing the
outcome/event itself (`matchatr_impute_outcome`). Engine (`mice`/`smcfcs`) and
Rubin pooling are reused; the design-aware *setup* is matchatr's.
