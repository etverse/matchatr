# matchatr

Causal inference for **(matched) case-control, nested case-control (NCC), and
case-cohort** study designs. Provides design-faithful classical estimators
(conditional logistic regression, risk-set-sampled / weighted Cox, Prentice /
Self-Prentice / Borgan case-cohort) and marginal causal effects via **case-control
weighting** (the Rose & van der Laan g-formula / IPW / AIPW / TMLE family) and
**design-based inclusion weighting** (Samuelsen, Borgan). Part of the
[etverse](https://github.com/etverse) ecosystem; delegates estimation to `causatr`
(g-comp / IPW / AIPW + sandwich/bootstrap variance) and `survatr` (causal survival
on person-period data) wherever possible.

> **Status: classical odds-ratio engines landing.** The PHASE_1 foundation
> (design taxonomy, unified `matchatr_design` S3 object + six constructors, the
> `matcha()` fit verb, the `(design, estimator)` dispatch + validation layer) is
> in place. **PHASE_2 (unmatched case-control) is complete**: `matcha()` fits
> the conditional-OR logistic via `stats::glm` (or a pluggable `model_fn` such as
> `mgcv::gam`) for binary / continuous / categorical / ordinal-trend exposures,
> and `estimator = "mh"` computes the Mantel-Haenszel stratified OR with
> Robins-Breslow-Greenland variance. **PHASE_3 (matched case-control) is
> complete**: `matcha(design = matched_cc(...), estimator = "clogit")` fits the
> conditional likelihood via `survival::clogit` and reports the conditional OR
> through the shared `conditional_or_result()` assembly, `estimator =
> "mcnemar"` computes the 1:1 matched-pair OR = n10/n01 with Var(log OR) =
> 1/n10 + 1/n01 in closed form (rejecting M:1 / richer matching toward
> `clogit`), and `effect_modifier = "m"` fits `outcome ~ exposure * m + ... +
> strata(set)` so `contrast(type = "or")` reports the exposure's
> **stratum-specific** conditional OR per modifier level (β_x at the reference,
> β_x + β_{x:level} elsewhere) via `stratum_specific_or_result()`. M:1 and
> variable-ratio matching need no special handling (the conditional likelihood
> treats any matched-set composition uniformly). `contrast(type = "or")` reports
> the OR(s); RD/RR are rejected as unidentified without q0. **PHASE_4 (multiple
> case/control groups) is complete**: `matcha(design = unmatched_cc(),
> estimator = "polytomous", reference = ...)` fits the baseline-category
> multinomial logistic via `nnet::multinom` for a ≥3-group outcome, and
> `contrast(type = "or")` reports each non-reference subtype's exposure odds
> ratio versus the reference (one OR per subtype × exposure-coefficient,
> information-matrix Wald interval) while `tidy()` adds a `y.level` column; a
> two-group / numeric / logical outcome is rejected toward the binary estimators
> (`matchatr_bad_outcome`). The dispatch gained an `outcome_kind` axis so
> `matcha()` resolves the multi-group outcome (`resolve_polytomous_outcome()`)
> instead of the binary one. **PHASE_4 Chunk 2** adds `test_homogeneity(fit)`:
> for each exposure term it runs the canonical **Wald** test of whether the
> exposure odds ratio is constant across the subtypes (H0: β₁ = … = β_M,
> W = (C b)′(C V C′)⁻¹(C b) ~ χ²₍M−1₎) and reports the efficient **GLS-pooled
> common OR** — both computed on the unconstrained `nnet::multinom` fit + its
> information matrix (reusing `multinom_exposure_or()`), so there is no
> constrained refit and continuous confounders are handled directly (chosen over
> a Poisson-surrogate or `VGAM` LRT; matches `riskclustr::eh_test_subtype`).
> Non-polytomous / non-estimated fits are rejected (`matchatr_bad_input` /
> `matchatr_not_estimated`). **PHASE_5 (nested case-control) Chunk 1 is
> complete**: `matcha(design = nested_cc(...), estimator = "clogit")` fits the
> risk-set conditional partial likelihood through the *same* `clogit` engine as
> the matched design (a sampled risk set and a matched set are the same stratum
> construction), and `contrast()` reports the exposure's **hazard ratio** via a
> new contrast scale `type = "hr"` — OR = HR exactly under risk-set
> (incidence-density) sampling, with no rare-disease caveat (Prentice & Breslow
> 1978). The design fixes the scale: `default_contrast_type()` is design-aware
> (`"hr"` for nested, `"or"` for matched) and the shared `conditional_or_result()`
> / `stratum_specific_or_result()` assemblies carry a `type` label (the
> arithmetic is identical). Each conditional design identifies exactly one scale,
> so requesting an OR from a risk-set design (or an HR from a matched design) is
> `matchatr_unidentified_estimand`; the design's `time` column records the
> sampling but the conditional likelihood reads the risk set from `strata`.
> **PHASE_5 Chunk 2 is complete**: the exported `sample_ncc(cohort, time, event,
> m, match, entry)` (`R/risk_set_sampling.R`) generates an analysis-ready NCC
> dataset from a cohort by risk-set (incidence-density) control sampling —
> appending `set`, the per-set `case` indicator, and `risk_time` — with optional
> population-stratum matching and delayed entry. The sampler is native (base
> R/data.table, deterministically seedable via the ambient RNG); `Epi::ccwc` is a
> test oracle, not a runtime dependency (matchatr delegates only to Imports-tier
> *estimation* engines and hand-rolls sampling / closed forms). A case left with
> no eligible control aborts `matchatr_empty_risk_set` (a generation-path failure,
> unlike an uninformative analysis stratum, which `clogit` drops); a late failure
> time with fewer than `m` eligible controls yields a smaller set, not an error.
> The test-only `sample_ncc_riskset()` fixture now delegates to `sample_ncc()`.
> Counter-matching remains a later chunk; PHASE_6+ remain `Status: DESIGN`.

## Guide files

- `FEATURE_COVERAGE_MATRIX.md` — **single source of truth for "what works".** Every
  PR that changes a feature MUST update this file. Records the PHASE_1 design/API
  layer and the PHASE_2 Chunk 1 logistic conditional OR; the remaining estimator
  cells stay pending until their phases land.
- `PHASE_*.md` — per-phase design docs in the project root (see roadmap below). They
  follow the `implement-feature` Step-1b 10-point structure.

## Project structure

This is an R package: `R/` (source), `tests/testthat/` (tests, `test-foo.R` mirrors
`R/foo.R`), `man/` (generated — do not edit), `NAMESPACE` (generated — do not edit),
`vignettes/` (long-form docs). The website is built with `altdoc` + Quarto
(`altdoc/quarto_website.yml`, `lumen` theme) to match the other etverse packages.

### R/ layout (created as phases land)

- **Design + API layer (PHASE_1, implemented):** `cc_design.R` (six design
  constructors + `new_matchatr_design()`), `matcha.R` (the `matcha()` fit verb;
  runs the resolved engine via `run_engine()`), `dispatch.R` (the
  `(design, estimator)` → engine table + `resolve_engine()` + `run_engine()`),
  `contrast.R` (the second-step `contrast()` verb; dispatches per engine),
  `constructors.R` (`new_matchatr_fit()` / `new_matchatr_result()`), the
  validator layer (`checks.R` primitive argument validators, `checks_design.R`
  analysis-role / exposure / effect-modifier / strata-informativeness checks,
  `resolve.R` outcome / exposure / event column coercions — all classed
  `matchatr_*` errors), `print.R`, `tidy.R`, `summary.R`.
  `risk_set_sampling.R` (PHASE_5 Chunk 2 — `sample_ncc()`: native risk-set
  control sampling from a cohort, with population-stratum matching, delayed
  entry, and the `matchatr_empty_risk_set` hard error; counter-matching to
  come). Still to come: `weights_cc.R` (case-control / q₀ weights),
  `weights_design.R` (Samuelsen / Borgan inclusion-probability weights).
- **Classical estimators:** `unconditional.R` (PHASE_2 — `fit_logistic_cc()`
  wraps `stats::glm` / pluggable `model_fn`, plus the conditional-OR contrast and
  the `matchatr_unidentified_estimand` rejection), `mantel_haenszel.R` (PHASE_2
  Chunk 3 — `fit_mh()` closed-form stratified OR + Robins-Breslow-Greenland
  variance), `coef_extract.R` (fitter-agnostic coefficient / variance extraction
  shared across engines — `term_assign()`, `estimable_vcov()`,
  `exposure_coef_index()`, `parametric_positions()`, plus the shared
  `conditional_or_result()` exp(beta)-assembly used by the logistic and clogit
  engines — odds ratio or, for the nested risk-set design, hazard ratio),
  `clogit.R` (PHASE_3 + PHASE_5 Chunk 1 — `fit_clogit()` wraps `survival::clogit`
  for the matched case-control AND nested case-control conditional partial
  likelihood; `contrast_clogit()` reports the conditional OR (matched) or hazard
  ratio (nested, `type = "hr"`), with `reject_offdesign_conditional_scale()`
  enforcing one scale per design), `mcnemar.R` (PHASE_3 Chunk 2 —
  `fit_mcnemar()` closed-form 1:1 matched-pair OR = n10/n01 + McNemar variance,
  mirroring the `mantel_haenszel.R` closed-form precedent),
  `effect_modification.R` (PHASE_3 Chunk 3 — `stratum_specific_or_result()`
  assembles the per-modifier-level conditional OR from the `exposure * modifier`
  clogit fit via a contrast matrix `C V C'`, plus `interaction_coef_index()`),
  `polytomous.R` (PHASE_4 Chunk 1 — `fit_polytomous()` wraps `nnet::multinom`
  for the unmatched ≥3-group multinomial logistic; `contrast_polytomous()` /
  `multinom_exposure_or()` assemble each subtype's exposure OR by term position +
  the `level:predictor` `vcov()` names, and `tidy_multinom()` renders the
  per-equation `y.level` table; the multi-group outcome is resolved by
  `resolve_polytomous_outcome()` in `checks.R`),
  `homogeneity.R` (PHASE_4 Chunk 2 — `test_homogeneity()` runs the per-exposure
  Wald homogeneity test + the GLS-pooled common OR from the stacked subtype
  log-ORs / covariance via `homogeneity_one_term()`, reusing
  `multinom_exposure_or()` and the `C V C'` pattern; `print` / `tidy` methods for
  the `matchatr_homogeneity` class),
  `weighted_cox.R` (NCC IPW Cox), `case_cohort.R` (`cch` wrappers).
- **Causal layer:** `ccw.R` (case-control-weighted dispatch into causatr), `tmle_ccw.R`
  (the NEW targeting step — causatr has no TL), `causal_survival_sampled.R` (design-
  weighted survatr).
- **Inference:** lean on causatr/survatr variance engines; matchatr adds only the
  sampling-variance corrections (`variance_self_prentice.R`, `variance_samuelsen.R`,
  `variance_ccw.R`).
- **S3 + support:** `print.R`, `summary.R`, `tidy.R`, `plot.R`, `coef.R`, `confint.R`,
  `data.R`, `matchatr-package.R`, `zzz.R`.

## Two-step API (PHASE_1; `contrast(type = "or")` computes the unmatched-CC conditional OR as of PHASE_2 Chunk 1, marginal contrasts await the causal phases)

The verb mirrors the siblings (`causatr::causat()`, `survatr::surv_fit()`):

```r
# Unmatched case-control -> conditional OR (implemented)
fit <- matcha(data, outcome = "case", exposure = "x",
              design = unmatched_cc(),
              confounders = ~ age + smoke, estimator = "logistic")
summary(fit)                      # OR table, Wald CIs
contrast(fit, type = "or")        # exposure conditional OR + CI
# contrast(fit, type = "difference")  # -> matchatr_unidentified_estimand (need q0)

# Matched case-control -> conditional OR
fit <- matcha(data, outcome = "case", exposure = "x",
              design = matched_cc(strata = "set"),
              confounders = ~ age + smoke, estimator = "clogit")

# Nested case-control -> risk-set HR (conditional partial likelihood)
fit <- matcha(data, outcome = "case", exposure = "x",
              design = nested_cc(strata = "set", time = "t"), estimator = "clogit")

# Marginal causal effect from a case-control sample (Rose & van der Laan)
fit <- matcha(data, outcome = "case", exposure = "x",
              design = unmatched_cc(prevalence = 0.02),   # q0
              confounders = ~ age + smoke, estimator = "ccw_gformula")

result <- contrast(fit, type = "difference", ci_method = "sandwich")  # marginal RD
```

`matcha()` returns a `matchatr_fit`; `contrast()` (reused) returns a `matchatr_result`.

## Development commands

```r
devtools::load_all()     # load for dev
devtools::test()         # run tests
devtools::check()        # R CMD check
devtools::document()     # regenerate roxygen
```

Shell: `air format .` (format all R files).

## Code style

- `pkg::fun()` for external calls (no bare `library()`)
- `rlang::abort()` / `rlang::warn()` / `rlang::inform()` (not `stop()` / `warning()` /
  `message()`); rejection errors are classed `matchatr_*`
- `data.table` internally; return `data.table` from user-facing functions
- Roxygen on every function, including `@noRd` helpers
- Generous inline comments for math, design rationale, subtle invariants
- Do NOT remove existing comments unless the related code is also removed
- Exported functions at top of files, internal helpers below

## Testing rules

- `expect_snapshot(error = TRUE)` for error conditions
- NEVER delete/mock failing tests — fix the source
- Truth-based simulation tests mandatory for new features (cohort DGP with known truth,
  then sample a CC / NCC / case-cohort from it)
- External oracle cross-checks: `survival::clogit` / `survival::cch` (classical),
  `Epi::ccwc` (risk-set sampling), `multipleNCC` (NCC IPW), `causatr`/`survatr` on
  the explicitly reweighted pseudo-cohort (CCW), R `tmle`/`tmle3` (CCW-TMLE)
- Python cross-language oracles via `statsmodels` for every implemented classical
  estimator (committed data + result CSVs under `tests/testthat/fixtures/python/`,
  compared in `test-python-oracle.R`, `skip_if(!file.exists())`-guarded so CI
  needs no Python). `statsmodels` for the classical MLEs; `delicatessen` is
  reserved for the causal / sandwich estimands of the later CCW phases.
- Update `FEATURE_COVERAGE_MATRIX.md` in the same PR as test changes

## Cost discipline

- **Targeted tests** with `devtools::test(filter = "foo")` during development; full
  `devtools::test()` only before committing (a hook blocks unfiltered runs).
- **Foreground** test/check commands with `timeout: 600000`; never `run_in_background`
  for `devtools::test()` / `check()` (a hook enforces this).
- **Batch R scripts** — combine diagnostics into one `Rscript -e '...'` call.
- **Model awareness** — Sonnet for routine work (formatting, edits, tests, git); Opus
  for variance derivations / subtle debugging / new-feature design.

## Constraints

- Run `devtools::test()` before committing
- Do not modify `man/` or `NAMESPACE` directly
- Run `devtools::document()` after changing roxygen comments

## Scope

matchatr owns the **sampling-design + weighting layer** for case-control-type designs
and the marginal causal contrasts they support. It **delegates** point estimation and
variance to `causatr` (g-comp / IPW / AIPW) and `survatr` (causal survival), to
`survival` (clogit / coxph / cch), `nnet` (polytomous), `multipleNCC` (NCC IPW), and
`survey` (two-phase / calibration). The ONE genuinely new estimator engine is
**CCW-TMLE** (targeting step), because the etverse has no targeted-learning code.
NOT in scope: genetics designs (handbook Ch23-28), measurement-error correction
(Ch10), case-crossover (Ch7, possible future module).

## R ecosystem integration

| Need | Package | Relationship |
|---|---|---|
| Causal engine (g-comp/IPW/AIPW + variance) | `causatr` | **Imports** (delegated) |
| Causal survival (person-period) | `survatr` | **Imports** (delegated) |
| Conditional logistic / weighted Cox / case-cohort | `survival` | **Imports** |
| Sandwich variance | `sandwich` | **Imports** |
| Numerical derivatives | `numDeriv` | **Imports** |
| Bootstrap | `boot` | **Imports** |
| NCC IPW weights + weighted Cox | `multipleNCC` | **Suggests** (engine + oracle) |
| NCC risk-set sampling | `Epi` (`ccwc`) | **Suggests** |
| Two-phase / calibration | `survey` | **Suggests** |
| Multiple imputation (missing-by-design) | `mice` | **Suggests** |
| Firth penalized likelihood | `logistf` | **Suggests** |
| Polytomous logistic | `nnet` | (via `survival`/base) |
| CCW-TMLE oracle | `tmle` | **Suggests** (test-only) |

## Phase roadmap (handbook chapter -> phase)

All `Status: DESIGN`. Reference: *Handbook of Statistical Methods for Case-Control
Studies* (Borgan, Breslow, Chatterjee, Gail, Scott, Wild, 2018).

**Foundation** — `PHASE_1_DESIGN_TAXONOMY` (Ch2: design object, S3, two-step API).

**Classical estimators** — `PHASE_2_UNMATCHED_CC` (Ch3) · `PHASE_3_MATCHED_CC` (Ch4) ·
`PHASE_4_MULTIPLE_GROUPS` (Ch5).

**Time-to-event sampling designs** — `PHASE_5_NESTED_CC` (Ch16,18) ·
`PHASE_6_CASE_COHORT` (Ch16,17) · `PHASE_7_IPW_NCC` (Ch19).

**Causal layer** — `PHASE_8_CAUSAL_STRATEGY` (strategy + Rose & van der Laan) ·
`PHASE_9_CCW_CONTRASTS` (CCW g-formula/IPW/AIPW/TMLE via causatr + new targeting) ·
`PHASE_10_CAUSAL_SURVIVAL_SAMPLED` (design-weighted survatr).

**Efficiency & advanced** — `PHASE_11_TWO_PHASE` (Ch12) · `PHASE_12_CALIBRATION` (Ch13) ·
`PHASE_13_MULTIPLE_IMPUTATION` (Ch20) · `PHASE_14_SEMIPARAMETRIC_MLE` (Ch21) ·
`PHASE_15_SMALL_SAMPLE` (Ch8) · `PHASE_16_POWER` (Ch9) · `PHASE_17_ALT_RISK_MODELS` (Ch11) ·
`PHASE_18_SECONDARY_ANALYSIS` (Ch14).

**Extensions** — `PHASE_19_SCCS` (Ch22) · `PHASE_20_RESPONSE_SELECTIVE` (Ch15).

## Key design decisions (carried from the roadmap)

- **`estimator =` selects the analysis; `design =` selects the sampling structure.**
  Two orthogonal axes. The design object carries strata, time, prevalence q₀, and
  inclusion weights; the estimator decides conditional vs marginal, OR vs HR vs RD/RR.
- **Case-control weights ≠ design weights** (see `hard-rules.md`). Both are observation
  weights into the engine, never data columns.
- **Marginal causal effects reuse causatr/survatr** by passing the q₀ / inclusion
  weights as `weights`; only the sampling-variance correction is matchatr's own.
- **CCW-TMLE is the single new estimator engine** (targeting/fluctuation) — the rest is
  delegation + a weighting/design layer.
- **Conditional likelihood for matched CC and NCC** (CMLE / partial likelihood), never
  unconditional MLE on matched-set indicators.
