# matchatr hard rules

Project-specific rules that override / extend the etverse-wide rules at
`~/Documents/personal/software/etverse/.claude/skills/*/SKILL.md`. Read by the
`implement-feature` and `critical-review-loop` skills before they do anything.

## Project conventions

- **Design-doc pattern.** `PHASE_*.md` at the repo root. When a skill says "read
  the design doc", read the numbered phase doc that covers the feature. All phases
  are currently `Status: DESIGN` (no estimator code shipped yet).
- **Feature coverage file.** `FEATURE_COVERAGE_MATRIX.md`. Every PR that adds,
  removes, or changes a feature MUST update this file.
- **Error-class prefix.** `matchatr_*` for all `rlang::abort()` calls.
- **Repro-script prefix.** `/tmp/matchatr_repro_<slug>.R`.
- **Test-log paths.** `/tmp/matchatr-test-<file>.txt` for per-file runs,
  `/tmp/matchatr-test-results.txt` for the full suite.
- **etverse engine reuse.** matchatr `Imports: causatr, survatr`. Estimation is
  delegated wherever possible; matchatr owns the *sampling-design + weight* layer,
  not a new variance engine. Do not reimplement g-comp / IPW / AIPW / sandwich /
  bootstrap that already live in causatr/survatr.

## Supported dimensions (for combination audits)

| Dimension | Values |
|---|---|
| **Design** | unmatched CC, matched CC, nested CC (NCC), case-cohort, two-phase, counter-matched |
| **Estimator** | conditional logistic (`survival::clogit`), unconditional logistic / Mantel-Haenszel, polytomous (`nnet::multinom`), risk-set / weighted Cox (`survival::coxph`), case-cohort (`survival::cch`: Prentice / Self-Prentice / Borgan I+II), CCW-g-formula, CCW-IPW, CCW-AIPW, CCW-TMLE, design-weighted survatr |
| **Weight type** | none, case-control weights (q₀-based, Rose & van der Laan), design / inclusion-probability weights (Samuelsen KM / GLM / GAM, Borgan), survey / calibration weights |
| **Outcome** | binary (case-control), time-to-event (NCC / case-cohort), polytomous (multiple case groups) |
| **Estimand** | conditional OR, conditional HR, marginal RD, marginal RR, marginal OR, absolute risk F(t), RMST |
| **Contrast** | OR, HR, RD, RR, marginal OR |
| **Variance** | model-based (information matrix), robust sandwich (Lin-Wei), Self-Prentice / Samuelsen asymptotic, CCW influence function, bootstrap |
| **Small-sample** | none, Firth penalized (`logistf`), exact / mid-p |
| **Missing data** | complete, missing-by-design (MI via `mice`) |

## Hard rules (appended to the skill's generic rules)

### Architecture invariants — DO NOT flag these as bugs without a numerical reproducer

- **Case-control weights and design weights are DISTINCT objects with distinct
  variance consequences.** Case-control weights (Rose & van der Laan) reweight a
  CC sample back to the source population using the marginal prevalence q₀.
  Design / inclusion-probability weights (Samuelsen, Borgan) reweight a sampled
  cohort by the inverse probability of being selected as a control / subcohort
  member. They answer different questions and their sandwich corrections differ.
  Never conflate them or share a code path that assumes one variance form.
- **Weights are observation weights into causatr/survatr fits, NEVER a data
  column.** Mirror causatr's invariant ("weights live in `fit$details$weights`").
  Both CC weights and design weights enter the engine as the `weights` argument;
  the design object carries them, the data frame does not.
- **Matched case-control requires CONDITIONAL likelihood (CMLE), not
  unconditional MLE.** Fitting matched-set indicators by ordinary logistic
  regression biases the OR (for 1:1 matching the MLE → OR² in large samples;
  Pike et al. 1980, Breslow & Day 1980). matchatr fits matched CC via
  `survival::clogit` (conditional partial likelihood). The unconditional path is
  only valid for few observable confounders adjusted as covariates.
- **NCC conditional partial likelihood = sampled risk set as a stratum.** The
  classical NCC analysis (Thomas 1977; Ch16, 18) is `survival::clogit`/`coxph`
  with each sampled risk set (case + its m controls) as a stratum. This is the
  same conditional likelihood form as matched CC. Do not "simplify" it to an
  unstratified model.
- **OR = HR exactly under proper risk-set (incidence-density) matching** — no
  rare-disease assumption needed (Miettinen 1976; Prentice & Breslow 1978). Do
  not add a rare-disease caveat to NCC/risk-set-matched contrasts.
- **causatr has NO targeted-learning machinery.** It ships g-comp / IPW / AIPW
  only. CCW-g-formula, CCW-IPW, and CCW-AIPW reuse causatr directly via the
  `weights` argument. CCW-**TMLE** (the targeting / fluctuation step) is genuinely
  NEW code that matchatr must implement itself — it is not a causatr reuse. Do not
  assume a TMLE entry point exists anywhere in the etverse.
- **Case-cohort pseudo-likelihood is NOT a true likelihood.** Prentice / Borgan
  estimators reuse controls across failure times, so the score factors are
  dependent: standard errors do NOT come from the information matrix and LR
  statistics are not χ². Use the Self-Prentice / Borgan asymptotic variance or a
  robust sandwich (conservative). Delegate to `survival::cch`, which implements
  these correctly; do not hand-roll the naive information-matrix SE.
- **Robust sandwich is conservative for stratified case-cohort (Borgan II).**
  Ch16.4.4 warns it substantially overestimates the variance there; the plug-in
  asymptotic estimator is preferred. Do not flag the plug-in vs sandwich gap as a
  bug.

### Invariants to enforce in code (tests must exercise, not flag)

- **Reject non-binary outcome for case-control designs** with a classed error
  (`matchatr_bad_outcome`).
- **Reject matching/conditional estimators with strata containing 0 cases or
  0 controls** (uninformative; `survival::clogit` drops them) — warn with a
  classed condition (`matchatr_uninformative_stratum`).
- **q₀ (prevalence) is required for case-control weighting.** Either supplied as
  known, or estimated from the full cohort. Missing q₀ for a CCW estimator is a
  classed error (`matchatr_missing_prevalence`).
- **The case indicator must be a genuine binary with both classes present.**
  `resolve_binary_outcome()` accepts logical / two-level factor / numeric 0/1,
  coerces to 0/1, then requires both a case and a control to occur — an
  all-cases / all-controls / all-NA column is `matchatr_bad_outcome` for *every*
  encoding (the check is uniform, not numeric-only).
- **A column may hold only one analysis role.** The outcome and exposure must
  each be distinct from the confounders and from the design-referenced columns
  (strata/time/subcohort/phase) — `matchatr_bad_input`. BUT confounders and
  design columns MAY overlap: frequency-matching on a variable while also
  adjusting for residual confounding by it is valid — do NOT flag that overlap.

### Review-time heuristics

- **Before flagging a CCW variance issue**, derive whether the q₀ weights are
  treated as fixed (known prevalence) or estimated (cohort-estimated) — the IF
  has an extra term in the estimated-q₀ case. Run sandwich vs bootstrap.
- **Cross-check NCC/case-cohort against `multipleNCC` and `survival::cch`**, and
  CCW marginal contrasts against `causatr` run on the explicitly reweighted
  pseudo-cohort, before claiming a numerical bug.
- **`survival`, `multipleNCC`, `survey`, `Epi` are oracles / delegated engines**,
  not things to reimplement.
