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
- **One conditional scale per design (the `clogit` engine).** The conditional
  partial likelihood gives `exp(beta)`, whose meaning is fixed by the sampling
  design: a matched case-control design reports the conditional odds ratio
  (`type = "or"`), a nested case-control (risk-set-sampled) design the hazard
  ratio (`type = "hr"`). Both resolve to the same `clogit` engine, so
  `default_contrast_type()` is design-aware and `contrast_clogit()` rejects the
  off-design scale (an OR asked of a risk-set design, or an HR asked of a matched
  design) with `matchatr_unidentified_estimand` via
  `reject_offdesign_conditional_scale()`. The number would be identical, but the
  estimand the design *targets* is not, and reporting the wrong label would
  mislead. Do NOT flag this rejection as overly strict, do NOT "generalise" the
  engine to emit both scales, and do NOT relabel the NCC result as an odds ratio.
  The design's `time` column records how controls were sampled; the conditional
  likelihood reads the risk set from `strata` and does not enter `time` (it feeds
  the later inclusion-weight / weighted-Cox designs).
- **causatr has NO targeted-learning machinery.** It ships g-comp / IPW / AIPW
  only. CCW-g-formula, CCW-IPW, and CCW-AIPW reuse causatr directly via the
  `weights` argument. CCW-**TMLE** (the targeting / fluctuation step) is genuinely
  NEW code that matchatr must implement itself — it is not a causatr reuse. Do not
  assume a TMLE entry point exists anywhere in the etverse.
- **A CCW (`ccw_*`) fit's `model` is a `causatr_fit`, which has NO `coef()` /
  `vcov()` method.** It reports a *marginal* effect whose scale is chosen at the
  contrast step (RD / RR / marginal OR), not a single conditional coefficient
  table, so `tidy.matchatr_fit` and `summary.matchatr_fit` deliberately branch on
  `inherits(model, "causatr_fit")` and surface the marginal contrast (the risk
  difference, the default ccw scale) instead of building an odds-ratio table via
  `estimable_vcov()`. Do NOT route a ccw fit through the `coef()`/`vcov()` path,
  and do NOT flag `tidy(fit)` returning the marginal contrast (rather than a
  coefficient table) as inconsistent — it is the analogue of the conditional-OR
  table the other engines report.
- **A CCW marginal contrast records `ci_method = "sandwich"` regardless of the
  requested label.** A marginal g-formula contrast has no information-matrix
  variance distinct from causatr's influence-function / sandwich one, so
  `ci_method = "model"` and `"sandwich"` both yield it and `contrast_ccw()`
  records what causatr actually computed (`"sandwich"`). `"bootstrap"` is the
  within-stratum percentile bootstrap (`ccw_bootstrap_ci()`, `R/variance_ccw.R`):
  it resamples cases and controls separately so the design (n1 / n0) is preserved,
  which keeps the q0 weights constant across replicates (known q0 is fixed; its
  sampling variability is a separate estimated-q0 IF term), refits the engine per
  replicate, and reports the percentile interval while keeping the analytic point
  estimate. Do NOT "fix" it to resample the whole sample (that mixes the strata
  and breaks the design), and do NOT flag the constant per-replicate weights as a
  missing reweighting. `fit_ccw()` fits the outcome model with `family = "quasibinomial"`
  (the right family for fractional case-control weights — identical mean model
  and sandwich to binomial, but silent on the spurious `non-integer #successes`
  warning a binomial fit raises), so do NOT switch it back to `"binomial"`. The
  single `family` argument governs the weighted outcome / marginal-mean fit for
  all three estimators (g-formula and AIPW outcome models, IPW's weighted marginal
  mean); causatr auto-detects the propensity family from the binary treatment, so
  passing it on the IPW path is correct, not a stray argument. `fit_ccw()` is
  parameterized over `fit$estimator` (`ccw_gformula`/`ccw_ipw`/`ccw_aipw` →
  `causat(estimator = "gcomp"/"ipw"/"aipw")`) and all three require `confounders`
  (the outcome / propensity / both adjustment models). Do NOT flag the `"model"` →
  `"sandwich"` relabeling or the confounders requirement as bugs.
- **An estimated q0 (`unmatched_cc(prevalence = q0, prevalence_n = N)`) widens the
  CCW interval; the point estimate is unchanged.** q0 estimated from N cohort
  members carries sampling uncertainty Var(q̂0) = q0(1−q0)/N into ψ̂ through the
  weights. The analytic (sandwich / EIF) path adds the delta-method term
  (∂ψ/∂q0)²·Var(q̂0) — `ccw_estimated_q0_term()` computes ∂ψ/∂q0 by a central
  finite difference (refit at q0±h) on the reported scale (linear for RD, log for
  RR/OR) via `ccw_apply_estimated_q0()`; the bootstrap path redraws
  q0* ~ Binomial(N, q0)/N per replicate. The two agree and both collapse onto the
  known-q0 interval as N → ∞. **Critical invariant:** `ccw_boot_point()` (the
  per-replicate / finite-difference point) MUST strip `prevalence_n` before
  computing the contrast, because the estimated-q0 and bootstrap variance branches
  both call back into it — leaving `prevalence_n` set causes runaway recursion. Do
  NOT remove that `fit$design$prevalence_n <- NULL`, and do NOT compute the
  estimated-q0 term for a known q0.
- **The CCW estimators run on `unmatched_cc` and `matched_cc`, never on
  `nested_cc`.** Case-control weighting maps a case-control sample to the source
  cohort via the q0 case/control reweighting; a matched CC is still a case-control
  sample (the matching variable is a baseline covariate, NOT a conditioning
  stratum — it must be in `confounders` so the marginal effect is standardized over
  its distribution; matched sets are ignored, Rose & van der Laan 2009). A **nested**
  CC is risk-set / incidence-density sampled, so its controls are not a case-control
  sample and binary q0 reweighting does not identify a marginal estimand —
  `matcha(design = nested_cc(...), estimator = "ccw_*")` aborts `matchatr_bad_estimator`
  toward `ipw_cox`, and this guard fires in `matcha()` BEFORE the missing-prevalence
  check. Do NOT add a `prevalence` arg to `nested_cc()`, and do NOT condition a CCW
  fit on the matched `set` column.
- **The whole CCW family complete-cases missing data once, in `ccw_prepare()`.**
  Rows with a missing outcome, exposure, or confounder are dropped (listwise
  deletion) with a classed `matchatr_dropped_rows` warning, and `cc_weights()` is
  computed on the complete-case sample so the weighted case fraction still equals
  q0 (do NOT compute the weights on the full sample then drop — that breaks the
  Rose–van der Laan mapping). This is deliberate and matches matchatr's classical
  engines and causatr's default `na.action`: it unifies what was inconsistent
  (CCW-g-formula complete-cased via causatr, CCW-IPW/AIPW errored on a confounder
  NA, the hand-rolled CCW-TMLE crashed on misaligned vectors). Do NOT move the
  complete-casing into each engine, do NOT switch CCW to a reject-on-NA error, and
  do NOT flag the listwise deletion as a bug. Multiple imputation (missing
  confounders) and an outcome-missingness / IPCW extended-TMLE (missing outcomes)
  are the deferred principled alternatives — PHASE_13, not this layer. A
  non-converging CCW-TMLE fluctuation warns `matchatr_tmle_convergence` and reverts
  to the untargeted initial fit (a defensive NaN guard, not a bug).
- **The CCW double-robustness test asserts recovery with `expect_lt(abs(est -
  truth), BAND)`, NOT `expect_equal`.** The marginal risk difference is
  small-magnitude (~0.05), below the level at which `all.equal()` / waldo switch
  from a relative to an absolute tolerance, which makes `expect_equal(est, truth,
  tolerance =)` behave unpredictably (a loose absolute band that a biased estimator
  also passes, or a too-tight relative one the consistent estimator fails). The
  absolute-error band `abs(est - truth) < BAND` IS a two-sided check against the
  analytical oracle (the g-formula truth), and the DGP gives a clean ~10x gap
  (consistent estimators within ~0.003, the misspecified one outside ~0.03), so
  `BAND = 0.01` separates them robustly. This is the same idiom causatr's own DR
  tests use (`expect_lt(abs(est - truth), tol)`). Do NOT flag this as a
  forbidden `expect_lt`-on-a-point-estimate — the ban is on direction tests
  (`expect_gt(est, 0)`), not on absolute-error bands around a known truth.
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
- **Never positionally index an SE vector against `coef()`.** `stats::vcov()`
  keeps aliased (rank-deficient) coefficients as `NA` rows, but
  `sandwich::sandwich()` *drops* them, so `sqrt(diag(.))` lengths differ between
  the two variance sources and recycling silently corrupts SEs. Always reduce
  variance to the estimable coefficients and align by NAME (see
  `estimable_vcov()`); aliased terms get `NA` SEs. A non-estimable exposure
  (constant / collinear) is `matchatr_unestimable_exposure`.
- **The McNemar estimator (`estimator = "mcnemar"`) is the 1:1 binary-exposure
  special case ONLY.** It computes OR = n10/n01 with Var(log OR) = 1/n10 + 1/n01
  in closed form (no `survival::clogit`) over the discordant pairs; pairs
  concordant on exposure cancel. A matched set with more than one case or more
  than one control is M:1 (or richer) matching, which has no two-cell closed
  form — it is rejected with `matchatr_not_one_to_one` and rerouted to
  `estimator = "clogit"`. Do NOT flag this rejection as overly strict, and do
  NOT "generalise" McNemar to M:1. A one-sided / empty set of discordant pairs
  (n10 = 0 or n01 = 0) is a boundary OR of 0 / Inf and is
  `matchatr_unestimable_exposure`, not a silent 0 / Inf. The case/control
  exposures are aligned per set by `order()` on the droplevels'd stratum factor;
  this is correct under scrambled row order and non-sequential set ids (verified
  against `clogit`), so do not flag the `order()`-index alignment as a bug.
- **The odds-ratio interval is Wald on the log scale, exponentiated** (asymmetric
  on the OR scale) — this is correct, not a bug. The OR-scale `se` in a result's
  `contrasts` is the delta-method `OR * SE(log OR)`; it does NOT reconstruct the
  CI. The reconstructable log-scale SE lives in the result's `estimates`. Do not
  "fix" the asymmetry or flag the se↔CI mismatch.
- **`nnet::multinom` does NOT alias a rank-deficient predictor to `NA`** the way
  `stats::glm` does — it splits the coefficient across the collinear columns and
  returns a finite, silently attenuated odds ratio (verified: a confounder
  `dup = x` halves the exposure log OR). The polytomous engine guards this at fit
  time in `reject_collinear_exposure()` (`R/polytomous.R`) by the design-matrix
  rank: it aborts `matchatr_unestimable_exposure` only when dropping the exposure
  column(s) fails to lower `qr(model.matrix)$rank` by the number of exposure
  columns. Collinearity confined to the confounders (the exposure still adds full
  rank) is deliberately NOT rejected — the exposure OR is still identified. Do
  not replace this with a raw constant-column check (it misses confounder
  collinearity), do not rely on an `anyNA(coef)` guard (multinom never returns
  `NA`, so that branch is dead), and do not flag the non-rejection of
  confounder-only collinearity as a bug.

### Review-time heuristics

- **Before flagging a CCW variance issue**, derive whether the q₀ weights are
  treated as fixed (known prevalence) or estimated (cohort-estimated) — the IF
  has an extra term in the estimated-q₀ case. Run sandwich vs bootstrap.
- **Cross-check NCC/case-cohort against `multipleNCC` and `survival::cch`**, and
  CCW marginal contrasts against `causatr` run on the explicitly reweighted
  pseudo-cohort, before claiming a numerical bug.
- **`survival`, `multipleNCC`, `survey`, `Epi` are oracles / delegated engines**,
  not things to reimplement.
- **`fit_ipw_cox()` fits the weighted Cox with `ties = "breslow"` (NOT the coxph
  default Efron).** This is deliberate: the IPW Breslow cumulative baseline hazard
  used for absolute risk (`ipw_breslow_ncc()`) is a plain Breslow step, which is
  inconsistent with Efron coefficients at tied event times. With `ties = "breslow"`
  the partial-likelihood β and the Breslow baseline agree (verified to machine
  precision against `survival::survfit` under heavy ties). Under incidence-density
  NCC sampling failure times are typically distinct, so it equals Efron there. Do
  NOT "restore" the Efron default or flag this as a non-standard choice.
- **Absolute-risk linear predictors are built from the FITTED model's terms, not
  a re-derived formula.** `ar_lp_from_newdata()` builds the design for the
  `ipw_cox`/coxph engine via `model.matrix(delete.response(terms(model)),
  model.frame(..., xlev = model$xlevels))` so a data-dependent confounder basis
  (`poly`/`ns`/`bs`/`scale`) is reproduced from the fit's `predvars`, not
  recomputed from `newdata` (which would give a different basis with identical
  coefficient names — a silent, catastrophic LP error). The `cch` engine
  deliberately uses the original-formula path instead, because `survival::cch`'s
  non-standard internal formula has no `predvars` map aligned to its coefficient
  names; data-dependent transforms are therefore not reproduced for `cch` (its
  designs do not use them). Do NOT unify these two paths.
- **For a `clogit`/`coxph` fit, `model$n` is the rows used, `nobs()` is the
  event count.** The analysis size and the missing-data count
  (`n_dropped = nrow(data) - model$n`) must read `model$n`, never `nobs()`.
  `model$n` still counts the rows of an *uninformative* stratum that `clogit`
  drops from the likelihood, so a dropped stratum does NOT inflate `n_dropped`
  and triggers no `matchatr_dropped_rows` warning — do not "fix" this.
- **For a `survival::cch` fit, `model$n` ≠ nrow(subset).** `cch$n` equals
  `subcohort_size + n_events`, double-counting subjects who are both subcohort
  members and cases. `nrow(subset) - cch$n` is therefore always negative and
  cannot detect NA-dropped rows. Do NOT attempt `n_dropped = nrow(data) - cch$n`;
  NA handling is delegated to `survival::cch` via its `na.action`. There is no
  simple `model$n`-based NA-drop check for `cch` fits.
- **The only oracle that validates the matched-CC conditional VARIANCE
  independently of `survival::clogit` is the 1:1 McNemar closed form**:
  OR = n10/n01, Var(log OR) = 1/n10 + 1/n01 over the discordant pairs. Comparing
  a clogit wrapper's SE to `clogit`'s own `vcov()` only checks forwarding. A
  truth-DGP recovery test must use a SE-scaled band (the estimator's sampling SD
  ≈ the reported SE), not a fixed absolute tolerance below one SD.
