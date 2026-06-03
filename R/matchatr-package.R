#' matchatr: Causal Inference for Case-Control, Nested Case-Control, and
#' Case-Cohort Studies
#'
#' Classical and causal estimation for (matched) case-control, nested
#' case-control, and case-cohort study designs. matchatr supplies the
#' sampling-design and weighting layer — case-control weights (the Rose and
#' van der Laan g-formula / IPW / AIPW / TMLE family) and design-based
#' inclusion weights (Samuelsen, Borgan) — and delegates point estimation and
#' variance to the etverse engines \pkg{causatr} (g-computation / IPW / AIPW)
#' and \pkg{survatr} (causal survival on person-period data), and to
#' \pkg{survival} (conditional logistic, weighted Cox, case-cohort).
#'
#' The package is under active development; see the `PHASE_*.md` design docs and
#' `CLAUDE.md` at the repository root for the implementation roadmap.
#'
#' @keywords internal
#' @importFrom survival clogit coxph strata Surv
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
