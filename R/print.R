#' Print a matchatr design
#'
#' @description
#' Displays a compact summary of a [matchatr_design][unmatched_cc] object: the
#' sampling structure, the columns it references, the matching ratio and
#' prevalence when present, and the weighting scheme it will use.
#'
#' @param x A `matchatr_design` object.
#' @param ... Unused; present for S3 consistency.
#' @returns Invisibly returns `x`.
#' @examples
#' print(matched_cc(strata = "set", ratio = 2))
#' print(unmatched_cc(prevalence = 0.02))
#' @seealso [unmatched_cc()], [matched_cc()], [nested_cc()], [case_cohort()]
#' @export
print.matchatr_design <- function(x, ...) {
  cat("<matchatr_design>\n")
  cat(" Type:       ", design_label(x$type), "\n", sep = "")
  if (!is.null(x$strata)) {
    cat(" Strata:     ", paste(x$strata, collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$time)) {
    cat(" Time:       ", x$time, "\n", sep = "")
  }
  if (!is.null(x$subcohort)) {
    cat(" Subcohort:  ", x$subcohort, "\n", sep = "")
  }
  if (!is.null(x$phase1)) {
    cat(" Phase 1:    ", paste(x$phase1, collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$phase2)) {
    cat(" Phase 2:    ", x$phase2, "\n", sep = "")
  }
  if (!is.null(x$ratio)) {
    cat(" Ratio:      ", x$ratio, ":1\n", sep = "")
  }
  if (!is.null(x$prevalence)) {
    cat(" Prevalence: ", x$prevalence, " (q0)\n", sep = "")
  }
  cat(" Weights:    ", weight_label(x$weight_spec$kind), "\n", sep = "")
  invisible(x)
}

#' Print a matchatr fit
#'
#' @description
#' Displays a compact summary of a [matchatr_fit][matcha] object: the sampling
#' design, the resolved estimator and engine, the outcome / exposure /
#' confounder roles, and the case / control counts.
#'
#' @param x A `matchatr_fit` object.
#' @param ... Unused; present for S3 consistency.
#' @returns Invisibly returns `x`.
#' @examples
#' df <- data.frame(case = c(1, 0, 1, 0), x = c(1, 0, 1, 0))
#' print(matcha(df, outcome = "case", exposure = "x", design = unmatched_cc()))
#' @seealso [matcha()]
#' @export
print.matchatr_fit <- function(x, ...) {
  cat("<matchatr_fit>\n")
  cat(" Design:     ", design_label(x$design$type), "\n", sep = "")
  cat(
    " Estimator:  ",
    x$estimator,
    "  (engine: ",
    x$engine,
    ")\n",
    sep = ""
  )
  cat(" Outcome:    ", x$outcome, "\n", sep = "")
  cat(" Exposure:   ", x$exposure, "\n", sep = "")
  conf_label <- if (is.null(x$confounders)) {
    "none"
  } else {
    deparse1(x$confounders)
  }
  cat(" Confounders: ", conf_label, "\n", sep = "")
  if (!is.null(x$effect_modifier)) {
    cat(" Modifier:   ", x$effect_modifier, "\n", sep = "")
  }
  if (identical(x$details$outcome_kind, "polytomous")) {
    # A multi-group outcome has no single case / control split; report the
    # per-group counts with the reference (baseline) group flagged.
    counts <- x$details$group_counts
    labels <- vapply(
      names(counts),
      function(g) {
        if (identical(g, x$details$reference)) paste0(g, "*") else g
      },
      character(1)
    )
    cat(
      " N:          ",
      nrow(x$data),
      "  (groups: ",
      paste0(labels, ": ", as.integer(counts), collapse = ", "),
      ")\n",
      sep = ""
    )
    cat("             * reference group\n", sep = "")
  } else {
    cat(
      " N:          ",
      nrow(x$data),
      "  (cases: ",
      x$details$n_cases,
      ", controls: ",
      x$details$n_controls,
      ")\n",
      sep = ""
    )
  }
  invisible(x)
}

#' Print a matchatr contrast result
#'
#' @description
#' Displays a compact summary of a `matchatr_result`: the
#' estimator and engine that produced it, the estimand and contrast scale, the
#' confidence-interval method and sample size, followed by the contrasts table.
#'
#' @param x A `matchatr_result` object returned by [contrast()].
#' @param ... Unused; present for S3 consistency.
#' @returns Invisibly returns `x`.
#' @examples
#' set.seed(1)
#' df <- data.frame(case = rep(c(1, 0), each = 100), x = rbinom(200, 1, 0.4))
#' fit <- matcha(df, outcome = "case", exposure = "x", design = unmatched_cc())
#' print(contrast(fit, type = "or"))
#' @seealso [contrast()]
#' @export
print.matchatr_result <- function(x, ...) {
  cat("<matchatr_result>\n")
  cat(
    " Estimator:  ",
    x$estimator,
    "  (engine: ",
    x$engine,
    ")\n",
    sep = ""
  )
  cat(" Estimand:   ", x$estimand, "\n", sep = "")
  cat(" Contrast:   ", contrast_label(x$type), "\n", sep = "")
  cat(" CI method:  ", x$ci_method, "\n", sep = "")
  cat(" N:          ", x$n, "\n", sep = "")
  cat("\nContrasts:\n")
  print(x$contrasts)
  invisible(x)
}

#' Human-readable label for a contrast scale
#'
#' @param type Character scalar contrast type.
#' @returns A character scalar label; the raw `type` if unrecognised so the
#'   print method never errors on a malformed object.
#' @family printing
#' @noRd
contrast_label <- function(type) {
  switch(
    type,
    difference = "Risk difference",
    ratio = "Risk ratio",
    or = "Odds ratio",
    hr = "Hazard ratio",
    type
  )
}

#' Human-readable label for a design type
#'
#' @param type Character scalar design type.
#' @returns A character scalar label; the raw `type` if unrecognised so the
#'   print methods never error on a malformed object.
#' @family printing
#' @noRd
design_label <- function(type) {
  switch(
    type,
    unmatched_cc = "Unmatched case-control",
    matched_cc = "Matched case-control",
    nested_cc = "Nested case-control",
    case_cohort = "Case-cohort",
    two_phase = "Two-phase",
    counter_matched = "Counter-matched",
    type
  )
}

#' Human-readable label for a weighting scheme
#'
#' @param kind Character scalar weight-spec kind.
#' @returns A character scalar label; the raw `kind` if unrecognised.
#' @family printing
#' @noRd
weight_label <- function(kind) {
  switch(
    kind,
    none = "none",
    case_control = "case-control (q0)",
    inclusion = "inclusion-probability",
    design = "design (two-phase)",
    counter_match = "counter-matching",
    kind
  )
}
