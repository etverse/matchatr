# Epi::ccwc risk-set sampling oracle for sample_ncc(). Epi::ccwc ("cohort
# case-control, weighted complement") performs the same incidence-density
# control sampling, so it is an independent external cross-check of HOW the risk
# set (the eligible control pool) is defined -- not an exact-draw oracle, since
# the two implementations consume separate random streams. Guarded by
# skip_if_not_installed("Epi") at every call site.
#
# Returns one row per sampled subject with standardised columns: `set` (matched-
# set id), `case` (1 for the failing subject, 0 for a control), `risk_time` (the
# set's failure time), and `row` (the index of the member in the input cohort, so
# the test can map a sampled control back to its cohort row). Epi's `Map` column
# is exactly that input-row index.
epi_ccwc_riskset <- function(cohort, time, event, m) {
  # Epi::ccwc evaluates `exit` / `fail` non-standardly inside `data`
  # (eval(substitute(exit), data)), so the string-named columns are copied to
  # fixed literal names the call can reference as bare symbols.
  d <- as.data.frame(cohort)
  d$.ccwc_exit <- d[[time]]
  d$.ccwc_fail <- d[[event]]
  cc <- Epi::ccwc(
    exit = .ccwc_exit,
    fail = .ccwc_fail,
    controls = m,
    data = d,
    silent = TRUE
  )
  data.frame(
    set = as.integer(cc$Set),
    case = as.integer(cc$Fail),
    risk_time = as.numeric(cc$Time),
    row = as.integer(cc$Map),
    stringsAsFactors = FALSE
  )
}
