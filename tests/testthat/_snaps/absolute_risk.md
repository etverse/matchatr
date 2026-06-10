# absolute_risk snapshots: rejection error messages

    Code
      absolute_risk(fit_fake, newdata = data.frame(x = 1), times = 1)
    Condition
      Error in `absolute_risk()`:
      ! `absolute_risk()` is not implemented for the `clogit` engine.
      i Supported engines: the case-cohort (`cch`), IPW nested case-control weighted Cox (`ipw_cox`), and IPW nested case-control Weibull AFT (`ipw_aft`) survival fits.

