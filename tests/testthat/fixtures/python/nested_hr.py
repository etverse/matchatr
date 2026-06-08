"""Independent oracle for matchatr's nested case-control hazard ratio.

Reads ``nested_hr_data.csv`` (a risk-set-sampled NCC dataset: per-set case
indicator, exposure ``x``, confounder ``z``, sampled-risk-set id ``set``) and
fits the conditional partial likelihood with each sampled risk set as a stratum
-- the same conditional logistic likelihood as the matched design. Under
risk-set (incidence-density) sampling exp(beta) IS the hazard ratio (OR = HR
exactly; Prentice & Breslow 1978), so the reported quantity is labelled a hazard
ratio. Writes ``nested_hr_results.csv``: per term, the log-HR, its SE, and the
HR with a 95% Wald interval (log scale, exponentiated).

Run from this directory:  python3 nested_hr.py
"""

import os

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> None:
    data = pd.read_csv(os.path.join(HERE, "nested_hr_data.csv"))
    # Same conditional partial likelihood as the matched design; the sampled risk
    # set is the stratum. No intercept (conditioned out).
    fit = sm.ConditionalLogit(
        data["case"], data[["x", "z"]], groups=data["set"]
    ).fit(disp=0)
    z = stats.norm.ppf(0.975)
    beta = fit.params
    se = fit.bse
    out = pd.DataFrame(
        {
            "term": beta.index,
            "estimate": beta.to_numpy(),  # log hazard ratio (= log OR under risk-set sampling)
            "std_error": se.to_numpy(),
            "hazard_ratio": np.exp(beta.to_numpy()),
            "conf_low": np.exp(beta.to_numpy() - z * se.to_numpy()),
            "conf_high": np.exp(beta.to_numpy() + z * se.to_numpy()),
        }
    )
    out.to_csv(os.path.join(HERE, "nested_hr_results.csv"), index=False)


if __name__ == "__main__":
    main()
