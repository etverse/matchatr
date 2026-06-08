"""Independent oracle for matchatr's polytomous (multinomial) subtype ORs.

Reads ``polytomous_or_data.csv`` (a three-group outcome ``g`` with levels
control / caseA / caseB, binary exposure ``x``, continuous confounder ``age``),
fits the baseline-category multinomial logistic model ``g ~ x + age`` with
``control`` as the reference, and writes ``polytomous_or_results.csv``: for each
non-reference subtype, the exposure log-odds-ratio versus the reference, its SE,
and the OR with a 95% Wald interval (log scale, exponentiated).

Run from this directory:  python3 polytomous_or.py
"""

import os

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> None:
    data = pd.read_csv(os.path.join(HERE, "polytomous_or_data.csv"))
    # Encode the outcome so "control" is code 0 -- MNLogit takes the smallest
    # code as the baseline category and contrasts the others against it. The
    # formula interface cannot be used here: it would expand the 3-level factor
    # into a multi-column endog, which MNLogit rejects.
    g = pd.Categorical(data["g"], categories=["control", "caseA", "caseB"])
    endog = g.codes  # control=0 (baseline), caseA=1, caseB=2
    exog = sm.add_constant(data[["x", "age"]])
    fit = sm.MNLogit(endog, exog).fit(disp=0)
    z = stats.norm.ppf(0.975)
    # params / bse are (predictor x non-baseline-equation) frames; columns are
    # the caseA, caseB equations in code order. Pull the exposure (x) row.
    levels = ["caseA", "caseB"]
    rows = []
    for j, lev in enumerate(levels):
        b = fit.params.loc["x"].iloc[j]
        se = fit.bse.loc["x"].iloc[j]
        rows.append(
            {
                "y_level": lev,
                "term": "x",
                "estimate": b,  # subtype log odds ratio vs control
                "std_error": se,
                "odds_ratio": np.exp(b),
                "conf_low": np.exp(b - z * se),
                "conf_high": np.exp(b + z * se),
            }
        )
    pd.DataFrame(rows).to_csv(
        os.path.join(HERE, "polytomous_or_results.csv"), index=False
    )


if __name__ == "__main__":
    main()
