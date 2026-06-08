"""Independent oracle for matchatr's Mantel-Haenszel summary odds ratio.

Reads ``mh_or_data.csv`` (a stratified case-control sample: case indicator,
binary exposure ``x``, stratum ``agegrp``), forms the per-stratum 2x2 exposure x
case table, and computes the Mantel-Haenszel pooled odds ratio with the
Robins-Breslow-Greenland confidence interval via statsmodels'
``StratifiedTable`` -- the same pooled estimator and RBG variance matchatr uses.
Writes ``mh_or_results.csv``: the pooled OR and its 95% interval.

Run from this directory:  python3 mh_or.py
"""

import os

import numpy as np
import pandas as pd
import statsmodels.api as sm

HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> None:
    data = pd.read_csv(os.path.join(HERE, "mh_or_data.csv"))
    # One 2x2 table per stratum, rows = exposure (1, 0), cols = case (1, 0), so
    # the table OR is (a*d)/(b*c) = (exposed-case * unexposed-control) /
    # (exposed-control * unexposed-case).
    tables = []
    for _, sub in data.groupby("agegrp"):
        a = int(((sub.x == 1) & (sub.case == 1)).sum())
        b = int(((sub.x == 1) & (sub.case == 0)).sum())
        c = int(((sub.x == 0) & (sub.case == 1)).sum())
        d = int(((sub.x == 0) & (sub.case == 0)).sum())
        tables.append(np.array([[a, b], [c, d]]))
    st = sm.stats.StratifiedTable(tables)
    lo, hi = st.oddsratio_pooled_confint()  # Robins-Breslow-Greenland
    out = pd.DataFrame(
        {
            "term": ["x"],
            "odds_ratio": [st.oddsratio_pooled],
            "conf_low": [lo],
            "conf_high": [hi],
        }
    )
    out.to_csv(os.path.join(HERE, "mh_or_results.csv"), index=False)


if __name__ == "__main__":
    main()
