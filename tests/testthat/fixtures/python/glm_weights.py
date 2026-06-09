"""
Cross-language oracle for matchatr's GLM working-model NCC inclusion probabilities.

Implements the Borgan, Samuelsen & Aastveit (2003) logistic working model for
the probability that a cohort subject is selected as a control at each event
time.  The steps mirror matchatr's compute_ncc_weights(method = "glm"):

1. Build the augmented selection dataset: one row per (eligible j, event k)
   pair where j is in the risk set at t_k and is not the case.
2. Fit logistic regression: selected ~ risk_time (time-only working model).
3. Apply the product formula:
     pi_j = 1 - prod_{k: j in R(t_k)} (1 - p_hat_jk)
4. Output the ipw_weight (1/pi_j) for each unique cohort subject in the NCC.

The results are compared against matchatr's compute_ncc_weights() in
test-ipw_ncc.R (tolerance 1e-4 in ipw_weight space, with cases forced to 1).

Usage:
    cd tests/testthat/fixtures/python
    python3 glm_weights.py

Requires: numpy, pandas, statsmodels (>= 0.14)
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm

# ── load shared data ──────────────────────────────────────────────────────────
cohort = pd.read_csv("glm_weights_cohort.csv")
ncc    = pd.read_csv("glm_weights_ncc.csv")

n_cohort = len(cohort)

# ── build augmented selection dataset ────────────────────────────────────────
# For each event time t_k (identified by each unique (set, risk_time) pair):
#   - find all cohort subjects at risk at t_k: cohort.t >= t_k, not the case
#   - mark which were selected as controls
events = (
    ncc[["set", "risk_time"]]
    .drop_duplicates()
    .sort_values("set")
    .reset_index(drop=True)
)

records = []
for _, evt in events.iterrows():
    risk_time_k = float(evt["risk_time"])
    set_k       = int(evt["set"])

    set_mask   = ncc["set"] == set_k
    # .cohort_row in R is 1-indexed; convert to 0-indexed for Python arrays.
    case_row_k = int(ncc.loc[set_mask & (ncc["case"] == 1), ".cohort_row"].values[0]) - 1
    ctrl_rows  = set(
        ncc.loc[set_mask & (ncc["case"] == 0), ".cohort_row"].values - 1
    )

    # Eligible pool: at risk at t_k, not the case (same logic as R eligible_controls).
    elig = np.where(cohort["t"].values >= risk_time_k)[0]
    elig = elig[elig != case_row_k]

    for j in elig:
        records.append({
            ".cohort_row": j + 1,          # back to 1-indexed to match R output
            "risk_time":   risk_time_k,
            "selected":    1 if j in ctrl_rows else 0,
        })

aug = pd.DataFrame(records)

# ── fit logistic selection model: selected ~ risk_time ────────────────────────
# Exactly the default selection_formula = ~ risk_time used by matchatr.
X     = sm.add_constant(aug["risk_time"].values.astype(float))
model = sm.Logit(aug["selected"].values.astype(float), X).fit(disp=False)
aug["p_hat"] = model.predict(X)

# ── compute per-subject inclusion probabilities (product formula) ─────────────
# pi_j = 1 - prod_{k: j in R(t_k)} (1 - p_hat_jk)
#       = -expm1(sum_{k: j in R(t_k)} log1p(-p_hat_jk))
p_clamped  = np.minimum(aug["p_hat"].values, 1.0 - np.finfo(float).eps)
log_contrib = np.log1p(-p_clamped)

log_surv = np.zeros(n_cohort)
for row_1idx, grp in aug.groupby(".cohort_row"):
    j = int(row_1idx) - 1   # 0-indexed
    log_surv[j] = grp["log_contrib"].values.sum() if "log_contrib" in grp.columns \
                  else grp.apply(lambda r: np.log1p(-min(r["p_hat"], 1.0 - np.finfo(float).eps)), axis=1).sum()

# Recompute cleanly using grouped log_contrib
log_surv = np.zeros(n_cohort)
aug["log_contrib"] = log_contrib
for row_1idx, grp in aug.groupby(".cohort_row"):
    j = int(row_1idx) - 1
    log_surv[j] = grp["log_contrib"].sum()

pi_j = -np.expm1(log_surv)

# ── ipw weights: 1/pi_j for controls; 1 for cases ────────────────────────────
ipw = np.where(pi_j > 0, 1.0 / pi_j, 0.0)
case_0idx = np.where(cohort["d"].values == 1)[0]
ipw[case_0idx] = 1.0

# ── write results: one row per unique cohort subject in the NCC ───────────────
unique_ncc = (
    ncc[["id", ".cohort_row", "case"]]
    .sort_values(".cohort_row")
    .drop_duplicates(subset=".cohort_row")
    .reset_index(drop=True)
)
# Force case weight = 1 (matches R: ncc_out$ipw_weight[as.logical(ncc_out$case)] <- 1.0)
unique_ncc["ipw_weight"] = [
    1.0 if int(row["case"]) == 1
    else float(ipw[int(row[".cohort_row"]) - 1])
    for _, row in unique_ncc.iterrows()
]

unique_ncc[["id", ".cohort_row", "case", "ipw_weight"]].to_csv(
    "glm_weights_results.csv", index=False
)
print(f"Written glm_weights_results.csv ({len(unique_ncc)} rows)")
print(unique_ncc[unique_ncc["case"] == 0]["ipw_weight"].describe())
