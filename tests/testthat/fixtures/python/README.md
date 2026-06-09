# Python cross-language oracles

Cross-language oracle fixtures for matchatr's classical estimators. Each
estimator is cross-checked against an independent **Python** implementation so a
bug shared between matchatr and its R engine (e.g. a `survival::clogit` wrapper
quirk) cannot hide.

## How it works

Each oracle is a triple of committed files in this directory:

- `<slug>_data.csv` — the shared dataset. **Both** the Python script and the R
  test read this exact file, so the two fit the *same* data and an exact
  numerical comparison is possible (no cross-language RNG mismatch). Generated
  once from the package's data-generating processes (`tests/testthat/helper-dgp.R`).
- `<slug>.py` — reads `<slug>_data.csv`, fits the equivalent estimator with
  `statsmodels`, and writes `<slug>_results.csv` (point estimate, SE, and Wald
  CI on the reported scale).
- `<slug>_results.csv` — the committed Python output. The R test reads it and
  asserts agreement; **tests never invoke Python at run time**, so CI needs no
  Python toolchain. Each comparison is guarded with
  `skip_if(!file.exists(test_path(...)))`.

## Oracle library

`statsmodels` anchors the classical maximum-likelihood estimators, which it
implements directly:

| matchatr estimator | estimand | statsmodels oracle |
|---|---|---|
| `logistic` (unmatched CC) | conditional OR | `smf.logit` (`Logit`) |
| `mh` | summary OR | `StratifiedTable` (Mantel–Haenszel) |
| `clogit` (matched CC) | conditional OR | `ConditionalLogit` |
| `clogit` (nested CC) | hazard ratio | `ConditionalLogit` (risk-set partial likelihood) |
| `polytomous` | subtype OR | `MNLogit` |
| `test_homogeneity` | Wald χ² + pooled OR | hand-computed from `MNLogit` params / cov |
| `compute_ncc_weights(method="glm")` | per-subject ipw_weight | `statsmodels.Logit` (same augmented dataset + product formula) |

The GLM-weights oracle (`glm_weights.*`) differs from the others: it uses a
*pair* of data files (`glm_weights_cohort.csv` + `glm_weights_ncc.csv`) and
the Python script builds the augmented selection dataset from scratch, mirroring
matchatr's `build_ncc_selection_dataset()`.  Both apply the identical product
formula on the same logistic fit, so agreement is within double-precision
rounding (< 1e-10 typically).

`delicatessen` (M-estimation + sandwich) is reserved for the **causal /
sandwich** estimands (the case-control-weighted g-formula / IPW / AIPW / TMLE
contrasts) when those phases land — that is where estimating-equation oracles
are the natural fit, not the conditional/multinomial MLEs here.

## Environment

- Python ≥ 3.11 (developed on 3.14)
- `numpy`, `scipy`, `pandas`, `statsmodels`

## Regenerating a results fixture

```sh
cd tests/testthat/fixtures/python
python3 <slug>.py        # rewrites <slug>_results.csv from <slug>_data.csv
```

Regenerate only when the estimator's reported quantity legitimately changes;
commit the updated `_results.csv` alongside the code change. The `_data.csv`
files are static inputs — do not regenerate them, or the oracle drifts.
