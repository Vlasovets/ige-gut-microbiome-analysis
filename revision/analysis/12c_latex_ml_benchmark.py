#!/usr/bin/env python
"""
Module 12c - supplementary tables for the classification benchmark (R3: random
forest / ensemble benchmark), with the main-paper sparse log-contrast models
(c-lasso) evaluated on the same folds.

Design (module 04): 221 sex/age/BMI-matched pairs for all sensitized (fewer for
subgroups), nested CV with outer StratifiedGroupKFold on pair (5 folds x 10
repeats = 50 test folds), tuning inside training folds, within-pair label
permutation null per model.

Sources per analysis (results/ml/):
  ml_<cell>_metrics.csv          August run (family/genus: L2-only null) and the
                                 species/ASV runs (all models incl. XGBoost)
  ml_<cell>_xgb_metrics.csv      XGBoost backfill, family/genus
  ml_<cell>_perm_metrics.csv     per-model nulls backfill, family/genus (same
                                 seeds/folds; AUCs asserted identical)
  ml_<cell>_classo_metrics.csv   c-lasso C1 log-contrast models (classo-venv)
The Dummy model is not tabulated: chance is judged against each model's own
permutation null (null AUC means 0.48-0.50).

Outputs:
  results/tables/latex/tableS_ml_benchmark.tex      all sensitized, 4 ranks
  results/tables/latex/tableS_ml_benchmark_all.tex  summary over all analyses
  results/tables/ml_benchmark_all_analyses.csv      source data
Usage:  python revision/analysis/12c_latex_ml_benchmark.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

REV = Path(__file__).resolve().parent.parent
ML = REV / "results" / "ml"
TAB = REV / "results" / "tables"
MODELS = [("classo-theoretical", r"Log-contrast, theoretical $\lambda$"),
          ("classo-CV", "Log-contrast, cross-validation"),
          ("classo-StabSel", "Log-contrast, stability selection"),
          ("L1-logistic", "L1-penalized logistic"),
          ("L2-logistic", "L2-penalized logistic"),
          ("RandomForest", "Random forest"),
          ("HistGradBoost", "Gradient boosting"),
          ("XGBoost", "XGBoost")]
# (level, prevalence, only-partition): None = all 15 IgE groups
CELLS = [("family", "prev10", None), ("genus", "prev10", None),
         ("species", "prev05", None), ("asv", "prev05", "any_ige")]
PARTS = ["any_ige", "single", "multiple", "high_score", "medium_score", "low_score"] + \
        [f"cluster_{c}" for c in "ABCDEFGHI"]
RLAB = {"family": "Family (10\\%)", "genus": "Genus (10\\%)",
        "species": "Species (5\\%)", "asv": "ASV (5\\%)"}
COLS = ["model", "auc_pooled", "auc_ci_low", "auc_ci_high", "mcr", "perm_p",
        "null_auc_mean", "n_selected_median"]


MAXDIFF = [0.0]   # largest August-vs-backfill AUC difference accepted


def read(stem):
    f = ML / f"{stem}_metrics.csv"
    return pd.read_csv(f) if f.exists() else None


def cell(level, prev, part):
    stem = f"ml_{level}_{prev}_{part}"
    base, xgb, perm, cl = (read(stem), read(stem + "_xgb"), read(stem + "_perm"),
                           read(stem + "_classo"))
    if base is None:
        return None
    d = base.copy()
    if xgb is not None:
        d = pd.concat([d[d.model != "XGBoost"], xgb[xgb.model == "XGBoost"]])
    if perm is not None:      # per-model nulls; AUCs must reproduce
        # 1e-4, not 1e-9: RandomForest (n_jobs=-1) sums trees in thread order, so
        # one of 26 cells differs at ~1e-5 - float noise, two orders below the
        # 2 decimals reported. Anything larger means a different model: stop.
        m = d.set_index("model").auc_pooled
        for _, r in perm.iterrows():
            diff = abs(m[r.model] - r.auc_pooled)
            assert diff < 1e-4, f"{stem} {r.model}: AUC not reproduced (diff {diff:.2e})"
            MAXDIFF[0] = max(MAXDIFF[0], diff)
        d = d.set_index("model")
        d.loc[perm.model, ["perm_p", "null_auc_mean"]] = perm.set_index("model")[["perm_p", "null_auc_mean"]]
        d = d.reset_index()
    if cl is not None:
        d = pd.concat([d, cl])
    d = d[d.model != "Dummy"].reindex(columns=COLS)
    d["level"], d["prevalence"], d["partition"] = level, prev, part
    d["source_complete"] = cl is not None
    return d


rows = []
for level, prev, only in CELLS:
    for part in ([only] if only else PARTS):
        c = cell(level, prev, part)
        if c is not None:
            rows.append(c)
t = pd.concat(rows, ignore_index=True)
t.to_csv(TAB / "ml_benchmark_all_analyses.csv", index=False)

print(f"August vs backfill AUCs: max |difference| = {MAXDIFF[0]:.2e}")
missing = sorted({(r.level, r.partition) for r in t.itertuples() if not r.source_complete})
if missing:
    print(f"WARNING: c-lasso results missing for {len(missing)} analyses, e.g. {missing[:3]}")
nan_p = t[t.perm_p.isna()]
if len(nan_p):
    print(f"WARNING: {len(nan_p)} model x analysis rows without a permutation p")


def fp(p):
    if pd.isna(p):
        return "n.d."
    return "$<$0.01" if p < 0.01 else f"{p:.2f}"


# ── Table 1: all sensitized, 4 ranks ──────────────────────────────────────
L = [r"\begin{table}[h!]", r"\centering", r"\footnotesize",
     r"\setlength{\tabcolsep}{4pt}",   # 6pt default overflows the iScience text width by ~13pt
     r"\begin{tabular}{l l c c c c c}",
     r"\textbf{Rank} & \textbf{Model} & \textbf{AUC} & \textbf{95\% CI} & \textbf{MCR} & "
     r"\textbf{Permutation $p$} & \textbf{Features selected} \\", r"\hline"]
p = t[t.partition == "any_ige"]
for level, _, _ in CELLS:
    sub = p[p.level == level].set_index("model")
    for i, (m, lab) in enumerate(MODELS):
        if m not in sub.index:
            continue
        r = sub.loc[m]
        nsel = "--" if pd.isna(r.n_selected_median) else f"{r.n_selected_median:.0f}"
        L.append(" & ".join([RLAB[level] if i == 0 else "", lab, f"{r.auc_pooled:.2f}",
                             f"{r.auc_ci_low:.2f}--{r.auc_ci_high:.2f}", f"{r.mcr:.2f}",
                             fp(r.perm_p), nsel]) + r" \\")
    L.append(r"\hline")
L += [r"\end{tabular}",
      r"\caption{Classification of IgE sensitization from gut microbial composition: "
      r"sparse log-contrast models (c-lasso, main analysis) compared with penalized logistic "
      r"regression and ensemble methods. All sensitized vs.\ non-sensitized participants, "
      r"221 sex-, age- and BMI-matched pairs. All models were evaluated on identical "
      r"cross-validation folds (5 outer folds $\times$ 10 repeats, both members of a pair "
      r"always in the same fold), with hyperparameters or the regularization parameter "
      r"selected within the training folds. AUC: pooled out-of-fold area under the ROC curve; "
      r"95\% CI: 2.5th--97.5th percentile across the 50 test folds; MCR: misclassification "
      r"rate; permutation $p$: proportion of within-pair label permutations (log-contrast "
      r"models 100, other models 200) reaching the observed AUC. Features selected: median "
      r"number of taxa with non-zero coefficient per training fold. Log-contrast models: "
      r"c-lasso, squared hinge loss (C1) with intercept and zero-sum constraint on "
      r"log-transformed counts (pseudo-count 1).\label{tab:ml_benchmark}}",
      r"\end{table}"]
(TAB / "latex" / "tableS_ml_benchmark.tex").write_text("\n".join(L) + "\n")

# ── Table 2: summary over every analysis ──────────────────────────────────
s = t.groupby("model").agg(n=("auc_pooled", "size"), med=("auc_pooled", "median"),
                           lo=("auc_pooled", "min"), hi=("auc_pooled", "max"),
                           ci_excl=("auc_ci_low", lambda x: int((x > 0.5).sum())),
                           p_sig=("perm_p", lambda x: int((x < 0.05).sum())),
                           p_n=("perm_p", lambda x: int(x.notna().sum())))
L = [r"\begin{table}[h!]", r"\centering", r"\footnotesize",
     r"\begin{tabular}{l c c c c}",
     r"\textbf{Model} & \textbf{Analyses} & \textbf{AUC, median (range)} & "
     r"\textbf{95\% CI excludes 0.5} & \textbf{Permutation $p<0.05$} \\", r"\hline"]
for m, lab in MODELS:
    if m in s.index:
        r = s.loc[m]
        # .loc on a mixed-dtype row upcasts the counts to float: cast back
        L.append(" & ".join([lab, f"{int(r.n)}", f"{r.med:.2f} ({r.lo:.2f}--{r.hi:.2f})",
                             f"{int(r.ci_excl)}", f"{int(r.p_sig)} of {int(r.p_n)}"]) + r" \\")
L += [r"\hline", r"\end{tabular}",
      r"\caption{Classification benchmark across all analyses: family and genus (10\% "
      r"prevalence filter) and species (5\%) for each of the 15 IgE groups, and ASV (5\%) for "
      r"all sensitized participants, each vs.\ matched non-sensitized controls. Design and "
      r"models as in Table~\ref{tab:ml_benchmark}. The largest AUCs occur in the smallest "
      r"groups and have confidence intervals spanning 0.5.\label{tab:ml_benchmark_all}}",
      r"\end{table}"]
(TAB / "latex" / "tableS_ml_benchmark_all.tex").write_text("\n".join(L) + "\n")

print(p.pivot_table(index="model", columns="level", values="auc_pooled").round(3).to_string())
print(s.to_string())
