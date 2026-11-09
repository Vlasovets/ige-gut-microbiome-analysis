#!/usr/bin/env python
"""
Module 04 - predictive benchmark: can microbial composition classify IgE status?

Answers R3 Major 6 / R2 Minor 9. The published analysis reported only a
misclassification rate, from a repeated 80/20 hold-out, with two sources of
optimism that this module removes:

  LEAK 1  Prevalence filtering was applied before the train/test split.
          Here every transform lives inside a Pipeline and is fitted on training
          folds only. (Note the CLR transform is a WITHIN-SAMPLE operation, so it
          cannot leak across samples; it is still placed in the pipeline for
          clarity. The upstream prevalence filter in 00b_aggregate_otu.py was
          computed on the full 2034-sample table independently of labels, so it is
          not label leakage either - but the in-fold scaler and selector are.)

  LEAK 2  Cross-validation folds could place the two members of a matched pair in
          different folds. Because pairs are matched on sex, age and BMI, that
          leaks the covariate profile of a test sample into training. Fixed with
          StratifiedGroupKFold grouped on pair_nb.

A property worth stating: every matched pair contributes exactly one sensitized
and one non-sensitized individual, so assigning WHOLE PAIRS to folds gives perfect
class balance in every fold automatically. The design does the stratification.

The permutation null permutes labels WITHIN each matched pair (an independent coin
flip per pair). That preserves both the matched structure and the exact class
balance, and is the same exchangeability used for PERMANOVA and the network tests.

Usage:
  python revision/analysis/04_predictive_ml.py \
    --level genus --prevalence prev10 --partition any_ige
  # sparse log-contrast (c-lasso) models only, same folds, writes *_classo_*:
  revision/env/classo-venv/bin/python revision/analysis/04_predictive_ml.py \
    --level genus --prevalence prev10 --partition any_ige --classo
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (average_precision_score, balanced_accuracy_score,
                             confusion_matrix, roc_auc_score, roc_curve)
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, StandardScaler

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

REVISION = Path(__file__).resolve().parent.parent
CONFIG = REVISION / "config" / "config.yaml"


def clr(x: np.ndarray) -> np.ndarray:
    """Centred log-ratio. A within-sample transform: no cross-sample leakage."""
    lx = np.log(x + 1.0)
    return lx - lx.mean(axis=1, keepdims=True)


def build_models(seed: int) -> dict:
    clr_t = FunctionTransformer(clr, validate=False)
    return {
        "Dummy": (Pipeline([("clr", clr_t), ("sc", StandardScaler()),
                            ("m", DummyClassifier(strategy="stratified",
                                                  random_state=seed))]), {}),
        "L1-logistic": (Pipeline([("clr", clr_t), ("sc", StandardScaler()),
                                  ("m", LogisticRegression(penalty="l1",
                                        solver="liblinear", max_iter=5000,
                                        random_state=seed))]),
                        {"m__C": [0.01, 0.1, 1.0]}),
        "L2-logistic": (Pipeline([("clr", clr_t), ("sc", StandardScaler()),
                                  ("m", LogisticRegression(penalty="l2",
                                        max_iter=5000, random_state=seed))]),
                        {"m__C": [0.01, 0.1, 1.0]}),
        "RandomForest": (Pipeline([("clr", clr_t),
                                   ("m", RandomForestClassifier(
                                       n_estimators=1000, min_samples_leaf=3,
                                       max_features="sqrt", n_jobs=-1,
                                       random_state=seed))]),
                         {"m__min_samples_leaf": [1, 3, 5]}),
        # Kept alongside XGBoost: it was the stand-in before xgboost existed
        # in any env, and their agreement is now a plumbing check.
        "HistGradBoost": (Pipeline([("clr", clr_t),
                                    ("m", HistGradientBoostingClassifier(
                                        random_state=seed))]),
                          {"m__max_leaf_nodes": [7, 15]}),
    }


def add_xgboost(models: dict, seed: int) -> dict:
    """Add XGBoost when the library is importable (revision/env/xgb-venv).

    Lazy on purpose: the script must keep running under plain slr-env, where
    xgboost is deliberately not installed (slr-env is frozen)."""
    try:
        from xgboost import XGBClassifier
    except ImportError:
        return models
    clr_t = FunctionTransformer(clr, validate=False)
    # Tiny grid on purpose: n=442 with <=1.5k features does not support a big
    # hyperparameter search, and n_jobs=1 leaves parallelism to the outer job.
    models["XGBoost"] = (Pipeline([("clr", clr_t),
                                   ("m", XGBClassifier(
                                       tree_method="hist", n_estimators=500,
                                       eval_metric="logloss", n_jobs=1,
                                       random_state=seed))]),
                         {"m__max_depth": [2, 3, 4],
                          "m__learning_rate": [0.05, 0.1]})
    return models


def evaluate(X, y, groups, seed_list, n_splits=5,
             models=None) -> tuple[pd.DataFrame, dict]:
    """Nested CV: outer folds grouped on pair, inner GridSearchCV for tuning."""
    rows, roc_store = [], {}
    available = add_xgboost(build_models(seed_list[0]), seed_list[0])
    if models:
        missing = [m for m in models if m not in available]
        if missing:
            raise SystemExit(f"requested model(s) not available here: {missing} "
                             f"(XGBoost needs revision/env/xgb-venv/bin/python)")
        available = {k: v for k, v in available.items() if k in models}
    for name, (pipe, grid) in available.items():
        oof_y, oof_p, fold_aucs = [], [], []
        for rep, seed in enumerate(seed_list):
            outer = StratifiedGroupKFold(n_splits=n_splits, shuffle=True,
                                         random_state=seed)
            for tr, te in outer.split(X, y, groups=groups):
                # Leakage assertion: no pair may appear on both sides.
                assert not (set(groups[tr]) & set(groups[te])), \
                    "pair leaked across the train/test split"
                est = pipe
                if grid:
                    inner = StratifiedGroupKFold(n_splits=4, shuffle=True,
                                                 random_state=seed)
                    est = GridSearchCV(pipe, grid, scoring="roc_auc",
                                       cv=inner.split(X[tr], y[tr], groups[tr]),
                                       n_jobs=1, refit=True)
                    est.fit(X[tr], y[tr])
                else:
                    est.fit(X[tr], y[tr])
                if hasattr(est, "predict_proba"):
                    p = est.predict_proba(X[te])[:, 1]
                else:
                    p = est.decision_function(X[te])
                oof_y.append(y[te]); oof_p.append(p)
                if len(np.unique(y[te])) == 2:
                    fold_aucs.append(roc_auc_score(y[te], p))
        yy = np.concatenate(oof_y); pp = np.concatenate(oof_p)
        fold_aucs = np.array(fold_aucs)
        # Primary CI: percentile across outer folds, which accounts for split
        # variability. A bootstrap over pooled predictions would ignore fold
        # dependence and give a falsely narrow interval.
        lo, hi = (np.percentile(fold_aucs, [2.5, 97.5])
                  if len(fold_aucs) > 3 else (np.nan, np.nan))
        pred = (pp > np.median(pp)).astype(int)
        tn, fp, fn, tp = confusion_matrix(yy, pred).ravel()
        rows.append({
            "model": name, "auc_pooled": roc_auc_score(yy, pp),
            "auc_mean_fold": fold_aucs.mean(), "auc_sd_fold": fold_aucs.std(),
            "auc_ci_low": lo, "auc_ci_high": hi,
            "pr_auc": average_precision_score(yy, pp),
            "balanced_acc": balanced_accuracy_score(yy, pred),
            "mcr": float((pred != yy).mean()),
            "tn": tn, "fp": fp, "fn": fn, "tp": tp,
            "n_folds": len(fold_aucs),
        })
        roc_store[name] = (yy, pp)
        print(f"    {name:14s} AUC={rows[-1]['auc_pooled']:.3f} "
              f"[{lo:.3f}, {hi:.3f}]  MCR={rows[-1]['mcr']:.3f}")
    return pd.DataFrame(rows), roc_store


# ── Sparse log-contrast classification (c-lasso), the main-paper method ──────
# Same settings as notebooks/04_classification.ipynb (C1 = squared hinge, intercept,
# zero-sum constraint, Path-Alg, Nlam 160, CV Nsubset 5, StabSel threshold 0.5,
# theoretical lambda) EXCEPT the feature transform: the notebook passed relative
# abundances (transform='rel'), which under a zero-sum constraint is not a
# log-contrast model. Here X = log(count + 1): with sum(beta) = 0 the per-sample
# log library size cancels, so this equals the log-contrast of relative abundance
# with pseudo-count 1. One c-lasso solve per fold yields all three selections.
CLASSO_MODES = {"classo-theoretical": "LAMfixed", "classo-CV": "CV",
                "classo-StabSel": "StabSel"}


def classo_scores(Xtr, ytr, Xte, seed, nlam=160) -> tuple[dict, dict]:
    """Fit c-lasso C1 on a training fold; return test decision scores per
    selection mode and the selected-feature mask per mode (intercept dropped)."""
    import contextlib, io
    from classo import classo_problem
    Ltr, Lte = np.log(Xtr + 1.0), np.log(Xte + 1.0)
    pb = classo_problem(Ltr, np.where(ytr == 1, 1.0, -1.0), np.ones((1, Ltr.shape[1])))
    f, ms = pb.formulation, pb.model_selection
    f.classification, f.concomitant, f.huber, f.intercept = True, False, False, True
    ms.LAMfixed, ms.CV, ms.StabSel, ms.PATH = True, True, True, True
    ms.CVparameters.seed, ms.CVparameters.Nsubset = seed, 5
    ms.PATHparameters.Nlam = nlam
    ms.PATHparameters.numerical_method = "Path-Alg"
    ms.LAMfixedparameters.lam = "theoretical"
    ms.StabSelparameters.method, ms.StabSelparameters.threshold = "first", 0.5
    ms.StabSelparameters.seed = seed
    with contextlib.redirect_stdout(io.StringIO()):
        pb.solve()
    sol = pb.solution
    betas = {"classo-theoretical": sol.LAMfixed.beta, "classo-CV": sol.CV.beta,
             "classo-StabSel": sol.StabSel.refit}
    scores, selected = {}, {}
    for name, b in betas.items():
        b = np.asarray(b, float).ravel()          # [intercept, beta_1..beta_p]
        sel = np.abs(b[1:]) > 1e-5                # c-lasso's own selection cut-off
        # Score with the SELECTED coefficients only: sub-threshold numerical
        # residue otherwise yields a non-constant score (and a spurious AUC)
        # from a model that selected nothing.
        scores[name] = b[0] + Lte @ np.where(sel, b[1:], 0.0)
        selected[name] = sel
    return scores, selected


def evaluate_classo(X, y, groups, seed_list, n_splits=5, taxa=None):
    """Nested evaluation identical in folds to evaluate(): same seeds, same
    StratifiedGroupKFold on pair_nb. Lambda is chosen inside each training fold
    by c-lasso itself (theoretical / internal CV / stability selection)."""
    store = {m: ([], [], []) for m in CLASSO_MODES}     # oof_y, oof_p, fold_aucs
    sel_count = {m: np.zeros(X.shape[1]) for m in CLASSO_MODES}
    n_sel = {m: [] for m in CLASSO_MODES}
    n_fits = 0
    for seed in seed_list:
        outer = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        for tr, te in outer.split(X, y, groups=groups):
            assert not (set(groups[tr]) & set(groups[te])), "pair leaked across the split"
            sc, sel = classo_scores(X[tr], y[tr], X[te], seed)
            n_fits += 1
            for m in CLASSO_MODES:
                oy, op, fa = store[m]
                oy.append(y[te]); op.append(sc[m])
                if len(np.unique(y[te])) == 2:
                    fa.append(roc_auc_score(y[te], sc[m]))
                sel_count[m] += sel[m]; n_sel[m].append(int(sel[m].sum()))
    rows, roc_store = [], {}
    for m, (oy, op, fa) in store.items():
        yy, pp, fa = np.concatenate(oy), np.concatenate(op), np.array(fa)
        lo, hi = np.percentile(fa, [2.5, 97.5]) if len(fa) > 3 else (np.nan, np.nan)
        pred = (pp > np.median(pp)).astype(int)
        tn, fp, fn, tp = confusion_matrix(yy, pred).ravel()
        rows.append({
            "model": m, "auc_pooled": roc_auc_score(yy, pp),
            "auc_mean_fold": fa.mean(), "auc_sd_fold": fa.std(),
            "auc_ci_low": lo, "auc_ci_high": hi,
            "pr_auc": average_precision_score(yy, pp),
            "balanced_acc": balanced_accuracy_score(yy, pred),
            "mcr": float((pred != yy).mean()),
            # the main paper's rule: predict sign(score); ties (beta = 0) -> +1
            "mcr_sign": float(((pp >= 0).astype(int) != yy).mean()),
            "n_selected_median": float(np.median(n_sel[m])),
            "tn": tn, "fp": fp, "fn": fn, "tp": tp, "n_folds": len(fa)})
        roc_store[m] = (yy, pp)
        print(f"    {m:18s} AUC={rows[-1]['auc_pooled']:.3f} [{lo:.3f}, {hi:.3f}]  "
              f"MCR={rows[-1]['mcr']:.3f}  sign-MCR={rows[-1]['mcr_sign']:.3f}  "
              f"median selected={rows[-1]['n_selected_median']:.0f}")
    freq = pd.DataFrame({m: sel_count[m] / n_fits for m in CLASSO_MODES}, index=taxa)
    freq.index.name = "taxon"
    return pd.DataFrame(rows), roc_store, freq


def permutation_null_classo(X, y, groups, seed, n_perm, n_splits) -> pd.DataFrame:
    """Within-pair label permutations (as permutation_null); one c-lasso solve per
    fold serves all three selection modes."""
    rng = np.random.default_rng(seed)
    uniq = np.unique(groups)
    rows = []
    for _ in range(n_perm):
        yp = y.copy()
        for g in uniq:
            idx = np.where(groups == g)[0]
            if len(idx) == 2 and rng.random() < 0.5:
                yp[idx] = yp[idx][::-1]
        oy, op = [], {m: [] for m in CLASSO_MODES}
        outer = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        try:
            for tr, te in outer.split(X, yp, groups=groups):
                sc, _ = classo_scores(X[tr], yp[tr], X[te], seed)
                oy.append(yp[te])
                for m in CLASSO_MODES:
                    op[m].append(sc[m])
        except Exception:
            continue
        for m in CLASSO_MODES:
            rows.append({"model": m, "null_auc": roc_auc_score(np.concatenate(oy),
                                                               np.concatenate(op[m]))})
    return pd.DataFrame(rows)


def permutation_null(X, y, groups, seed, n_perm=200, n_splits=5,
                     model="L2-logistic") -> np.ndarray:
    """Null AUC from permuting labels WITHIN each matched pair.

    Computed per model: a null from one model does not license a p-value for
    another, and the letter quotes model-specific AUCs."""
    rng = np.random.default_rng(seed)
    pipe = add_xgboost(build_models(seed), seed)[model][0]
    uniq = np.unique(groups)
    aucs = []
    for _ in range(n_perm):
        yp = y.copy()
        # One independent coin flip per pair: swap that pair's two labels or not.
        for g in uniq:
            idx = np.where(groups == g)[0]
            if len(idx) == 2 and rng.random() < 0.5:
                yp[idx] = yp[idx][::-1]
        oy, op = [], []
        outer = StratifiedGroupKFold(n_splits=n_splits, shuffle=True,
                                     random_state=seed)
        try:
            for tr, te in outer.split(X, yp, groups=groups):
                pipe.fit(X[tr], yp[tr])
                op.append(pipe.predict_proba(X[te])[:, 1]); oy.append(yp[te])
            aucs.append(roc_auc_score(np.concatenate(oy), np.concatenate(op)))
        except Exception:
            continue
    return np.array(aucs)


def detectable_auc(n_pos: int, n_neg: int, power=0.80, alpha=0.05) -> float:
    """Smallest AUC detectable at the given power (Hanley-McNeil variance)."""
    from scipy.stats import norm
    za, zb = norm.ppf(1 - alpha / 2), norm.ppf(power)
    for auc in np.arange(0.50, 0.999, 0.001):
        q1 = auc / (2 - auc); q2 = 2 * auc ** 2 / (1 + auc)
        v = (auc * (1 - auc) + (n_pos - 1) * (q1 - auc ** 2)
             + (n_neg - 1) * (q2 - auc ** 2)) / (n_pos * n_neg)
        if (auc - 0.5) / np.sqrt(v) >= za + zb:
            return float(auc)
    return float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", default="genus")
    ap.add_argument("--prevalence", default="prev10")
    ap.add_argument("--partition", default="any_ige")
    ap.add_argument("--repeats", type=int, default=10)
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--models", default=None,
                    help="comma list; default = every available model")
    ap.add_argument("--classo", action="store_true",
                    help="run only the c-lasso log-contrast models "
                         "(needs revision/env/classo-venv); writes *_classo_*")
    ap.add_argument("--n-perm-classo", type=int, default=100,
                    help="within-pair permutations for the c-lasso null")
    ap.add_argument("--tag", default="",
                    help="output-stem suffix so backfill runs never clobber "
                         "the CSVs behind already-quoted numbers")
    args = ap.parse_args()

    cfg = yaml.safe_load(CONFIG.read_text())
    pipe_root = Path(cfg["paths"]["pipeline_results"])
    seed = cfg["statistics"]["seed"]
    out = REVISION / "results" / "ml"; out.mkdir(parents=True, exist_ok=True)
    figs = REVISION / "results" / "figures"; figs.mkdir(parents=True, exist_ok=True)

    src = pipe_root / args.level / args.prevalence / args.partition / "source"
    oi = pd.read_csv(src / "otu_ige.csv", index_col=0)
    oc = pd.read_csv(src / "otu_non_ige.csv", index_col=0)
    mi = pd.read_csv(src / "sample_data_ige.csv", index_col=0)
    mo = pd.read_csv(src / "sample_data_non_ige.csv", index_col=0)

    taxa = [c for c in oi.columns if c in set(oc.columns)]
    Xdf = pd.concat([oi[taxa], oc[taxa]])
    y = np.r_[np.ones(len(oi), int), np.zeros(len(oc), int)]
    groups = np.r_[mi["pair_nb"].values, mo["pair_nb"].values]

    # Keep only intact pairs, so within-pair permutation is well defined.
    vc = pd.Series(groups).value_counts()
    keep = np.isin(groups, vc[vc == 2].index)
    Xdf, y, groups = Xdf[keep], y[keep], groups[keep]
    X = Xdf.values.astype(float)

    print(f"=== ML benchmark: {args.level} / {args.prevalence} / {args.partition} ===")
    print(f"  n={len(y)} ({y.sum()} sensitized, {(1-y).sum()} not), "
          f"features={X.shape[1]}, complete pairs={len(np.unique(groups))}")
    if len(np.unique(groups)) < 10:
        print("  too few pairs; skipping")
        return 0

    n_splits = min(5, len(np.unique(groups)) // 2)
    seeds = [seed + i for i in range(args.repeats)]
    model_list = ([m.strip() for m in args.models.split(",")]
                  if args.models else None)
    print(f"\n  Nested CV: {n_splits} outer folds x {args.repeats} repeats, "
          f"grouped on pair_nb")

    if args.classo:
        # Sparse log-contrast models only (revision/env/classo-venv), same folds.
        import classo
        print(f"  c-lasso {classo.__version__} from {classo.__file__}")
        res, _, freq = evaluate_classo(X, y, groups, seeds, n_splits, taxa=taxa)
        null = permutation_null_classo(X, y, groups, seed, args.n_perm_classo, n_splits)
        stats = null.groupby("model").null_auc.agg(["mean", "std"])
        res["null_auc_mean"] = res.model.map(stats["mean"])
        res["null_auc_sd"] = res.model.map(stats["std"])
        res["perm_p"] = [float((null.null_auc[null.model == m] >= a).mean())
                         for m, a in zip(res.model, res.auc_pooled)]
        res["n_perm"] = null.groupby("model").size().reindex(res.model).values
        for _, r in res.iterrows():
            print(f"    {r.model:18s} null mean={r.null_auc_mean:.4f} (must be ~0.500) "
                  f"empirical p={r.perm_p:.3f}")
        res["level"], res["prevalence"] = args.level, args.prevalence
        res["partition"] = args.partition
        res["classo_version"] = classo.__version__
        res["detectable_auc_80pct"] = detectable_auc(int(y.sum()), int((1 - y).sum()))
        stem = f"ml_{args.level}_{args.prevalence}_{args.partition}_classo"
        res.to_csv(out / f"{stem}_metrics.csv", index=False)
        freq.to_csv(out / f"{stem}_selection.csv")
        null.to_csv(out / f"{stem}_null.csv", index=False)
        print(f"\nWrote {out / f'{stem}_metrics.csv'}")
        return 0

    res, roc_store = evaluate(X, y, groups, seeds, n_splits, models=model_list)

    # Per-model permutation nulls. Every model whose AUC may be quoted gets its
    # own within-pair null; the Dummy is excluded (its null IS the null).
    null_models = [m for m in res.model if m != "Dummy"]
    print(f"\n  Permutation null ({args.n_perm} within-pair label permutations, "
          f"per model: {', '.join(null_models)})")
    null_rows, null_stats = [], {}
    for m in null_models:
        nm = permutation_null(X, y, groups, seed, args.n_perm, n_splits, model=m)
        obs_m = float(res.loc[res.model == m, "auc_pooled"].iloc[0])
        pm = float((nm >= obs_m).mean()) if len(nm) else np.nan
        null_stats[m] = (nm.mean() if len(nm) else np.nan,
                         nm.std() if len(nm) else np.nan, pm)
        null_rows += [{"model": m, "null_auc": v} for v in nm]
        print(f"    {m:14s} null mean={null_stats[m][0]:.4f} "
              f"sd={null_stats[m][1]:.4f} (must be ~0.500)  "
              f"observed={obs_m:.3f}  empirical p={pm:.3f}")
    null = np.array([r["null_auc"] for r in null_rows
                     if r["model"] == "L2-logistic"])
    if not len(null) and null_rows:
        null = np.array([r["null_auc"] for r in null_rows])
    obs = (float(res.loc[res.model == "L2-logistic", "auc_pooled"].iloc[0])
           if (res.model == "L2-logistic").any()
           else float(res.auc_pooled.iloc[0]))
    p_emp = float((null >= obs).mean()) if len(null) else np.nan

    da = detectable_auc(int(y.sum()), int((1 - y).sum()))
    print(f"\n  Power: with {int(y.sum())} vs {int((1-y).sum())}, the smallest AUC "
          f"detectable at 80% power is {da:.3f}")
    print(f"  -> the null is INFORMATIVE: an effect of AUC >= {da:.2f} would have "
          f"been detected.")

    # Feature importance from the honest source: permutation on held-out folds.
    print("\n  Permutation importance (RandomForest, held-out folds)")
    imp_rows = []
    outer = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    rf = build_models(seed)["RandomForest"][0]
    for tr, te in outer.split(X, y, groups=groups):
        rf.fit(X[tr], y[tr])
        r = permutation_importance(rf, X[te], y[te], n_repeats=5,
                                   random_state=seed, scoring="roc_auc")
        imp_rows.append(r.importances_mean)
    imp = pd.DataFrame({"taxon": taxa,
                        "perm_importance_mean": np.mean(imp_rows, axis=0),
                        "perm_importance_sd": np.std(imp_rows, axis=0)})
    imp = imp.sort_values("perm_importance_mean", ascending=False)
    print(imp.head(8).to_string(index=False))

    # ── Figures ────────────────────────────────────────────────────────
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    for name, (yy, pp) in roc_store.items():
        fpr, tpr, _ = roc_curve(yy, pp)
        ax[0].plot(fpr, tpr, lw=1.4,
                   label=f"{name} (AUC={roc_auc_score(yy, pp):.3f})")
    ax[0].plot([0, 1], [0, 1], "k--", lw=0.8)
    ax[0].set_xlabel("False positive rate", fontsize=10)
    ax[0].set_ylabel("True positive rate", fontsize=10)
    ax[0].set_title(f"{args.partition} — {args.level}, {args.prevalence}", fontsize=10)
    ax[0].legend(fontsize=8, loc="lower right")

    if len(null):
        ax[1].hist(null, bins=30, color="#95A5A6", edgecolor="white")
        ax[1].axvline(obs, color="#C0392B", lw=2,
                      label=f"observed {obs:.3f}\nempirical p={p_emp:.3f}")
        ax[1].axvline(0.5, color="black", ls="--", lw=0.8, label="chance")
        ax[1].set_xlabel("AUC under within-pair label permutation", fontsize=10)
        ax[1].set_ylabel("count", fontsize=10)
        ax[1].set_title("Permutation null", fontsize=10)
        ax[1].legend(fontsize=8)
    fig.tight_layout()
    stem = f"ml_{args.level}_{args.prevalence}_{args.partition}"
    if args.tag:
        stem += f"_{args.tag}"
    for ext in ("pdf", "png"):
        fig.savefig(figs / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    res["level"] = args.level; res["prevalence"] = args.prevalence
    res["partition"] = args.partition
    res["null_auc_mean"] = res.model.map(lambda m: null_stats.get(m, (np.nan,)*3)[0])
    res["null_auc_sd"]   = res.model.map(lambda m: null_stats.get(m, (np.nan,)*3)[1])
    res["perm_p"]        = res.model.map(lambda m: null_stats.get(m, (np.nan,)*3)[2])
    res["perm_p_L2"] = p_emp        # kept for compatibility with earlier runs
    res["models_requested"] = ",".join(model_list) if model_list else "all"
    res["detectable_auc_80pct"] = da
    res.to_csv(out / f"{stem}_metrics.csv", index=False)
    imp.to_csv(out / f"{stem}_importance.csv", index=False)
    if null_rows:
        pd.DataFrame(null_rows).to_csv(out / f"{stem}_null.csv", index=False)
    print(f"\nWrote {out / f'{stem}_metrics.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
