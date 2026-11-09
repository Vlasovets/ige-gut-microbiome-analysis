#!/usr/bin/env python
"""
Module 12d - supplementary tables of the log-contrast classification per IgE
group (replace the submitted out-of-sample MCR tables tab:out_sample_mcr and
tab:out_sample_mcr_clusters, which came from relative abundances, sex-only
matching and seven ranks).

Same layout as the submitted tables: one block per lambda selection scheme, one
row per rank, one column per group. Cell: pooled out-of-fold AUC (median number
of taxa selected per training fold); bold = within-pair permutation p < 0.05.
Source: results/ml/ml_<rank>_<prev>_<group>_classo_metrics.csv (module 04
--classo). Column headers give the number of matched pairs.

Outputs: results/tables/latex/tableS_classo_auc_burden.tex   (tab:classo_auc_burden)
         results/tables/latex/tableS_classo_auc_clusters.tex (tab:classo_auc_clusters)
Usage:  python revision/analysis/12d_latex_classo_groups.py
"""

from pathlib import Path

import pandas as pd

REV = Path(__file__).resolve().parent.parent
ML = REV / "results" / "ml"
OUT = REV / "results" / "tables" / "latex"
SCHEMES = [("classo-theoretical", r"Theoretical $\lambda_{\text{fixed}}$"),
           ("classo-CV", r"Cross-validation $\lambda_{\text{CV}}$"),
           ("classo-StabSel", r"Stability selection $\lambda_{\text{SS}}$")]
RANKS = [("family", "prev10", "Family"), ("genus", "prev10", "Genus"),
         ("species", "prev05", "Species"), ("asv", "prev05", "ASV")]
BURDEN = [("any_ige", "Any IgE"), ("single", "Single"), ("multiple", "Multiple"),
          ("high_score", "High score"), ("medium_score", "Medium score"),
          ("low_score", "Low score")]
CLUSTERS = [(f"cluster_{c}", c) for c in "ABCDEFGHI"]


def load(rank, prev, part):
    f = ML / f"ml_{rank}_{prev}_{part}_classo_metrics.csv"
    return pd.read_csv(f).set_index("model") if f.exists() else None


def table(groups, label, caption):
    pairs = {}
    for part, _ in groups:
        m = load("family", "prev10", part)
        pairs[part] = int(round((m.tn.iloc[0] + m.fp.iloc[0]) / 10))   # 10 repeats
    ncol = len(groups)
    head = " & ".join(rf"\textbf{{{lab}}} ({pairs[p]})" for p, lab in groups)
    L = [r"\begin{sidewaystable}[h!]", r"\small", r"\centering",
         r"\begin{tabular}{l l " + "c " * ncol + "}",
         r" & \textbf{Rank} & " + head + r" \\", r"\hline"]
    for model, slab in SCHEMES:
        L.append(rf" & \textbf{{{slab}}} " + "& " * ncol + r"\\")
        for rank, prev, rlab in RANKS:
            cells, any_value = [], False
            for part, _ in groups:
                m = load(rank, prev, part)
                if m is None:
                    cells.append("--"); continue
                r = m.loc[model]; any_value = True
                v = f"{r.auc_pooled:.2f}"
                if r.perm_p < 0.05:
                    v = rf"\textbf{{{v}}}"
                cells.append(f"{v} ({r.n_selected_median:.0f})")
            if any_value:
                L.append(f"& {rlab} & " + " & ".join(cells) + r" \\")
        L.append(r"\hline")
    L += [r"\end{tabular}", rf"\caption{{{caption}\label{{{label}}}}}", r"\end{sidewaystable}"]
    return "\n".join(L) + "\n"


COMMON = (r"Out-of-sample AUC of the log-contrast classification models, with the median "
          r"number of taxa selected per training fold in parentheses. Each group was compared "
          r"with its sex-, age- and BMI-matched non-sensitized controls (number of matched pairs "
          r"in parentheses in the column header); models were evaluated by cross-validation "
          r"(5 folds $\times$ 10 repeats, both members of a pair in the same fold). Bold: "
          r"permutation $p<0.05$ (100 within-pair label permutations). Prevalence filter: 10\% "
          r"(family, genus), 5\% (species, ASV)")
(OUT / "tableS_classo_auc_burden.tex").write_text(table(
    BURDEN, "tab:classo_auc_burden",
    r"Log-contrast classification of the IgE frequency and score groups across taxonomic "
    r"ranks. " + COMMON + r"; the ASV rank was analyzed for all sensitized individuals only."))
(OUT / "tableS_classo_auc_clusters.tex").write_text(table(
    CLUSTERS, "tab:classo_auc_clusters",
    r"Log-contrast classification of the IgE clusters (A--I) across taxonomic ranks. "
    + COMMON + r". Clusters A and F have 12 and 10 pairs, too few for reliable classification."))
print("wrote", OUT / "tableS_classo_auc_burden.tex", "and", OUT / "tableS_classo_auc_clusters.tex")
