#!/usr/bin/env python
"""
Module 02g - consensus DA under the two confounder-control designs, side by side.

  matched  1:1 pairs (sex exact, age +-5, BMI +-4) + covariates   (modules 02/02d)
  full     all sensitized vs all 233 controls + covariates         (02f + --design full)

Only the five covariate-adjusted methods are compared (LinDA, ANCOM-BC2, ALDEx2,
ZicoSeq, limma-voom). The paired Wilcoxon exists only in the matched design, so
it is left out of both counts to keep them comparable.

Consensus is scored exactly as in 02b/02e: q_within < 0.05 (BH within one
analysis), and a taxon counts as consensus when >= 2 methods call it with the
same sign as every other method that returned an estimate.

Outputs (results/tables/):
  da_design_method_rates.csv    % significant per method x level x design
  da_design_consensus.csv       every taxon significant in >= 1 method, either design
  da_design_focus_taxa.csv      the families/taxa the submitted paper highlighted,
                                plus the P. copri aggregate, any_ige, all prevalences

Usage:  python revision/analysis/02g_full_vs_matched.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

REVISION = Path(__file__).resolve().parent.parent
RES = REVISION / "results"
TABLES = RES / "tables"
ALPHA = 0.05
METHODS = ["LinDA", "ANCOMBC2", "ALDEx2", "ZicoSeq", "limma_voom"]
KEY = ["level", "prevalence", "partition", "taxon"]
FOCUS = {"family": ["Bacteroidaceae", "Veillonellaceae", "Lachnospiraceae",
                    "Oscillospiraceae", "Prevotellaceae"],
         "genus": ["Prevotella"],
         "species": ["Prevotella_copri"]}
COLS = KEY + ["method", "effect_log2", "p_raw", "q_within", "n_ige", "n_ctrl"]


def load(design):
    da, ext = ("da_full", "da_extended_full") if design == "full" else ("da", "da_extended")
    files = sorted((RES / da).glob("da_long_*.csv")) + sorted((RES / ext).glob("daext_long_*.csv"))
    d = pd.concat([pd.read_csv(f, usecols=COLS) for f in files], ignore_index=True)
    d = d[d.method.isin(METHODS)].copy()
    d["design"] = design
    return d


d = pd.concat([load("matched"), load("full")], ignore_index=True)
cells = d.groupby(["design", "method"])[["level", "prevalence", "partition"]] \
         .apply(lambda x: len(x.drop_duplicates())).unstack("method")
print("Analyses per design x method (expect 180 each):\n", cells, sep="")
assert (cells.values == 180).all(), "incomplete grid - do not interpret"

d["sig"] = d.q_within < ALPHA
rates = (d.groupby(["design", "level", "method"])
          .agg(n_tests=("sig", "size"), n_sig=("sig", "sum")).reset_index())
rates["pct_sig"] = (100 * rates.n_sig / rates.n_tests).round(2)
rates.to_csv(TABLES / "da_design_method_rates.csv", index=False)


def score(g):
    eff = g.effect_log2.dropna()
    sig = g[g.sig]
    concordant = bool(len(eff)) and (np.sign(eff).nunique() == 1)
    return pd.Series({
        "n_methods": len(eff), "n_sig": len(sig),
        "sig_methods": ";".join(sorted(sig.method)),
        "sign_concordant": concordant,
        "direction": "up" if eff.median() > 0 else "down",
        "consensus": len(sig) >= 2 and concordant,
    })


hit = d.groupby(["design"] + KEY).sig.transform("any")
s = d[hit].groupby(["design"] + KEY).apply(score, include_groups=False).reset_index()
wide = d[hit].pivot_table(index=["design"] + KEY, columns="method",
                          values=["effect_log2", "q_within"])
wide.columns = [f"{m}_{v.replace('effect_log2', 'effect')}" for v, m in wide.columns]
s = s.merge(wide.reset_index(), on=["design"] + KEY, how="left")
s.sort_values(["consensus", "n_sig"], ascending=False).to_csv(
    TABLES / "da_design_consensus.csv", index=False)

foc = d[(d.partition == "any_ige") &
        d.apply(lambda r: r.taxon in FOCUS.get(r.level, []), axis=1)]
foc = foc.pivot_table(index=["taxon", "level", "prevalence", "method"],
                      columns="design", values=["effect_log2", "q_within"])
foc.columns = [f"{v}_{des}" for v, des in foc.columns]
foc = foc.reset_index()
foc.to_csv(TABLES / "da_design_focus_taxa.csv", index=False)

# ── console summary ──────────────────────────────────────────────────────
n = d.drop_duplicates(["design", "partition"]).set_index(["design", "partition"])
print("\nany_ige sample sizes:",
      {k: (int(v.n_ige), int(v.n_ctrl)) for k, v in
       n.xs("any_ige", level="partition").iterrows()})
print("\n% significant, all ranks pooled:")
print(rates.groupby(["design", "method"])[["n_sig", "n_tests"]].sum()
      .assign(pct=lambda x: (100 * x.n_sig / x.n_tests).round(2)).pct.unstack("design"))
for des in ["matched", "full"]:
    c = s[(s.design == des) & s.consensus]
    print(f"\nConsensus (>=2 of 5 methods, concordant sign), {des}: {len(c)} taxon-analyses")
    if len(c):
        print(c[KEY + ["n_sig", "sig_methods", "direction"]].to_string(index=False))
print("\nFocus taxa, any_ige (min q across 5 methods; direction = median effect sign):")
g = d[(d.partition == "any_ige") & d.apply(lambda r: r.taxon in FOCUS.get(r.level, []), axis=1)]
print(g.groupby(["taxon", "prevalence", "design"])
       .agg(min_q=("q_within", "min"), n_sig=("sig", "sum"),
            n_up=("effect_log2", lambda x: int((x > 0).sum())))
       .unstack("design").round(3).to_string())
