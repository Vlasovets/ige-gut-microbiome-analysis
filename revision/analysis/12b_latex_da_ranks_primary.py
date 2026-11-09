#!/usr/bin/env python
"""
Module 12b - LaTeX table tab:da_ranks_primary: differential abundance of the
primary contrast across ranks, prevalence filters and methods.

Primary contrast: all 275 sensitized vs the 233 IgE-negative controls, full
cohort, ~ IgE status + sex + age + BMI (02f + --design full). One row per
analysis (rank x prevalence filter). Prevalence "None" = the submitted
manuscript's unfiltered feature set (every feature present in >= 1 of the 508
participants; 02f --unfiltered + slurm/02h), shown when its results exist.

Counting rules (unchanged from the hand-built first version of this table,
which this script reproduces exactly for the filtered rows):
  Features   features in the analysis = the denominator for every method;
             features a method did not test count as not significant
  method %   features with q_within < 0.05 (BH within the analysis)
  Hits       features significant in >= 1 method
  >=2        hits significant in >= 2 methods (% of hits)
Methods that call nothing in any row are dropped and named in the caption.

Output: results/tables/latex/tableS_da_ranks_primary.tex
Usage:  python revision/analysis/12b_latex_da_ranks_primary.py
"""

from pathlib import Path

import pandas as pd

REV = Path(__file__).resolve().parent.parent
RES = REV / "results"
OUT = RES / "tables" / "latex" / "tableS_da_ranks_primary.tex"
ALPHA, PART = 0.05, "any_ige"
METHODS = ["ANCOMBC2", "limma_voom", "ZicoSeq", "LinDA", "ALDEx2"]
HEAD = {"ANCOMBC2": r"\textbf{ANCOM-BC2}",
        "limma_voom": r"\textbf{limma-voom}",
        "ZicoSeq": r"\textbf{ZicoSeq}", "LinDA": r"\textbf{LinDA}", "ALDEx2": r"\textbf{ALDEx2}"}
NAME = {"ANCOMBC2": "ANCOM-BC2", "limma_voom": "limma-voom", "ZicoSeq": "ZicoSeq",
        "LinDA": "LinDA", "ALDEx2": "ALDEx2"}
LEVELS = [("family", "Family"), ("genus", "Genus"), ("species", "Species"), ("asv", "ASV")]
PREVS = [("prev00", "None"), ("prev01", r"1\%"), ("prev05", r"5\%"), ("prev10", r"10\%")]
WORDS = {2: "two", 3: "three", 4: "four", 5: "five"}


def load(level, prev):
    f1 = RES / "da_full" / f"da_long_{level}_{prev}_{PART}.csv"
    f2 = RES / "da_extended_full" / f"daext_long_{level}_{prev}_{PART}.csv"
    if not (f1.exists() and f2.exists()):
        return None
    d = pd.concat([pd.read_csv(f1), pd.read_csv(f2)], ignore_index=True)
    return d[d.method.isin(METHODS)]


rows = []
for level, _ in LEVELS:
    for prev, _ in PREVS:
        d = load(level, prev)
        if d is None:
            if prev != "prev00":
                raise FileNotFoundError(f"{level}/{prev} missing")
            continue
        n_feat = d[d.method == "LinDA"].taxon.nunique()
        d = d.assign(sig=d.q_within < ALPHA)
        per = d.groupby("method").agg(n_sig=("sig", "sum"), n_tested=("taxon", "nunique"))
        nsig_taxon = d.groupby("taxon").sig.sum()
        r = {"level": level, "prev": prev, "features": n_feat,
             "hits": int((nsig_taxon >= 1).sum()), "ge2": int((nsig_taxon >= 2).sum()),
             "max_methods": int(nsig_taxon.max())}
        for m in METHODS:
            r[f"{m}_sig"] = int(per.n_sig.get(m, 0))
            r[f"{m}_tested"] = int(per.n_tested.get(m, 0))
        rows.append(r)
t = pd.DataFrame(rows)
t.to_csv(RES / "tables" / "da_ranks_primary.csv", index=False)

shown = [m for m in METHODS if t[f"{m}_sig"].sum() > 0]
silent = [m for m in METHODS if m not in shown]
has_none = (t.prev == "prev00").any()
plab = dict(PREVS)

# ── LaTeX ─────────────────────────────────────────────────────────────────
ncol = 3 + len(shown) + 2
L = [r"\begin{table}[h!]", r"\centering", r"\footnotesize", r"\setlength{\tabcolsep}{2.5pt}",
     r"\begin{tabular}{l r r " + "r " * len(shown) + "r r}",
     r" & & & \multicolumn{" + str(len(shown)) + r"}{c}{\textbf{Features significant (\%)}} & & \\",
     r"\textbf{Rank} & \textbf{Prevalence} & \textbf{Features} & "
     + " & ".join(HEAD[m] for m in shown)
     + r" & \textbf{Hits} & \textbf{$\geq$2 methods} \\",
     r"\hline"]
for level, llab in LEVELS:
    sub = t[t.level == level]
    for i, (_, r) in enumerate(sub.iterrows()):
        cells = [llab if i == 0 else "", plab[r.prev], f"{r.features:,}"]
        cells += [f"{100 * r[f'{m}_sig'] / r.features:.1f}" for m in shown]
        ge2 = f"{r.ge2} ({100 * r.ge2 / r.hits:.1f}\\%)" if r.hits else "0 (0.0\\%)"
        cells += [f"{r.hits:,}", ge2]
        L.append(" & ".join(cells) + r" \\")
    L.append(r"\hline")
L.append(r"\end{tabular}")

# caption pieces derived from the data, never hand-typed numbers
fewer = {}
for m in METHODS:
    ps = [plab[p] for p in t.prev.unique()
          if (t[t.prev == p][f"{m}_tested"] < t[t.prev == p].features).any()]
    if ps:
        fewer[m] = ps
fewer_txt = ""
if fewer:
    ms = [NAME[m] for m in METHODS if m in fewer]
    who = ms[0] if len(ms) == 1 else ", ".join(ms[:-1]) + " and " + ms[-1]
    fewer_txt = (f"{who} apply their own internal filtering and test fewer features in some "
                 f"analyses; features a method did not test are counted as not significant. ")
print("tested fewer than all features:", fewer)
silent_txt = ""
if silent:
    names = " and ".join(NAME[m] for m in silent)
    silent_txt = (f"{names} called no feature significant in any of these analyses and "
                  f"{'is' if len(silent) == 1 else 'are'} not shown. ")
# ZicoSeq permutation floor: min attainable BH q = m / (nperm + 1), nperm = 4999
NPERM = 4999
floor_rows = t[t.ZicoSeq_tested / (NPERM + 1) > ALPHA]
floor_txt = ""
if len(floor_rows):
    parts = []
    for lv, sub in floor_rows.groupby("level", sort=False):
        if len(sub) == (t.level == lv).sum():
            parts.append(f"all {dict(LEVELS)[lv]} rows")
        else:
            parts += [f"{dict(LEVELS)[lv]} {('no filter' if p == 'prev00' else plab[p] + ' filter')}"
                      for p in sub.prev]
    where = ", ".join(parts)
    floor_txt = (f"ZicoSeq uses {NPERM:,} permutations, so with more than "
                 f"{int(ALPHA * (NPERM + 1))} tested features it cannot reach $q < 0.05$ "
                 f"({where}). ")
silent_txt += floor_txt
mx = int(t.max_methods.max())
max_txt = (f"none by three or more" if mx <= 2 else
           f"at most {WORDS.get(mx, mx)} methods agreed on any feature")
none_txt = (r" Prevalence ``None'': no prevalence filter, every feature present in at least "
            r"one of the 508 participants, as in the submitted analysis (6,958 ASVs)."
            if has_none else "")
cap = (r"\caption{Differential abundance of IgE sensitization across taxonomic ranks, prevalence "
       r"filters and methods for the primary contrast: all 275 sensitized individuals vs.\ the 233 "
       r"IgE-negative controls (full cohort), adjusted for sex, age and BMI. Each row is one analysis."
       + none_txt +
       r" Features: number of features in the analysis, used as the denominator for every method. "
       + fewer_txt.replace("%", r"\%").replace(r"\\%", r"\%") + silent_txt +
       r"A test is significant at Benjamini--Hochberg $q < 0.05$ within its analysis. Hits: features "
       r"called significant by at least one method; $\geq$2 methods: hits also called significant by "
       r"a second method (" + max_txt + r").\label{tab:da_ranks_primary}}")
L += [cap, r"\end{table}"]
OUT.write_text("\n".join(L) + "\n")
print("wrote", OUT)
print(t[["level", "prev", "features", "hits", "ge2", "max_methods"]
        + [f"{m}_sig" for m in METHODS]].to_string(index=False))
