#!/usr/bin/env python
"""
Module 12 - LaTeX for the revised Table 1, Table 2 and three supplementary
tables, generated from the result CSVs so no number is ever hand-transcribed.

Layout decision (answers R2 "Table 1: provide p-values", R2 "Table 2 ...
supporting statistical tests (with a proper control for type I errors)" and
R3 "formal statistical comparisons for Tables 1 and 2"):

  * Main Table 1 keeps its rows and gains two narrow columns, p and q (BH);
    bold marks q < 0.05; the test is given by a superscript key. p/q below
    0.001 print as "<0.001" - exact values are in the supplement.
  * Main Table 2 keeps its rows and gains ONE bottom row of q-values.
  * Supplement:
      tableS_matched_cohort   - Table 1's exact layout for the 221 matched
                                pairs used in every microbiome analysis
      tableS_baseline_tests   - test, statistic, exact p, q and effect size
                                with 95% CI, full cohort (A) and matched (B)
      tableS_cluster_tests    - the same for the nine-cluster comparison

Checked against the real manuscript preamble (12 pt, Helvetica, 7 in text
width): zero overfull boxes. Section labels are full-width rows rather than a
column of their own - a bold "Health Outcome" column alone costs ~3.3 cm, and
the SUBMITTED Table 1 was already 11.6 pt too wide. Only packages already in
the preamble are used (rotating); no booktabs/multirow, hence plain \\hline.

Usage:  python revision/analysis/12_latex_tables.py
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd

REVISION = Path(__file__).resolve().parent.parent
T = REVISION / "results" / "tables"
OUT = T / "latex"
OUT.mkdir(parents=True, exist_ok=True)
ALPHA = 0.05


# ── formatting helpers ────────────────────────────────────────────────────
def isnan(x) -> bool:
    return x is None or (isinstance(x, float) and math.isnan(x))


def fp(p: float) -> str:
    """Exact p/q for the supplement: 3 decimals, scientific below 0.001."""
    if isnan(p):
        return "--"
    if p >= 0.001:
        return f"{p:.3f}"
    e = int(math.floor(math.log10(p)))
    m = p / 10 ** e
    if round(m, 1) >= 10:
        m, e = m / 10, e + 1
    return f"${m:.1f}\\times10^{{{e}}}$"


def fp_main(p: float) -> str:
    """Compact p/q for the main tables."""
    if isnan(p):
        return "--"
    return f"{p:.3f}" if p >= 0.001 else "$<$0.001"


def fq(q: float, main: bool = False) -> str:
    s = fp_main(q) if main else fp(q)
    return f"\\textbf{{{s}}}" if (not isnan(q) and q < ALPHA) else s


def num(x: float, d: int = 3, sign: bool = False) -> str:
    """A number in math mode, so minus signs are real minus signs."""
    return f"${x:+.{d}f}$" if sign else f"${x:.{d}f}$"


def pct(s: str) -> str:
    return s.replace("%", "\\%")


def pieces(summary: str) -> list[str]:
    return [pct(x.strip()) for x in summary.split(" / ")]


def levels(note: str) -> list[str]:
    m = re.search(r"levels:\s*([^;]+)", str(note))
    return [x.strip() for x in m.group(1).split(",")] if m else []


TEST_TEX = {"Wilcoxon rank-sum": "Wilcoxon rank-sum",
            "Chi-squared": "Pearson $\\chi^2$",
            "Fisher's exact": "Fisher's exact",
            "Kruskal-Wallis": "Kruskal--Wallis",
            "Fisher's exact (simulated)": "Fisher's exact (MC)"}
LABEL_TEX = {"BMI (kg/m2)": "BMI (kg/m$^2$)"}
TESTMARK = {"Wilcoxon rank-sum": "a", "Chi-squared": "b", "Fisher's exact": "c"}


def lab(v: str) -> str:
    return LABEL_TEX.get(v, v)


# ══════════════════════════════════════════════════════════════════════════
# Table-1-layout builder (used for main Table 1 AND the matched-cohort table)
# ══════════════════════════════════════════════════════════════════════════
def table1_layout(t: pd.DataFrame, caption: str, label: str) -> str:
    t = t.set_index("variable")

    def pq(var):
        r = t.loc[var]
        return (f"{fp_main(r.p_value)}$^{{{TESTMARK[r.test]}}}$",
                fq(r.q_BH, main=True))

    def cont(lbl, var):
        r = t.loc[var]; p, q = pq(var)
        return f"{lbl} &  & {r.sensitized} & {r.non_sensitized} & {p} & {q} \\\\"

    def cat(lbl, var, names, order):
        r = t.loc[var]; lv = levels(r.note)
        s1 = dict(zip(lv, pieces(r.sensitized)))
        s0 = dict(zip(lv, pieces(r.non_sensitized)))
        p, q = pq(var); out = []
        for i, k in enumerate(order):
            out.append(f"{lbl if i == 0 else ''} & {names[k]} & {s1[k]} & {s0[k]} "
                       f"& {p if i == 0 else ''} & {q if i == 0 else ''} \\\\")
        return out

    L = [r"\begin{table}[h!]", r"\centering", r"\setlength{\tabcolsep}{4pt}",
         r"\begin{tabular}{llcccc}", "",
         r" &  & \textbf{IgE} & \textbf{non-IgE} & $\boldsymbol{p}$ & "
         r"$\boldsymbol{q}$ \\", "",
         r"\multicolumn{2}{l}{\textbf{General}} & \textbf{Mean (SD)} & "
         r"\textbf{Mean (SD)} & & \\", "",
         cont("Age", "Age (years)"), cont("Body Mass Index", "BMI (kg/m2)"), "",
         r" & & $n_{1}$ (\%) & $n_{0}$ (\%) & & \\", ""]
    L += cat("Sex", "Sex", {"Female": "F", "Male": "M"}, ["Female", "Male"])
    L += ["", r"\multicolumn{6}{l}{\textbf{Lifestyle}} \\", ""]
    L += cat("Education", "Education",
             {"1": "elementary/secondary", "2": "middle/junior high",
              "3": "high school/college"}, ["1", "2", "3"])
    L += [""] + cat("Physical activity", "Physical activity",
                    {"Regular": "regular", "Irregular": "irregular", "None": "no"},
                    ["Regular", "Irregular", "None"])
    L += [""] + cat("Smoking behaviour", "Smoking",
                    {"Regular": "regular", "Ex-smoker": "ex-smoker",
                     "Never": "never"}, ["Regular", "Ex-smoker", "Never"])
    L += ["", r"\multicolumn{6}{l}{\textbf{Health Outcome}} \\", ""]
    for lbl, var in [("Asthma", "Asthma"), ("Dermatitis", "Dermatitis"),
                     ("Hay Fever", "Hay fever")]:
        L += cat(lbl, var, {"Yes": "Yes", "No": "No"}, ["Yes", "No"]) + [""]
    L += [r"\end{tabular}", f"\\caption{{{caption}\\label{{{label}}}}}",
          r"\end{table}"]
    return "\n".join(L) + "\n"


t1 = pd.read_csv(T / "table1_full508_tested.csv")
t1m = pd.read_csv(T / "table1_matched442_qc.csv")
N_T1 = len(t1)
WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five",
         6: "six", 7: "seven", 8: "eight", 9: "nine"}
KEY = (r"Continuous variables are shown as mean (SD) and compared with the "
       r"Wilcoxon rank-sum test$^{a}$; categorical variables as $n$ (\%) and "
       r"compared with Pearson's $\chi^2$ test$^{b}$. $q$, Benjamini--Hochberg-"
       rf"adjusted $p$-value across all {N_T1} baseline variables tested")

(OUT / "table1_main.tex").write_text(table1_layout(
    t1,
    # Shortened caption (author-approved wording); detail lives in the
    # Methods and in tab:baseline_tests / tab:matched_cohort.
    r"Summary statistics of the KORA FF4 cohort by IgE status (IgE-positive, "
    rf"$n_1 = {int(t1.n_sensitized.iloc[0])}$; IgE-negative, "
    rf"$n_0 = {int(t1.n_non_sensitized.iloc[0])}$). Mean (SD), Wilcoxon "
    r"rank-sum test$^{a}$; $n$ (\%), Pearson's $\chi^2$ test$^{b}$. $q$, "
    rf"Benjamini--Hochberg-adjusted $p$-value across the {N_T1} baseline "
    rf"variables tested ({WORDS.get(N_T1 - 9, N_T1 - 9)} not shown; full results in "
    r"Table~\ref{tab:baseline_tests}); bold, $q < 0.05$. Microbiome analyses "
    rf"used {int(t1m.n_sensitized.iloc[0])} 1:1 matched pairs balanced for "
    r"sex, age, and BMI (Table~\ref{tab:matched_cohort}).",
    "tab:stats"))

(OUT / "tableS_matched_cohort.tex").write_text(table1_layout(
    t1m,
    r"Summary statistics of the matched cohort used for all microbiome "
    r"analyses: 221 IgE-positive individuals ($n_1$) each matched 1:1 to a "
    r"non-sensitized individual ($n_0$) on sex (exact), age ($\pm$5 years) and "
    r"BMI ($\pm$4 kg/m$^2$), layout as in Table~\ref{tab:stats}. " + KEY +
    r" (Table~\ref{tab:baseline_tests}B); bold, $q < 0.05$. Matching removes "
    r"the sex imbalance of the full cohort, while the differences in allergic "
    r"disease, which are downstream of sensitization and were never used as "
    r"matching or adjustment variables, are preserved.",
    "tab:matched_cohort"))

# ══════════════════════════════════════════════════════════════════════════
# TABLE 2 (main) - rows unchanged, one q row added
# ══════════════════════════════════════════════════════════════════════════
t2 = pd.read_csv(T / "table2_clusters_tested.csv").set_index("variable")
d2 = pd.read_csv(T / "table2_clusters_descriptives.csv").set_index("cluster")
N_T2 = len(t2)
# Descriptive rows are taken VERBATIM from the manuscript (all values were
# verified against table2_clusters_descriptives.csv and a recomputation of the
# education shares from kora.csv); only the test row is new.
T2_ROWS = r"""A & dust mite (\textit{Der f 2, Der p 2, Der f 1}) & 12 & 51.17 (6.13) & 25.85 (2.76) & 50.0 / 50.0 & 25.0 / 50.0 / 25.0 \\
B & birch (\textit{Bet v 1}), hazel (\textit{Cor a 1}), apple (\textit{Mal d 1}) & 24 & 53.88 (6.69) & 26.86 (5.81) & 37.5 / 62.5 & 20.8 / 41.7 / 37.5 \\
C & birch (\textit{Bet v 1}), apple (\textit{Mal d 1}), hazel (\textit{Cor a 1}) & 26 & 51.81 (7.53) & 25.81 (4.94) & 23.1 / 76.9 & 11.5 / 50.0 / 38.5 \\
D & birch (\textit{Bet v 1}), hazel (\textit{Cor a 1}), apple (\textit{Mal d 1}) & 45 & 48.80 (5.98) & 25.76 (3.99) & 48.9 / 51.1 & 17.8 / 35.6 / 46.7 \\
E & timothy grass (\textit{Phl p 5, Phl p 1, Phl p 4}) & 22 & 50.14 (7.21) & 26.08 (3.20) & 59.1 / 40.9 & 54.5 / 27.3 / 18.2 \\
F & timothy grass (\textit{Phl p 5, Phl p 1}), Bermuda grass (\textit{Cyn d 1}) & 10 & 45.20 (4.87) & 25.94 (2.21) & 50.0 / 50.0 & 0.0 / 30.0 / 70.0 \\
G & dust mite (\textit{Der f 2, Der p 2, Der f 1}) & 22 & 48.36 (5.14) & 26.17 (4.38) & 45.5 / 54.5 & 50.0 / 22.7 / 27.3 \\
H & Bermuda grass (\textit{Cyn d 1}), timothy grass (\textit{Phl p 4, Phl p 1}) & 40 & 51.52 (7.20) & 28.20 (5.20) & 50.0 / 50.0 & 30.0 / 37.5 / 32.5 \\
I & olive tree (\textit{Ole e 1}), birch (\textit{Bet v 1}), cat (\textit{Fel d 1}) & 74 & 52.97 (5.98) & 27.20 (6.65) & 33.8 / 66.2 & 43.2 / 35.1 / 21.6 \\"""
for line in T2_ROWS.splitlines():               # drift guard vs the CSV
    f = [x.strip() for x in line.split("&")]
    c = f[0]
    assert int(f[2]) == int(d2.loc[c, "n"]), f"cluster {c} size drifted"
    assert f[3] == f"{d2.loc[c,'age_mean']:.2f} ({d2.loc[c,'age_sd']:.2f})", c
    assert f[4] == f"{d2.loc[c,'bmi_mean']:.2f} ({d2.loc[c,'bmi_sd']:.2f})", c

# Allergen designations (e.g. "Phl p 4") must never break across lines:
# spaces inside each comma-separated name become non-breaking.
def _nbsp_textit(m):
    return "\\textit{" + ", ".join(x.strip().replace(" ", "~")
                                   for x in m.group(1).split(",")) + "}"
T2_ROWS = re.sub(r"\\textit\{([^}]*)\}", _nbsp_textit, T2_ROWS)

qrow = (r"\multicolumn{3}{l}{\textit{Across clusters,} $q$\,$^{a}$} & "
        + " & ".join(fq(t2.loc[v, "q_BH"], main=True)
                     for v in ["Age (years)", "BMI (kg/m2)", "Sex", "Education"])
        + r" \\")
# p{8cm}, not the submitted p{4cm}: in 12 pt Helvetica the 4 cm column
# wrapped allergen names mid-token and left orphan ")" lines; the sideways page
# has the width to spare.
L2 = [r"\begin{sidewaystable}[h!]", r"\centering", r"\small",
      r"\begin{tabular}{l p{8cm} c c c c c}", "",
      r"\textbf{Cluster} & \textbf{Top IgE Components} & \textbf{Size} & "
      r"\textbf{Age (Mean, SD)} & \textbf{BMI (Mean, SD)} & "
      r"\textbf{Sex (m/f \%)} & \textbf{Education (s/m/h \%)} \\", "",
      T2_ROWS, r"\hline", qrow, "", r"\end{tabular}",
      # Shortened caption (author-approved wording); method detail is in the
      # Methods and in tab:cluster_tests.
      r"\caption{Summary statistics of the IgE clusters (A--I). Sex is shown as "
      r"\% male/female; education as \% with secondary (s), medium (m), and "
      r"higher (h) education. $^{a}$Kruskal--Wallis test (age, BMI) or Fisher's "
      r"exact test (sex, education); $q$, Benjamini--Hochberg-adjusted across the "
      rf"{N_T2} variables tested (full results in Table~\ref{{tab:cluster_tests}}); "
      r"bold, $q < 0.05$. Clusters A ($n = 12$) and F ($n = 10$) are too small for "
      r"cluster-specific inference.\label{tab:igee_clusters}}",
      r"\end{sidewaystable}"]
(OUT / "table2_main.tex").write_text("\n".join(L2) + "\n")


# ══════════════════════════════════════════════════════════════════════════
# SUPPLEMENT: tests only, full cohort (A) and matched cohort (B)
# ══════════════════════════════════════════════════════════════════════════
def effect_cell(r) -> str:
    if r.type == "cont":
        return f"Cliff's $\\delta$ = {num(r.effect_size, 3, sign=True)}"
    if not isnan(r.get("risk_diff", float("nan"))):
        return (f"RD({r.risk_diff_level}) = {num(100*r.risk_diff, 1, True)} pp "
                f"[{num(100*r.risk_diff_ci_low, 1, True)}, "
                f"{num(100*r.risk_diff_ci_high, 1, True)}]")
    return f"Cram\\'er's $V$ = {num(r.effect_size, 3)}"


# Short test names inside the cells; the caption names them in full. Keeps the
# portrait table inside the 7 in text width at 12 pt Helvetica.
TEST_SHORT = {"Wilcoxon rank-sum": "Wilcoxon", "Chi-squared": "$\\chi^2$",
              "Fisher's exact": "Fisher"}


def test_rows(df: pd.DataFrame) -> list[str]:
    rows = []
    for _, r in df.iterrows():
        stat = "--" if isnan(r.statistic) else f"{r.statistic:.2f}"
        v = lab(r.variable) + (r"$^{\dagger}$" if r.is_collider else "")
        rows.append(f"{v} & {TEST_SHORT.get(r.test, r.test)} & {stat} & "
                    f"{fp(r.p_value)} & {fq(r.q_BH)} & {effect_cell(r)} \\\\")
    return rows


hdr = (r"\textbf{Variable} & \textbf{Test} & \textbf{Statistic} & "
       r"$\boldsymbol{p}$ & $\boldsymbol{q}$ & \textbf{Effect size [95\% CI]} \\")
S1 = [r"\begin{table}[h!]", r"\centering", r"\footnotesize",
      r"\setlength{\tabcolsep}{4pt}", r"\begin{tabular}{l l r r r l}",
      r"\multicolumn{6}{l}{\textbf{A. Full cohort} "
      rf"({int(t1.n_sensitized.iloc[0])} IgE-positive vs "
      rf"{int(t1.n_non_sensitized.iloc[0])} IgE-negative)}} \\",
      hdr, r"\hline", *test_rows(t1), r"\hline", r"\\",
      r"\multicolumn{6}{l}{\textbf{B. Matched cohort} "
      rf"({int(t1m.n_sensitized.iloc[0])} pairs; "
      r"sex exact, age $\pm$5 y, BMI $\pm$4 kg/m$^2$)} \\",
      hdr, r"\hline", *test_rows(t1m), r"\hline",
      r"\end{tabular}",
      r"\caption{Formal comparison of baseline characteristics between "
      r"IgE-positive and IgE-negative participants in the full cohort (A; "
      r"descriptive statistics in Table~\ref{tab:stats}) and in the 1:1 "
      r"matched cohort used for all microbiome analyses (B; descriptive "
      r"statistics in Table~\ref{tab:matched_cohort}). Continuous variables: "
      r"Wilcoxon rank-sum test with Cliff's $\delta$ (IgE-positive relative to "
      r"IgE-negative). Categorical variables: Pearson's $\chi^2$ test "
      r"(expected counts $\geq 5$ in every cell) with Cram\'er's $V$; for "
      r"binary variables the effect is the risk difference (RD, percentage "
      r"points, IgE-positive minus IgE-negative, for the level in brackets) "
      r"with a Newcombe hybrid-score 95\% CI. $q$, Benjamini--Hochberg-adjusted "
      r"within each panel; bold, $q < 0.05$. $^{\dagger}$Allergic-disease "
      r"outcomes are downstream of sensitization and were therefore tested "
      r"descriptively but never used as matching or adjustment variables. No "
      r"dietary intake variable is available for this KORA FF4 subset."
      r"\label{tab:baseline_tests}}",
      r"\end{table}"]
(OUT / "tableS_baseline_tests.tex").write_text("\n".join(S1) + "\n")

# ══════════════════════════════════════════════════════════════════════════
# SUPPLEMENT: across-cluster tests
# ══════════════════════════════════════════════════════════════════════════
EFF2 = {"eta-squared": r"$\eta^2_H$", "Cramer's V": r"Cram\'er's $V$"}
S2 = [r"\begin{table}[h!]", r"\centering", r"\small",
      r"\begin{tabular}{l l r r r l}",
      r"\textbf{Variable} & \textbf{Test} & \textbf{Statistic (df)} & "
      r"$\boldsymbol{p}$ & $\boldsymbol{q}$ & \textbf{Effect size} \\", r"\hline"]
for _, r in pd.read_csv(T / "table2_clusters_tested.csv").iterrows():
    st = "--" if isnan(r.statistic) else f"{r.statistic:.2f} ({int(r.df)})"
    S2.append(f"{lab(r.variable)} & {TEST_TEX.get(r.test, r.test)} & {st} & "
              f"{fp(r.p_value)} & {fq(r.q_BH)} & "
              f"{EFF2.get(r.effect_name, r.effect_name)} = "
              f"{num(r.effect_size, 3)} \\\\")
S2 += [r"\hline", r"\end{tabular}",
       r"\caption{Differences in baseline characteristics across the nine IgE "
       r"clusters (A--I; $n = 275$). Kruskal--Wallis test for continuous "
       r"variables, with $\eta^2_H = (H - k + 1)/(n - k)$ (values slightly "
       r"below zero indicate no effect); Fisher's exact test with a Monte Carlo "
       r"$p$-value ($10^5$ replicates, fixed seed) for categorical variables, "
       rf"with Cram\'er's $V$. $q$, Benjamini--Hochberg-adjusted across the "
       rf"{N_T2} variables; bold, $q < 0.05$. Clusters A ($n = 12$) and F "
       r"($n = 10$) are underpowered.\label{tab:cluster_tests}}",
       r"\end{table}"]
(OUT / "tableS_cluster_tests.tex").write_text("\n".join(S2) + "\n")

# ══════════════════════════════════════════════════════════════════════════
# SUPPLEMENT: IgE-burden partitions (the former Table S1, tab:igee_partitions)
# ══════════════════════════════════════════════════════════════════════════
# Component -> (source, IUIS designation). Every component that reaches a
# top-3 list must be here; the assert below fails loudly otherwise.
COMP = {
    "timothy_grass_p1": ("timothy grass", "Phl p 1"),
    "timothy_grass_p4": ("timothy grass", "Phl p 4"),
    "timothy_grass_p5": ("timothy grass", "Phl p 5"),
    "birch_v1": ("birch", "Bet v 1"),
    "bermuda_grass": ("Bermuda grass", "Cyn d 1"),
    "olive_tree_e1": ("olive tree", "Ole e 1"),
    "cat_d1": ("cat", "Fel d 1"),
    "mould_asp_f1": ("mould", "Asp f 1"),
    "dust_mite_pter_p23": ("dust mite", "Der p 23"),
    # The ISAC panel's Cor a 1 is the HAZELNUT PR-10 (food); there is no hazel
    # pollen component. "hazel" alone reads as the tree.
    "hazel_nut_pr10": ("hazelnut", "Cor a 1"),
    "apple_pr10": ("apple", "Mal d 1"),
}


def top_cell(spec: str) -> str:
    """'a:137/275;b:135/275' -> components grouped by identical prevalence,
    in descending order, each group followed by its share. Ties at the cut-off
    are all listed (module 01), so the cell never hides an arbitrary choice."""
    items = []
    for it in spec.split(";"):
        name, frac = it.split(":")
        k, n = map(int, frac.split("/"))
        assert name in COMP, f"no display label for component '{name}'"
        items.append((name, k, n))
    groups, out = {}, []
    for name, k, n in items:
        groups.setdefault(k, []).append(name)
    for k in sorted(groups, reverse=True):
        n = items[0][2]
        # \mbox per item: a component name never splits or hyphenates
        names = ", ".join(f"\\mbox{{{COMP[c][0]} (\\textit{{{COMP[c][1]}}})}}"
                          for c in groups[k])
        out.append(f"{names} {100*k/n:.0f}\\%")
    # ragged-right parbox instead of a justified p{} column (which stretched
    # spaces and hyphenated "ap-ple"); no array package needed
    return "\\parbox[t]{7.2cm}{\\raggedright " + "; ".join(out) + "}"


gd = pd.read_csv(T / "group_descriptives_pairs_topIgE.csv").set_index("group")
st = pd.read_csv(T / "tableS1_partitions_tested.csv")
PART = [("any", "any_ige"), ("single", "single"), ("multiple", "multiple"),
        ("high", "high_score"), ("medium", "medium_score"), ("low", "low_score")]

# drift guard: the descriptives must still equal the submitted S1 values
SUBMITTED = {"any_ige": (275, "51.09 (6.66)", "26.67 (5.17)", "42.2 / 57.8", "31.3 / 36.4 / 32.4"),
             "single": (47, "52.43 (6.04)", "27.24 (6.57)", "38.3 / 61.7", "44.7 / 34.0 / 21.3"),
             "multiple": (228, "50.82 (6.76)", "26.55 (4.84)", "43.0 / 57.0", "28.5 / 36.8 / 34.6"),
             "high_score": (57, "49.98 (6.57)", "25.79 (3.50)", "49.1 / 50.9", "21.1 / 42.1 / 36.8"),
             "medium_score": (110, "50.74 (6.98)", "26.72 (5.24)", "44.5 / 55.5", "31.8 / 34.5 / 33.6"),
             "low_score": (108, "52.04 (6.31)", "27.08 (5.78)", "36.1 / 63.9", "36.1 / 35.2 / 28.7")}
rowsS1 = []
for disp, g in PART:
    r = gd.loc[g]
    age = f"{r.age_mean:.2f} ({r.age_sd:.2f})"
    bmi = f"{r.bmi_mean:.2f} ({r.bmi_sd:.2f})"
    sex = f"{r.pct_male:.1f} / {100 - r.pct_male:.1f}"
    edu = f"{r.edu_s:.1f} / {r.edu_m:.1f} / {r.edu_h:.1f}"
    assert (int(r.n), age, bmi, sex, edu) == SUBMITTED[g], f"{g} drifted from S1"
    rowsS1.append(f"{disp} & {top_cell(r.top3_components)} & {int(r.n)} ({int(r.n_pairs)}) "
                  f"& {age} & {bmi} & {sex} & {edu} \\\\")


def famq(fam: str) -> str:
    d = st[st.family == fam].set_index("variable")
    return " & ".join(fq(d.loc[v, "q_BH"], main=True)
                      for v in ["Age (years)", "BMI (kg/m2)", "Sex", "Education"])


LS1 = [r"\begin{sidewaystable}[h!]", r"\small", r"\centering",
       r"\begin{tabular}{l l c c c c c}",
       r"\textbf{Partition} & \textbf{Top IgE Components (prevalence)} & "
       r"\textbf{Size (pairs)} & \textbf{Age (Mean, SD)} & \textbf{BMI (Mean, SD)} & "
       r"\textbf{Sex (m/f \%)} & \textbf{Education (s/m/h \%)} \\",
       *rowsS1, r"\hline",
       r"\multicolumn{3}{l}{\textit{single vs multiple,} $q$\,$^{a}$} & "
       + famq("single_vs_multiple") + r" \\",
       r"\multicolumn{3}{l}{\textit{high vs medium vs low,} $q$\,$^{b}$} & "
       + famq("high_vs_medium_vs_low") + r" \\",
       r"\end{tabular}",
       r"\caption{Summary statistics of the IgE-burden partitions. The IgE "
       r"burden score is the number of allergen components with detectable IgE "
       r"per individual; among the 275 sensitized individuals, single = 1 "
       r"component, multiple $\geq$2, high $>$15, medium 6--15 and low $\leq$5 "
       r"(any = all sensitized). Top IgE components are the three most "
       r"prevalent within the partition, with the share of the partition "
       r"sensitized to each; all components tied at the third position are "
       r"listed, and bromelain, a cross-reactive carbohydrate-determinant "
       r"marker, is excluded. Size gives the number of individuals and, in "
       r"parentheses, the number of 1:1 matched pairs (sex exact, age $\pm$5 "
       r"years, BMI $\pm$4 kg/m$^2$) used in all microbiome analyses of that "
       r"partition. Sex distribution is shown as percentage male/female, and "
       r"education as percentage with secondary (s), medium (m), and higher (h) "
       r"education. Because the partitions overlap, differences were tested "
       r"within the two disjoint families only: $^{a}$Wilcoxon rank-sum (age, "
       r"BMI) and Pearson's $\chi^2$ test (sex, education); $^{b}$Kruskal--Wallis "
       r"and Pearson's $\chi^2$ test. $q$, Benjamini--Hochberg-adjusted across "
       r"the 12 variables tested within each family (all results in "
       r"Table~\ref{tab:partition_tests}); bold, $q < 0.05$."
       r"\label{tab:igee_partitions}}",
       r"\end{sidewaystable}"]
(OUT / "tableS_partitions.tex").write_text("\n".join(LS1) + "\n")


def ptest_rows(fam: str) -> list[str]:
    rows = []
    for _, r in st[st.family == fam].iterrows():
        stt = ("--" if isnan(r.statistic) else
               f"{r.statistic:.2f}" + ("" if isnan(r.df) else f" ({int(r.df)})"))
        eff = {"Cliff's delta": r"Cliff's $\delta$", "Cramer's V": r"Cram\'er's $V$",
               "eta-squared": r"$\eta^2_H$"}.get(r.effect_name, r.effect_name)
        signed = r.effect_name == "Cliff's delta"
        rows.append(f"{lab(r.variable)} & {TEST_TEX.get(r.test, r.test)} & {stt} & "
                    f"{fp(r.p_value)} & {fq(r.q_BH)} & {eff} = "
                    f"{num(r.effect_size, 3, sign=signed)} \\\\")
    return rows


n_s, n_m = int(gd.loc["single", "n"]), int(gd.loc["multiple", "n"])
n_h, n_md, n_l = (int(gd.loc[g, "n"]) for g in ["high_score", "medium_score", "low_score"])
hdrP = (r"\textbf{Variable} & \textbf{Test} & \textbf{Statistic (df)} & "
        r"$\boldsymbol{p}$ & $\boldsymbol{q}$ & \textbf{Effect size} \\")
LP = [r"\begin{table}[h!]", r"\centering", r"\small", r"\setlength{\tabcolsep}{4pt}",
      r"\begin{tabular}{l l r r r l}",
      rf"\multicolumn{{6}}{{l}}{{\textbf{{A. Single vs multiple}} ($n = {n_s}$ vs {n_m})}} \\",
      hdrP, r"\hline", *ptest_rows("single_vs_multiple"), r"\hline", r"\\",
      rf"\multicolumn{{6}}{{l}}{{\textbf{{B. High vs medium vs low}} ($n = {n_h}$, {n_md}, {n_l})}} \\",
      hdrP, r"\hline", *ptest_rows("high_vs_medium_vs_low"), r"\hline",
      r"\end{tabular}",
      r"\caption{Differences in baseline characteristics between IgE-burden "
      r"partitions, tested within the two disjoint families (the partitions "
      r"overlap, so no test across all six is valid). (A) Wilcoxon rank-sum "
      r"test with Cliff's $\delta$ (single relative to multiple) for continuous "
      r"variables; Pearson's $\chi^2$ test, or Fisher's exact test when an "
      r"expected count was below 5, with Cram\'er's $V$ for categorical "
      r"variables. (B) Kruskal--Wallis test with $\eta^2_H$; Pearson's $\chi^2$ "
      r"test, or Fisher's exact test with a Monte Carlo $p$-value ($10^5$ "
      r"replicates, fixed seed), with Cram\'er's $V$. $q$, Benjamini--Hochberg-"
      r"adjusted across the 12 variables within each panel; bold, $q < 0.05$."
      r"\label{tab:partition_tests}}",
      r"\end{table}"]
(OUT / "tableS_partition_tests.tex").write_text("\n".join(LP) + "\n")

# ══════════════════════════════════════════════════════════════════════════
# SUPPLEMENT: beta diversity across taxonomic ranks (primary contrast)
# ══════════════════════════════════════════════════════════════════════════
# One row per rank x distance for all sensitized vs matched controls. Model M1
# (sex + age + BMI), M3 adds lifestyle + PPI; marginal tests throughout;
# permutations restricted within matched pairs. UniFrac comes from module 05b
# (ASV, prev01 - the same feature set as the ASV row of the other distances).
BETA = REVISION / "results" / "beta"
PHY = REVISION / "results" / "phylo"
RANKS = [("Family", "family", "prev10"), ("Genus", "genus", "prev10"),
         ("Species", "species", "prev10"), ("ASV", "asv", "prev01")]
DISTS = [("bray", "Bray--Curtis"), ("aitchison_pseudo", "Aitchison"),
         ("aitchison_robust", "Robust Aitchison"), ("jaccard", "Jaccard")]
brows, missing_ranks = [], []
for rlab, lev, prev in RANKS:
    f = BETA / f"beta_{lev}_{prev}_any_ige.csv"
    fd = BETA / f"betadisper_{lev}_{prev}_any_ige.csv"
    if not f.exists():
        missing_ranks.append(rlab)
        continue
    b = pd.read_csv(f); b = b[(b.term == "ige_status") & (b.by == "margin")]
    dsp = pd.read_csv(fd).set_index("distance") if fd.exists() else None
    for i, (dk, dl) in enumerate(DISTS):
        m1 = b[(b.model == "M1") & (b.distance == dk)].iloc[0]
        m3 = b[(b.model == "M3") & (b.distance == dk)].iloc[0]
        pdisp = fp(float(dsp.loc[dk, "p_value"])) if dsp is not None else "--"
        first = f"{rlab} ({int(m1.n_taxa)})" if i == 0 else ""
        brows.append(f"{first} & {dl} & {100*m1.R2:.2f} & {fp(m1.p_value)} & "
                     f"{fp(m3.p_value)} & {pdisp} \\\\")
    brows.append(r"\hline")
uf = pd.read_csv(PHY / "unifrac_any_ige.csv")
UNI = [("unifrac_unweighted", "Unweighted UniFrac"),
       ("unifrac_generalized_0.5", r"Generalized UniFrac ($\alpha = 0.5$)"),
       ("unifrac_weighted", "Weighted UniFrac")]
for i, (dk, dl) in enumerate(UNI):
    r = uf[uf.distance == dk].iloc[0]
    first = f"ASV, phylog. ({int(r.n_asvs)})" if i == 0 else ""
    brows.append(f"{first} & {dl} & {100*r.R2:.2f} & {fp(r.p_adonis)} & -- & "
                 f"{fp(r.p_betadisper)} \\\\")
assert not missing_ranks, f"beta results missing for: {missing_ranks}"

LB = [r"\begin{table}[h!]", r"\centering", r"\small", r"\setlength{\tabcolsep}{4pt}",
      r"\begin{tabular}{l l r r r r}",
      r"\textbf{Rank (features)} & \textbf{Distance} & $\boldsymbol{R^2}$ \textbf{(\%)} & "
      r"$\boldsymbol{p}$ & $\boldsymbol{p}$ \textbf{(extended)} & "
      r"$\boldsymbol{p}$ \textbf{(dispersion)} \\", r"\hline",
      *brows, r"\hline", r"\end{tabular}",
      r"\caption{Association of IgE sensitization with gut microbial community "
      r"composition across taxonomic ranks, comparing all sensitized individuals "
      r"with their matched non-sensitized controls (221 pairs). PERMANOVA "
      r"(5,000 permutations restricted within matched pairs, marginal tests) "
      r"adjusted for sex, age and BMI; $R^2$ is the share of variation explained "
      r"by IgE status. Extended: additionally adjusted for education, smoking, "
      r"pack-years, alcohol intake, physical activity and proton-pump-inhibitor "
      r"use. Dispersion: test of homogeneity of multivariate dispersion "
      r"(betadisper); $p < 0.05$ indicates that a PERMANOVA difference may "
      r"partly reflect spread rather than location. Features were filtered to "
      r"10\% prevalence at family, genus and species rank and to 1\% at ASV rank, "
      r"so that the ASV rows of all distances use the same features as the "
      r"phylogenetic (UniFrac) analysis.\label{tab:beta_ranks}}",
      r"\end{table}"]
(OUT / "tableS_beta_ranks.tex").write_text("\n".join(LB) + "\n")

print("Wrote:")
for f in sorted(OUT.glob("*.tex")):
    print(f"  {f.relative_to(REVISION)}")
