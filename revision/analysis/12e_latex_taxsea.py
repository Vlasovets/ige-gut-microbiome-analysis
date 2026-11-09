#!/usr/bin/env python
"""
Module 12e - supplementary table of the re-computed taxon set enrichment
(TaxSEA) on the full-cohort, covariate-adjusted ranking. Replaces the three
submitted TaxSEA tables (vitamin A, folic acid, SCFA), whose FDR column could
not be produced by Benjamini-Hochberg (FDR < p in four of five rows).

Source: results/tables/taxsea_recomputed_full.csv (module 07 --design full):
species ranked by the signed LinDA statistic of all sensitized vs all 233
non-sensitized individuals adjusted for sex, age and BMI; TaxSEA 1.2.0, bundled
database, set size 5-100. Shown: the primary analysis (5% prevalence filter,
LinDA statistic), the 10 sets with the smallest p-value plus the folic acid and
vitamin A producer sets; FDR as returned by TaxSEA within each sub-database. The
caption summarizes all six analyses (3 prevalence filters x 2 ranking statistics).

Output: results/tables/latex/tableS_taxsea_recomputed.tex (tab:taxsea_recomputed)
Usage:  python revision/analysis/12e_latex_taxsea.py
"""

from pathlib import Path

import pandas as pd

REV = Path(__file__).resolve().parent.parent
TAB = REV / "results" / "tables"
HEADLINE = ["MiMeDB_producers_of_Folic_acid", "MiMeDB_producers_of_Vitamin_A"]
DB = {"BacDive_bacterial_physiology": "Bacterial physiology", "Metabolite_producers": "Metabolite producers",
      "Health_associations": "Health associations", "Gut_Brain_Modules_VallesColomer2019": "Gut-brain modules"}

t = pd.read_csv(TAB / "taxsea_recomputed_full.csv")
sub = t[t.sub_database != "All_databases"]          # FDR within each sub-database
p = sub[(sub.prevalence == "prev05") & (sub.rank_by == "stat")]
show = pd.concat([p.nsmallest(10, "PValue"), p[p.taxonSetName.isin(HEADLINE)]]).drop_duplicates("taxonSetName")
n_sets = t[t.sub_database == "All_databases"].groupby(["prevalence", "rank_by"]).size()


def tex(s):
    return s.replace("_", " ").replace("&", r"\&")


def name(s):
    """Readable set name with its source, e.g. 'Acetate producers (gutMGene)'."""
    import re
    for pre, src in (("MiMeDB_producers_of_", "MiMeDB"), ("GutMGene_producers_of_", "gutMGene")):
        if s.startswith(pre):
            return tex(s[len(pre):]) + f" producers ({src})"
    if s.startswith("GMRepoV2_"):
        v = tex(s[len("GMRepoV2_"):])
        return v[0].upper() + v[1:]
    if s.startswith("BacDive_Enzyme_"):
        return "Enzyme: " + tex(s[len("BacDive_Enzyme_"):])
    if s.startswith("Valles-Colomer2019_"):
        v = re.sub(r"\.+", " ", s[len("Valles-Colomer2019_"):]).replace("( ", "(").replace(" )", ")")
        return tex(v.strip())
    return tex(s)


L = [r"\begin{table}[h!]", r"\centering", r"\footnotesize",
     r"\begin{tabular}{l l r r r}",
     r"\textbf{Taxon set} & \textbf{Database} & \textbf{Median rank} & \textbf{$p$} & \textbf{FDR} \\",
     r"\hline"]
for _, r in show.iterrows():
    L.append(" & ".join([name(r.taxonSetName), DB.get(r.sub_database, tex(r.sub_database)),
                         f"{r.median_rank_of_set_members:.3f}", f"{r.PValue:.3f}",
                         f"{r.FDR:.2f}"]) + r" \\")
L += [r"\hline", r"\end{tabular}",
      r"\caption{Taxon set enrichment analysis (TaxSEA) for all 275 sensitized vs.\ 233 "
      r"non-sensitized individuals. Species were ranked by the signed \texttt{LinDA} test "
      r"statistic from the differential abundance model adjusted for sex, age and BMI (5\% "
      r"prevalence filter, " + f"{int(p.n_input_species.iloc[0])} species, "
      f"{int(p.n_mapped.iloc[0])} mapped to the TaxSEA database" + r"). Shown are the 10 sets with "
      r"the smallest $p$-value and the folic acid and vitamin A producer sets. Median rank: median "
      r"ranking statistic of the set members (positive, higher in sensitized individuals). FDR: "
      r"Benjamini--Hochberg within each database, as returned by TaxSEA. Across all six analyses "
      r"(10\%, 5\% and 1\% prevalence filters; ranking by test statistic or log$_2$ fold change; "
      f"{n_sets.min()}--{n_sets.max()}" + r" sets each), no set reached FDR $< 0.1$ (smallest FDR "
      f"{sub.FDR.min():.2f}" + r").\label{tab:taxsea_recomputed}}",
      r"\end{table}"]
(TAB / "latex" / "tableS_taxsea_recomputed.tex").write_text("\n".join(L) + "\n")
print(show[["taxonSetName", "sub_database", "PValue", "FDR"]].round(3).to_string(index=False))
print("sets per analysis:", n_sets.to_dict(), "| min FDR", round(sub.FDR.min(), 3))
