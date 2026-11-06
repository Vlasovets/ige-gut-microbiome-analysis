#!/usr/bin/env python
"""
Module 03b - annotate the full-design genus networks with the prevalence core
and the paper's DA results, and quantify where the core sits in the networks.

WHY (R1: "overlay high-confidence core taxa onto the SPIEC-EASI networks")
  "Core" here is the PREVALENCE core (ecological core microbiome), not DA taxa:
  a genus is core if present in >= 80% of the cohort (sensitivity 70% / 90%).
  "High-confidence" core additionally requires >= 80% in BOTH groups of the
  partition, so core membership does not depend on IgE status.

  The DA columns are the paper's primary analysis, unchanged: full cohort, all
  sensitized vs all 233 controls, ~ IgE + sex + age + BMI, genus rank, 5%
  prevalence filter (results/da_full + da_extended_full). LinDA and ALDEx2 are
  the calibrated methods; ANCOM-BC2 calls are kept in the table but flagged,
  because the permutation check (module 02i/02j) showed them at the null level.

  The previous overlay (06b, networks/network_overlay_annotation.csv) marked the
  genus containing the matched-design consensus taxon. That consensus is not
  what the paper reports, so it is superseded by this table.

CAVEAT built into the stats: SPIEC-EASI estimates prevalent taxa from fewer
zeros, so they tend to get more edges. Degree vs prevalence (Spearman) is
reported next to the core comparison so the reader can see how much of any
"core taxa are central" pattern is this artifact.

INPUTS (explicit paths, nothing discovered)
  results/taxa_sets/genus_top50_cohort.csv                  cohort prevalence
  <kora_data>/taxonomy_clean.csv                            genus -> family
  results/full_cohort/inputs/genus/prev10/<part>/source/    group prevalence
  results/networks_full/<part>/edges_<part>.csv             module 03 --design full
  results/da_full/da_long_genus_prev05_<part>.csv
  results/da_extended_full/daext_long_genus_prev05_<part>.csv

OUTPUTS (results/networks_full/)
  network_node_annotation.csv   one row per partition x genus
  core_network_stats.csv        one row per partition x group
  core_differential_edges.csv   FDR-significant differential edges + core status

Usage:  python revision/analysis/03b_network_overlay.py
"""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import mannwhitneyu, spearmanr

REV = Path(__file__).resolve().parent.parent
RES = REV / "results"
NETF = RES / "networks_full"
cfg = yaml.safe_load((REV / "config" / "config.yaml").read_text())
KORA = Path(cfg["paths"]["kora_data"])
PARTS = cfg["networks"]["inference_partitions"]
SEED = cfg["statistics"]["seed"]

CORE_PRIMARY = 0.80
CORE_SENS = [0.70, 0.90]
DA_PREV = "prev05"          # the paper's primary DA filter
ALPHA = 0.05
CALIBRATED = ["LinDA", "ALDEx2"]
UNCALIBRATED = ["ANCOMBC2"]
N_SHUFFLE = 10000

# Reuse 00c's label cleaning so genus names match the frozen manifest exactly.
spec = importlib.util.spec_from_file_location("taxa_sets", REV / "analysis" / "00c_define_taxa_sets.py")
taxa_sets = importlib.util.module_from_spec(spec)
spec.loader.exec_module(taxa_sets)


def genus_to_family(genera):
    """Majority family over ASVs per genus; genera spanning families are flagged."""
    tax = pd.read_csv(KORA / "taxonomy_clean.csv", index_col=0)
    g = tax.genus.map(lambda x: taxa_sets.clean_label(x, "g__"))
    f = tax.family.map(lambda x: taxa_sets.clean_label(x, "f__"))
    t = pd.DataFrame({"genus": g, "family": f}).dropna(subset=["genus"])
    t = t[t.genus.isin(genera)]
    counts = t.groupby(["genus", "family"], dropna=False).size().rename("n").reset_index()
    counts = counts.sort_values(["genus", "n"], ascending=[True, False])
    top = counts.groupby("genus").first()
    n_fam = counts.groupby("genus").family.nunique()
    share = top.n / counts.groupby("genus").n.sum()
    return pd.DataFrame({"family": top.family, "n_families": n_fam,
                         "family_asv_share": share.round(3)})


def group_prevalence(part, genera):
    src = RES / "full_cohort" / "inputs" / "genus" / "prev10" / part / "source"
    out = {}
    for grp, fn in [("ige", "otu_ige.csv"), ("ctrl", "otu_non_ige.csv")]:
        o = pd.read_csv(src / fn, index_col=0)
        missing = set(genera) - set(o.columns)
        assert not missing, f"{part}/{fn} lacks {sorted(missing)}"
        out[f"prev_{grp}"] = (o[genera] > 0).mean(axis=0)
        out[f"n_{grp}"] = len(o)
    return out


def da_wide(part, genera):
    d = pd.concat([pd.read_csv(RES / "da_full" / f"da_long_genus_{DA_PREV}_{part}.csv"),
                   pd.read_csv(RES / "da_extended_full" / f"daext_long_genus_{DA_PREV}_{part}.csv")])
    missing = set(genera) - set(d.taxon)
    assert not missing, f"DA {part} lacks {sorted(missing)}"
    d = d[d.taxon.isin(genera)]
    w = d.pivot_table(index="taxon", columns="method",
                      values=["effect_log2", "p_raw", "q_within"], aggfunc="first")
    w.columns = [f"{m}_{v.replace('_within', '').replace('_raw', '')}" for v, m in w.columns]
    calib = d[d.method.isin(CALIBRATED) & (d.q_within < ALPHA)].groupby("taxon").size()
    uncal = d[d.method.isin(UNCALIBRATED) & (d.q_within < ALPHA)].groupby("taxon").size()
    w["n_sig_calibrated"] = calib.reindex(w.index, fill_value=0)
    w["sig_uncalibrated_only"] = (uncal.reindex(w.index, fill_value=0) > 0) & (w.n_sig_calibrated == 0)
    return w


def degrees(edges, genera, col):
    e = edges[edges[col]]
    deg = pd.concat([e.taxon_i, e.taxon_j]).value_counts()
    return deg.reindex(genera, fill_value=0)


def edge_mix(edges, col, core):
    e = edges[edges[col]]
    ci, cj = e.taxon_i.map(core).values, e.taxon_j.map(core).values
    return int((ci & cj).sum()), int((ci ^ cj).sum()), int((~ci & ~cj).sum())


def main():
    top = pd.read_csv(RES / "taxa_sets" / "genus_top50_cohort.csv", index_col=0)
    genera = list(top.index)
    fam = genus_to_family(genera)
    assert set(fam.index) == set(genera), set(genera) - set(fam.index)
    rng = np.random.default_rng(SEED)

    ann_rows, stat_rows, diff_rows = [], [], []
    for part in PARTS:
        ef = NETF / part / f"edges_{part}.csv"
        if not ef.exists():
            print(f"  {part}: no full-design network yet ({ef}), skipped")
            continue
        edges = pd.read_csv(ef)
        assert (edges.design == "full").all()
        assert set(edges.taxon_i) | set(edges.taxon_j) == set(genera)

        gp = group_prevalence(part, genera)
        ann = pd.DataFrame(index=pd.Index(genera, name="genus"))
        ann["partition"] = part
        ann["phylum"] = edges.drop_duplicates("taxon_i").set_index("taxon_i").phylum_i \
            .reindex(genera).fillna(edges.drop_duplicates("taxon_j").set_index("taxon_j").phylum_j)
        ann = ann.join(fam)
        ann["mean_rel_abundance"] = top.mean_rel_abundance
        ann["prev_cohort"] = top.prevalence
        ann["prev_ige"], ann["prev_ctrl"] = gp["prev_ige"], gp["prev_ctrl"]
        for t in [CORE_PRIMARY] + CORE_SENS:
            ann[f"core{int(t * 100)}"] = ann.prev_cohort >= t
        ann["core80_both_groups"] = (ann.prev_ige >= CORE_PRIMARY) & (ann.prev_ctrl >= CORE_PRIMARY)
        ann["degree_ige"] = degrees(edges, genera, "present_ige")
        ann["degree_ctrl"] = degrees(edges, genera, "present_non_ige")
        ann = ann.join(da_wide(part, genera))
        ann_rows.append(ann.reset_index())

        core = ann.core80
        for grp, col, dcol in [("ige", "present_ige", "degree_ige"),
                               ("ctrl", "present_non_ige", "degree_ctrl")]:
            dc, dn = ann.loc[core, dcol], ann.loc[~core, dcol]
            cc, cn, nn = edge_mix(edges, col, core)
            labels = core.values.copy()
            null_cc = np.empty(N_SHUFFLE)
            for b in range(N_SHUFFLE):
                rng.shuffle(labels)
                null_cc[b] = edge_mix(edges, col, pd.Series(labels, index=core.index))[0]
            rho, prho = spearmanr(ann.prev_cohort, ann[dcol])
            stat_rows.append({
                "partition": part, "group": grp,
                "n_ige": gp["n_ige"], "n_ctrl": gp["n_ctrl"],
                "core_def": f"cohort prevalence >= {CORE_PRIMARY:.0%}",
                "n_core": int(core.sum()), "n_core_both_groups": int(ann.core80_both_groups.sum()),
                "n_edges": cc + cn + nn,
                "mean_degree_core": dc.mean(), "mean_degree_noncore": dn.mean(),
                "mannwhitney_p": mannwhitneyu(dc, dn, alternative="two-sided").pvalue,
                "edges_core_core": cc, "edges_core_noncore": cn, "edges_noncore_noncore": nn,
                "expected_core_core": null_cc.mean(),
                "p_core_core_enrichment": (1 + (null_cc >= cc).sum()) / (N_SHUFFLE + 1),
                "spearman_degree_vs_prevalence": rho, "spearman_p": prho})

        sig = edges[edges.testable & (edges.q_adaptBH < ALPHA)]
        for _, e in sig.iterrows():
            diff_rows.append({"partition": part, "taxon_i": e.taxon_i, "taxon_j": e.taxon_j,
                              "asso_ige": e.asso_ige, "asso_ctrl": e.asso_non_ige,
                              "p_perm": e.p_perm, "q_adaptBH": e.q_adaptBH,
                              "core_i": bool(core[e.taxon_i]), "core_j": bool(core[e.taxon_j])})
        print(f"  {part}: core80 {int(core.sum())}/50 (both groups {int(ann.core80_both_groups.sum())}), "
              f"calibrated DA hits {int((ann.n_sig_calibrated > 0).sum())}, "
              f"uncalibrated-only {int(ann.sig_uncalibrated_only.sum())}, "
              f"differential edges {len(sig)}")

    if not ann_rows:
        raise SystemExit("no full-design networks found; run module 03 --design full first")
    pd.concat(ann_rows).to_csv(NETF / "network_node_annotation.csv", index=False)
    pd.DataFrame(stat_rows).to_csv(NETF / "core_network_stats.csv", index=False)
    pd.DataFrame(diff_rows, columns=["partition", "taxon_i", "taxon_j", "asso_ige", "asso_ctrl",
                                     "p_perm", "q_adaptBH", "core_i", "core_j"]) \
        .to_csv(NETF / "core_differential_edges.csv", index=False)
    print(f"wrote network_node_annotation.csv, core_network_stats.csv, "
          f"core_differential_edges.csv to {NETF}")


if __name__ == "__main__":
    main()
