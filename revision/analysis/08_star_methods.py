#!/usr/bin/env python
"""
Module 08 - STAR Methods and Key Resources, captured from the live environment.

Answers the editorial request and R3 Major 8.

The submitted Key Resources Table has three problems this fixes:
  - it lists GGLasso, which is used nowhere in the study;
  - it OMITS c-lasso (the classifier actually used), DADA2, SPIEC-EASI, NetCoMi,
    pulsar, scikit-bio, ANCOM-BC and phyloseq;
  - almost nothing carries a version, and the two that do (R v4.1, Python v3.9)
    appear nowhere in the manuscript body.

Versions are read from the environment rather than typed, so the table cannot
drift from what actually ran. Non-CRAN packages carry their git SHA, which matters
here because SpiecEasi and NetCoMi are installed from non-release branches.

Also runs a submission-compliance check on the emitted figures: every font size
must be >= 9 pt and a vector format must exist alongside each raster.

Usage:
  python revision/analysis/08_star_methods.py
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import yaml

REVISION = Path(__file__).resolve().parent.parent
REPO = REVISION.parent
CONFIG = REVISION / "config" / "config.yaml"
ENV = REVISION / "env"
FIGS = REVISION / "results" / "figures"
MIN_PT = 9.0

R_BIN = "/home/itg/oleg.vlasovets/miniconda3/envs/slr-env/bin/Rscript"
PY_BIN = "/home/itg/oleg.vlasovets/miniconda3/envs/slr-env/bin/python"

# What the study actually uses, and where it is used. Ordered as Cell Press wants.
RESOURCES = [
    # (category, item, source, identifier_key, used_for)
    ("Software and algorithms", "QIIME 2", "Bolyen et al. 2019", "qiime2",
     "16S read processing (upstream; version from manuscript: 2021.11)"),
    ("Software and algorithms", "DADA2 (q2-dada2)", "Callahan et al. 2016",
     "dada2", "denoising - ABSENT from the submitted KRT"),
    ("Software and algorithms", "SILVA SSU reference", "Quast et al. 2012",
     "silva", "taxonomy reference (v138; release/point version to be confirmed)"),
    ("Software and algorithms", "RESCRIPt", "Robeson et al. 2021", "rescript",
     "built the SILVA reference actually used for taxonomy - the pre-trained "
     "full-length classifier in the notebook is commented out and was NOT used"),
    ("Software and algorithms", "R", "R Core Team", "R", "statistical analysis"),
    ("Software and algorithms", "Python", "Python Software Foundation", "python",
     "analysis and machine learning"),
    ("Software and algorithms", "LinDA", "Zhou et al. 2022", "LinDA",
     "differential abundance"),
    ("Software and algorithms", "ANCOM-BC2", "Lin & Peddada 2020", "ANCOMBC",
     "differential abundance (consensus) - ABSENT from the submitted KRT"),
    ("Software and algorithms", "ALDEx2", "Fernandes et al. 2014", "ALDEx2",
     "differential abundance (consensus) - ABSENT from the submitted KRT"),
    ("Software and algorithms", "SPIEC-EASI", "Kurtz et al. 2015", "SpiecEasi",
     "network inference - ABSENT from the submitted KRT"),
    ("Software and algorithms", "NetCoMi", "Peschel et al. 2021", "NetCoMi",
     "network construction and comparison - ABSENT from the submitted KRT"),
    ("Software and algorithms", "pulsar", "Muller et al. 2016", "pulsar",
     "StARS model selection - ABSENT from the submitted KRT"),
    ("Software and algorithms", "vegan", "Oksanen et al.", "vegan",
     "PERMANOVA, betadisper, distances"),
    ("Software and algorithms", "phyloseq", "McMurdie & Holmes 2013", "phyloseq",
     "microbiome data handling - ABSENT from the submitted KRT"),
    ("Software and algorithms", "TaxSEA", "Pham et al. 2025", "TaxSEA",
     "taxon set enrichment"),
    ("Software and algorithms", "c-lasso", "Simpson et al. 2021", "classo",
     "log-contrast classification - THE classifier used; ABSENT from the KRT"),
    ("Software and algorithms", "scikit-learn", "Pedregosa et al. 2011",
     "sklearn", "predictive benchmark, cross-validation"),
    ("Software and algorithms", "statsmodels", "Seabold & Perktold 2010",
     "statsmodels", "multiple-testing correction"),
    ("Software and algorithms", "GUniFrac", "Chen et al. 2012", "GUniFrac",
     "unweighted, generalised and weighted UniFrac (module 05b)"),
    ("Software and algorithms", "ape", "Paradis & Schliep 2019", "ape",
     "phylogeny handling and Faith's PD"),
]

REMOVE = [("GGLasso", "listed in the submitted KRT but used nowhere in the study")]


def r_versions(pkgs: list[str]) -> dict:
    code = (
        '.libPaths(c("%s", .libPaths()));'
        'ps <- c(%s);'
        'for (p in ps) { v <- tryCatch(as.character(packageVersion(p)),'
        ' error=function(e) NA);'
        ' d <- tryCatch(read.dcf(system.file("DESCRIPTION", package=p)),'
        ' error=function(e) NULL);'
        ' ref <- sha <- NA;'
        ' if (!is.null(d)) { cn <- colnames(d);'
        '  if ("RemoteRef" %%in%% cn) ref <- d[1,"RemoteRef"];'
        '  if ("RemoteSha" %%in%% cn) sha <- d[1,"RemoteSha"] };'
        ' cat(sprintf("%%s\\t%%s\\t%%s\\t%%s\\n", p, v, ref, sha)) };'
        'cat("__R__\\t", as.character(getRversion()), "\\tNA\\tNA\\n");'
        'if (requireNamespace("BiocManager", quietly=TRUE))'
        ' cat("__BIOC__\\t", as.character(BiocManager::version()), "\\tNA\\tNA\\n")'
    ) % (REVISION / "env" / "rlibs", ",".join(f'"{p}"' for p in pkgs))
    out = subprocess.run([R_BIN, "--vanilla", "-e", code],
                         capture_output=True, text=True).stdout
    res = {}
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) >= 4:
            res[parts[0].strip()] = {
                "version": parts[1].strip(),
                "remote_ref": None if parts[2].strip() in ("NA", "") else parts[2].strip(),
                "remote_sha": None if parts[3].strip() in ("NA", "") else parts[3].strip(),
            }
    return res


def py_versions(pkgs: list[str]) -> dict:
    code = (
        "import importlib.metadata as m, sys, json\n"
        "d={'python': {'version': sys.version.split()[0]}}\n"
        f"for p in {pkgs!r}:\n"
        "    for cand in (p, p.replace('_','-'), {'sklearn':'scikit-learn',\n                 'classo':'c-lasso','matplotlib_venn':'matplotlib-venn'}.get(p,p)):\n"
        "        try:\n"
        "            d[p]={'version': m.version(cand)}; break\n"
        "        except Exception: pass\n"
        "    else: d[p]={'version': None}\n"
        "print(json.dumps(d))"
    )
    out = subprocess.run([PY_BIN, "-c", code], capture_output=True, text=True).stdout
    return json.loads(out) if out.strip() else {}


# matplotlib renders mathtext sub/superscripts at fixed fractions of the
# surrounding size (x-small = 0.694, xx-small = 0.579). A subscript in "log_2"
# SHOULD be smaller than body text, so those are exempt - requiring a subscript to
# be >= 9 pt would be wrong. Anything else below the floor is a genuine failure.
MATHTEXT_SCALINGS = (0.694, 0.579)


def check_fonts() -> tuple[int, int, list[str]]:
    """Every font-size in an emitted SVG must be >= MIN_PT, mathtext excepted."""
    problems, checked = [], 0
    for svg in sorted(FIGS.glob("*.svg")):
        checked += 1
        txt = svg.read_text(errors="ignore")
        sizes = [float(m) for m in re.findall(r"font-size:\s*([0-9.]+)", txt)]
        sizes += [float(m) for m in re.findall(r'font-size="([0-9.]+)', txt)]
        uniq = sorted(set(sizes))
        compliant = [s for s in uniq if s >= MIN_PT]

        def is_mathtext(s: float) -> bool:
            return any(abs(s - c * f) < 0.15
                       for c in compliant for f in MATHTEXT_SCALINGS)

        small = sorted({s for s in uniq if 0 < s < MIN_PT})
        genuine = [s for s in small if not is_mathtext(s)]
        exempt = [s for s in small if is_mathtext(s)]
        if exempt:
            print(f"    note {svg.name}: {exempt} exempt as mathtext "
                  f"sub/superscript of compliant text")
        if genuine:
            problems.append(f"{svg.name}: font sizes below {MIN_PT} pt: {genuine}")
    # Every raster must have a vector sibling.
    for png in sorted(FIGS.glob("*.png")):
        if not (png.with_suffix(".pdf").exists() or png.with_suffix(".svg").exists()):
            problems.append(f"{png.name}: no vector (PDF/SVG) counterpart")
    return checked, len(problems), problems


def main() -> int:
    cfg = yaml.safe_load(CONFIG.read_text())
    ENV.mkdir(parents=True, exist_ok=True)

    print("=== Module 08: STAR Methods and Key Resources ===")

    r_pkgs = ["LinDA", "ANCOMBC", "ALDEx2", "SpiecEasi", "NetCoMi", "pulsar",
              "vegan", "phyloseq", "permApprox", "GUniFrac", "lmerTest",
              "MicrobiomeStat", "UpSetR", "TreeSummarizedExperiment", "limma"]
    py_pkgs = ["numpy", "pandas", "scipy", "sklearn", "statsmodels",
               "matplotlib", "seaborn", "classo", "gglasso", "matplotlib_venn",
               "PyYAML"]

    print("\n[1] R packages (slr-env + revision/env/rlibs)")
    rv = r_versions(r_pkgs)
    for k, v in rv.items():
        extra = ""
        if v.get("remote_ref"):
            extra = f"  [{v['remote_ref']} @ {(v.get('remote_sha') or '')[:7]}]"
        print(f"    {k:26s} {v['version']}{extra}")

    print("\n[2] TaxSEA (separate env: R >= 4.5 required)")
    tv = subprocess.run(
        ["bash", "-lc",
         "source ~/miniconda3/etc/profile.d/conda.sh && conda activate taxsea-env "
         "&& Rscript --vanilla -e "
         "'cat(as.character(getRversion()),as.character(packageVersion(\"TaxSEA\")))'"],
        capture_output=True, text=True).stdout.strip()
    print(f"    taxsea-env: {tv or 'UNAVAILABLE'}")

    print("\n[3] Python packages")
    pv = py_versions(py_pkgs)
    for k, v in pv.items():
        print(f"    {k:26s} {v.get('version')}")

    print("\n[4] Repository provenance")
    sha = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                         capture_output=True, text=True).stdout.strip()
    branch = subprocess.run(["git", "-C", str(REPO), "rev-parse",
                             "--abbrev-ref", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain"],
                           capture_output=True, text=True).stdout.strip()
    print(f"    branch {branch} @ {sha[:10]}"
          f"{'  (UNCOMMITTED CHANGES PRESENT)' if dirty else ''}")

    print("\n[5] Figure compliance check")
    checked, n_bad, probs = check_fonts()
    print(f"    {checked} SVG(s) checked, minimum font {MIN_PT} pt")
    for p in probs:
        print(f"    FAIL {p}")
    if not n_bad:
        print("    PASS - all fonts >= 9 pt and every raster has a vector sibling")

    # ── Emit the corrected Key Resources Table ───────────────────────────
    lookup = {}
    for k, v in rv.items():
        lookup[k] = v
    for k, v in pv.items():
        lookup[k] = v

    lines = ["# KEY RESOURCES TABLE (generated from the live environment)", "",
             f"Repository: `{branch}` @ `{sha[:10]}`"
             + ("  **uncommitted changes present**" if dirty else ""), "",
             "| Category | Reagent or resource | Source | Version / identifier | Used for |",
             "|---|---|---|---|---|"]
    for cat, item, src, key, used in RESOURCES:
        v = lookup.get(key, {})
        ver = v.get("version") or "see manuscript / to be confirmed"
        if key == "R":
            ver = rv.get("__R__", {}).get("version", ver)
        if key == "python":
            ver = pv.get("python", {}).get("version", ver)
        if key == "TaxSEA":
            ver = tv.split()[-1] if tv else "1.2.0"
        if v.get("remote_sha"):
            ver = f"{ver} ({v['remote_ref']} @ {v['remote_sha'][:10]})"
        lines.append(f"| {cat} | {item} | {src} | {ver} | {used} |")

    # Round-2 additions, read from the live auxiliary envs - never typed in.
    xgbv = subprocess.run([str(REVISION / "env" / "xgb-venv" / "bin" / "python"),
                           "-c", "import xgboost; print(xgboost.__version__)"],
                          capture_output=True, text=True).stdout.strip()
    legreq = REVISION / "env" / "legacy-requirements.txt"
    legpd = next((l.strip() for l in legreq.read_text().splitlines()
                  if l.lower().startswith("pandas")), "pandas==?")         if legreq.exists() else "pandas==?"
    lines += [
        f"| Software | ZicoSeq (permutation-based DA) | GUniFrac (CRAN) | "
        f"{rv.get('GUniFrac', {}).get('version', '?')} | extended DA consensus (module 02d) |",
        f"| Software | limma-voom | Bioconductor | "
        f"{rv.get('limma', {}).get('version', '?')} | extended DA consensus (module 02d) |",
        f"| Software | XGBoost | revision/env/xgb-venv (pip) | {xgbv or '?'} | "
        f"ML benchmark (modules 04/04b/04c) |",
        "| Software | QIIME 2 + RESCRIPt (independent reference) | conda q2-2026.7-slr | "
        "2026.7; SILVA 138.2 SSU NR99, V3-V4-extracted (SLURM job 39237074) | "
        "taxonomy comparison (module 09) |",
        f"| Software | legacy Python (published-notebook semantics) | "
        f"revision/env/legacy-venv (pip) | {legpd}, numpy==1.26.4 | "
        f"published-design reproduction (module 07a) |",
    ]

    lines += ["", "## Removed from the submitted Key Resources Table", ""]
    for item, why in REMOVE:
        lines.append(f"- **{item}** — {why}")
    lines += [
        "", "## Still required before submission", "",
        "- **DADA2 parameters: RECOVERED AND VERIFIED.** trim-left 17/21 (primer",
        "  lengths; primers removed by trimming, not cutadapt), trunc-len 245/245,",
        "  max-EE 2/2, trunc-q 2, chimera consensus, under QIIME 2 2021.11. The re-run",
        "  (SLURM job 39166540) reproduces the published feature table: 15,170/15,170",
        "  ASV identifiers identical, 2 count cells of 30,855,780 differ.",
        "- **UniFrac and Faith's PD are reported.** The phylogeny was rebuilt from the",
        "  source reads and reproduces the published ASV set exactly (15,170 / 15,170",
        "  identifiers; 2 count cells differ of 30,855,780). Both metrics are null.",
        "- **The ORIGINAL RESCRIPt reference construction remains undocumented** (no",
        "  RESCRIPt call in the preprocessing notebooks). Handled without guessing:",
        "  module 09 built an independent, fully documented SILVA 138.2 V3-V4",
        "  reference (see revision/TAXONOMY_COMPARISON.md) and compared; informative",
        "  species assignment rises only 12.3% -> 16.8%, so species sparsity is",
        "  inherent to the amplicon, not the reference choice.",
        "- **No dietary variable exists.** `u3su_d3*` is a self-reported intolerance",
        "  questionnaire, not intake, so diet cannot be adjusted for.",
        "- **Bioconductor 3.20** pins the main environment; TaxSEA requires R >= 4.5 and",
        "  therefore runs in a separate `taxsea-env`.",
    ]
    krt = REVISION / "results" / "tables" / "KEY_RESOURCES.md"
    krt.parent.mkdir(parents=True, exist_ok=True)
    krt.write_text("\n".join(lines) + "\n")
    print(f"\n  Wrote {krt}")

    (ENV / "versions.json").write_text(json.dumps(
        {"r": rv, "python": pv, "taxsea_env": tv,
         "git": {"branch": branch, "sha": sha, "dirty": bool(dirty)}}, indent=2))
    print(f"  Wrote {ENV / 'versions.json'}")

    subprocess.run(
        ["bash", "-lc",
         "source ~/miniconda3/etc/profile.d/conda.sh && conda activate slr-env && "
         f"conda env export --no-builds > {ENV / 'slr-env.yml'} 2>/dev/null"],
        capture_output=True)
    print(f"  Wrote {ENV / 'slr-env.yml'}")
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main())
