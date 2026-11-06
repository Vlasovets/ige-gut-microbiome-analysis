# Taxonomy — SILVA 138.2 reclassification

## Assignment and agreement by rank

| rank | assigned_original_pct | assigned_rescript_pct | change_pp | n_both_assigned | agreement_raw_pct | agreement_with_renames_pct |
|---|---|---|---|---|---|---|
| domain | 100 | 100 | 0 | 15170 | 100 | 100 |
| phylum | 99.98 | 99.99 | 0.006592 | 15165 | 30.1 | 99.99 |
| class | 99.97 | 99.99 | 0.01978 | 15163 | 99.66 | 99.66 |
| order | 99.5 | 99.82 | 0.3164 | 15090 | 99.47 | 99.47 |
| family | 98.34 | 99.56 | 1.226 | 14903 | 92.34 | 92.34 |
| genus | 84.4 | 93.8 | 9.394 | 12581 | 70.16 | 76.62 |
| species | 12.43 | 25.02 | 12.59 | 1543 | 98.44 | 98.44 |


---

# Taxonomy: original vs a documented RESCRIPt V3–V4 reference

A **separate analysis**, not a reproduction. `taxonomy_clean.csv` came from a RESCRIPt
classification whose reference construction is recorded nowhere, so rather than guess
the original settings we built one with published best-practice parameters and
reclassified the same 15,170 ASVs.

**Two things differ from the original at once** — the reference construction *and* the
QIIME 2 version (2026.7 vs 2021.11) — so no single difference can be attributed to
either alone. `taxonomy_clean.csv` remains the record of what the manuscript reports.

Pipeline: `revision/slurm/00e_rescript_taxonomy.sh` (job 39237074, **5 h 32 m**, peak
322 GB). SILVA 138.2 SSU NR99 → cull (5 degenerates, 8-mer homopolymers) →
length-filter (Archaea 900 / Bacteria 1200 / Eukaryota 1400) → dereplicate `uniq` →
**extract to V3–V4 with 341F/785R** → dereplicate → naive-Bayes → classify.
Comparison: `revision/analysis/09_taxonomy_comparison.py`.

Cost note: `fit-classifier-naive-bayes` alone took **4 h 29 m** of the 5 h 32 m. The
`uniq` dereplication mode that preserves species-level information also produces a very
large class space, which is what makes training expensive.

## Assignment rates

| Rank | Original | RESCRIPt V3–V4 | Change |
|---|---|---|---|
| Domain / Phylum / Class | 100% | 100% | — |
| Order | 99.5% | 99.8% | +0.3pp |
| Family | 98.3% | 99.6% | +1.2pp |
| **Genus** | 84.4% | **93.8%** | **+9.4pp** |
| Species (raw) | 12.4% | 25.0% | +12.6pp |
| **Species (informative)** | **12.3%** | **16.8%** | **+4.5pp** |

⚠ **The raw species figure is misleading.** SILVA carries placeholder "species" labels
that are not species: `metagenome` (525 ASVs), `gut_metagenome` (521), `human_gut` (116),
`unidentified` (86) — **1,248 in total**, against 19 in the original. Counting them as
assignments inflates 16.8% to 25.0%. Only the informative figure should be quoted.

**Genus is where region extraction genuinely pays off: +9.4pp.** Species barely moves.

## Agreement — raw numbers are not interpretable here

| Rank | Both assigned | Raw agreement | After collapsing renames |
|---|---|---|---|
| Domain | 15,170 | 100.0% | 100.0% |
| **Phylum** | 15,165 | **30.1%** | **100.0%** |
| Class | 15,163 | 99.7% | 99.7% |
| Order | 15,090 | 99.5% | 99.5% |
| Family | 14,903 | 92.3% | 92.3% |
| **Genus** | 12,581 | 70.2% | 76.6% |
| Species | 1,543 | 98.4% | 98.4% |

**Phylum "disagreement" is entirely nomenclature.** SILVA 138.2 adopts the
ICNP-standardised names, so every phylum was renamed:

| Original | SILVA 138.2 | ASVs |
|---|---|---|
| Firmicutes | **Bacillota** | 9,093 |
| Proteobacteria | **Pseudomonadota** | 755 |
| Actinobacteriota | **Actinomycetota** | 428 |
| Desulfobacterota | Thermodesulfobacteriota | 182 |
| Cyanobacteria | Cyanobacteriota | 106 |
| Euryarchaeota | Methanobacteriota | 22 |
| Campilobacterota | Campylobacterota | 11 |
| Chloroflexi | Chloroflexota | 1 |

Exactly **one** ASV moves phylum for a non-nomenclatural reason
(Firmicutes → Bacteroidota). Reporting "30% phylum agreement" without this would be
badly misleading — nothing moved.

**Genus is different in kind, and the 76.6% figure should be read with care.** ~6pp of
the apparent disagreement is documented reclassification, but *Prevotella* →
*Segatella*/*Leyella* is a **split into eight genera**, not a rename (table below), so
"agreement" is not really the right question there — the newer taxonomy is finer, not
contradictory. The rename-collapsing in `09_taxonomy_comparison.py` treats
*Prevotella*-descendant genera as concordant, which is defensible for *nomenclature*
but must never be used to argue that a *Prevotella* result transfers to a *Segatella*
result: those are different feature sets. Most of the remaining ~23% is movement to or
from `Incertae_Sedis` — SILVA 138.2 being more conservative about uncertain placements
(`Clostridia_UCG-014` 450, `Muribaculaceae` 408, `UCG-010` 277,
`Clostridia_vadinBB60_group` 255, `RF39` 236).

## The headline signal: *Prevotella copri* → *Segatella*

This is the most consequential result, and it settles a question open since the audit.

**All 108 ASVs labelled `Prevotella_copri` in the original are assigned to genus
*Segatella* by SILVA 138.2.** At species level the current reference **declines the
call**:

| RESCRIPt species label for those 108 ASVs | n |
|---|---|
| *unassigned* | **101** |
| `metagenome` (a placeholder) | 4 |
| `Prevotella_copri` | 3 |

⚠ **But this does NOT license restating the finding at genus rank.** An earlier draft of
this file recommended exactly that; it was wrong, and the error is worth recording
because it is the same overclaiming the revision exists to correct.

**Two facts block the genus reading.**

**(1) The tested feature is not the genus.** The DA test was run on the 108-ASV
aggregate. Genus *Segatella* contains **494** ASVs. The tested feature is **21.9%** of
it. Nothing was ever tested about the other 386.

**(2) *Prevotella* does not rename to *Segatella* — it splits eight ways.** The 812
original *Prevotella* ASVs land as:

| SILVA 138.2 genus | ASVs |
|---|---|
| *Segatella* | 493 |
| *Leyella* | 193 |
| *Prevotella* (retained) | 47 |
| *Massiliprevotella* | 32 |
| Incertae_Sedis | 25 |
| *Hoylesella* | 9 |
| *unassigned* | 9 |
| *Xylanibacter* | 4 |

So "*Prevotella* → *Segatella*" is a **genus split**, not a rename, and the two genus
aggregates are not comparable features.

**(3) The genus-rank test that does exist is not significant.** From
`results/tables/da_consensus_tiers.csv`, genus *Prevotella* / `any_ige`:

| Prevalence | Tier | Methods sig | LinDA q | ANCOMBC2 q | ALDEx2 q |
|---|---|---|---|---|---|
| prev10 | 2 | ANCOMBC2 only | **0.379** | 3.8e-4 | 0.801 |
| prev05 | 2 | ANCOMBC2 only | **0.345** | 2.8e-4 | 0.837 |
| prev01 | 2 | ANCOMBC2 only | **0.453** | 4.9e-5 | 0.930 |

Only ANCOM-BC2 fires — and ANCOM-BC2 is the method independently shown to be
anti-conservative on this dataset (18.3% of taxa significant, vs 0.1% LinDA and 0.0%
ALDEx2). Compare the species-level feature, which is **Tier 1** on two methods:
LinDA q = 0.0215 / 0.0210 / 0.0245, ANCOMBC2 q ≈ 1e-13.

### What this means for the manuscript

1. **Scope the claim to the feature actually tested, at the rank actually tested.**
   Write: *"an aggregate of 108 ASVs previously reported as* Prevotella copri*, which
   current taxonomy (SILVA 138.2) places in the genus* Segatella*"*. Then report the
   species-level statistics, because that is the test that was run.
2. **Do not write "genus *Segatella* is depleted"** — that generalises to 494 ASVs from
   a test on 108, and the genus-rank test is null under two of three methods.
3. **Do not write "*Segatella copri*"** either — 101 of the 108 ASVs lose the species
   call against a current reference, so the binomial is not supported.
4. **The honest gap, and how to close it.** No test exists for a genuine *Segatella*
   feature. Creating one means re-aggregating ASV counts under the SILVA 138.2 taxonomy
   and re-running the DA — see `REVISION_PLAN.md` §6.4. Until that is run, the
   species-level result on the 108-ASV aggregate is the only defensible statement, and
   the taxonomy note is a **labelling** correction, not a change of the tested unit.
5. **The phylum names in the manuscript are the old nomenclature.** *Firmicutes*,
   *Proteobacteria* and *Actinobacteriota* are now *Bacillota*, *Pseudomonadota* and
   *Actinomycetota*. Either update throughout or state explicitly which nomenclature is
   used and why.

## Verdict on the species-framing question

The analysis was designed to test whether the sparse species assignment is an artifact
of reference choice. **It is not.** With placeholders excluded, a current,
best-practice, region-extracted reference lifts informative species assignment only from
12.3% to 16.8%.

The provisional species-level framing in `REVISION_PLAN.md` §6.4 therefore **stands, and
now rests on a tested alternative rather than an assumption** — a materially stronger
position for the response letter than the argument it replaces.

## Artifacts

`KORA_data/q2_rescript_taxonomy/`: `kora_taxonomy_rescript_v3v4.tsv`,
`silva-138.2-v3v4-nb-classifier.qza`, `silva-v3v4-derep-{seqs,tax}.qza`,
`silva-v3v4-extraction-stats.qza`.
Per-rank and per-ASV comparisons in `revision/results/tables/taxonomy_comparison_*.csv`.


---
_Generated 2026-09-26 by `revision/analysis/10_export_markdown.py`._
