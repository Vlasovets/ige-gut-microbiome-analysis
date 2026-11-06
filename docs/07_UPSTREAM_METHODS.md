# Upstream 16S processing — recovered provenance

Recovered on 2026-08-06 from `KORA_data/preprocessing/` (`KORA_rawseq_preprocessing.ipynb`
and `q2_dada2_KORA.ipynb`), the notebooks that produced `feature_table.tsv` and
`taxonomy_clean.csv`.

This closes the largest documentation gap in the submitted manuscript, which stated
QIIME 2 v2021.11, SILVA 138 and the V3–V4 primers but **gave no DADA2 parameter at all**
and described no phylogeny construction.

## Sequencing

Illumina MiSeq, paired-end **2 × 250 bp** (verified directly from the reads: mean length
250.0 bp in both directions). 2,034 samples, 4,068 FASTQ files.

Primers, V3–V4:

| | Sequence | Length |
|---|---|---|
| Forward, 341F | `CCTACGGGNGGCWGCAG` | 17 nt |
| Reverse, 785R | `GACTACHVGGGTATCTAATCC` | 21 nt |

## Import

```
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path <manifest> \
  --input-format PairedEndFastqManifestPhred33V2
```

Artifact: `paired-end-demux.qza`, UUID `c5a66171-3be9-426e-b6e0-c8113b0b5f8d`,
8.45 GB, created 2021-08-25, held at
`/lustre/groups/bds01/datasets/KORA/paired-end-demux.qza`.

## Denoising — the parameters that were missing

```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trim-left-f 17 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 245 \
  --p-trunc-len-r 245 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-trunc-q 2 \
  --p-n-threads 0
```

Two things worth stating explicitly in the Methods, because they are not obvious:

- **Primers were removed by positional trimming, not by cutadapt.** `--p-trim-left-f 17`
  and `--p-trim-left-r 21` are exactly the primer lengths. There is no separate
  primer-removal step, so the Methods should not imply one.
- **Chimera removal was left at the default**, `--p-chimera-method consensus`
  (de novo, per-sample, consensus across samples). The submitted text says nothing
  about chimera handling.

Yield: **15,170 ASVs**, matching `feature_table.tsv`. The notebook records that an
earlier RDP-based run by a collaborator gave 15,801 ASVs, the difference being the
classifier.

## Taxonomy — RESCRIPt, and what is still undocumented

**`taxonomy_clean.csv` came from the RESCRIPt-based classification** (confirmed by the
authors). The notebook's live cell tabulates
`KORA_rescript_classified_sequences/classification.qza`.

⚠ **Correction to an earlier draft of this file.** It previously attributed the taxonomy
to `classify-sklearn` against the pre-trained full-length `silva-138-99-nb-classifier.qza`,
and drew a caveat from that. **That was wrong.** Those commands appear in the notebook
but are *commented out*, and they write to a different directory
(`KORA_classified_sequences/`). The full-length-classifier caveat is withdrawn.

⚠ **The RESCRIPt commands are NOT in either notebook.** Only the output is consumed, so
the reference construction — which SILVA release, whether sequences were culled by
ambiguity/homopolymer content, whether the reference was dereplicated, and above all
**whether reads were extracted to the V3–V4 region** — cannot be recovered from the
material supplied. These are needed for the Methods and must come from whoever ran
RESCRIPt.

**Because they could not be recovered, they were not guessed at.** Instead a reference
was built from the published best-practice parameters and the same 15,170 ASVs were
reclassified as a *separate analysis* — `revision/slurm/00e_rescript_taxonomy.sh`,
compared in `revision/analysis/09_taxonomy_comparison.py`, written up in
[`TAXONOMY_COMPARISON.md`](TAXONOMY_COMPARISON.md). Headline results are folded into the
sections below. `taxonomy_clean.csv` remains the record of what the manuscript reports.

### Assignment rates actually achieved

Measured directly from `taxonomy_clean.csv` (15,170 ASVs):

| Rank | Assigned | Unassigned |
|---|---|---|
| Domain / Phylum / Class | 100% | 0% |
| Order | 99.5% | 0.5% |
| Family | 98.3% | 1.7% |
| Genus | 84.4% | 15.6% |
| **Species** | **12.4%** | **87.6%** |

Note this does not match Table S2, which reports 96% unknown at species level against
the 87.6% measured here; the discrepancy should be checked before the table is reused.

Sparse species assignment is expected regardless of classifier: SILVA's species
annotations are themselves incomplete, and a V3–V4 amplicon does not carry enough
information to separate many congeners. **This has now been tested rather than
assumed.** A current, culled, length-filtered, V3–V4-extracted SILVA 138.2 reference
raises *informative* species assignment only from **12.3% to 16.8%** (+4.5pp) — the
apparent jump to 25.0% is 1,248 placeholder labels (`metagenome`, `gut_metagenome`,
`human_gut`, `unidentified`) that are not species. Genus assignment does improve
materially, 84.4% → 93.8%.

So the sparsity is inherent to the amplicon and the reference, not an artifact of
reference choice. The species-framing policy in `REVISION_PLAN.md` §6.4 stands, and now
rests on a tested alternative rather than on amplicon-length reasoning alone.

### One detail that matters for the headline signal

***Prevotella_copri* is the single most common species label in the dataset, carried by
108 distinct ASVs.** So the species-level DA result is an aggregate over 108 sequences
rather than one, which makes it more robust than a single-ASV hit — but also means the
label may span several closely related *Prevotella*/*Segatella* lineages that V3–V4
cannot separate. Both points belong in the Results.

**The reclassification analysis confirms exactly this, and sharpens it.** Under SILVA
138.2, **all 108 of those ASVs move to genus *Segatella*** (cohort-wide: *Prevotella*
812 → 47, *Segatella* 0 → 494), but **101 of the 108 lose the species label entirely**
(4 → the placeholder `metagenome`, 3 retain *Prevotella_copri*). So the binomial is not
supported by a current reference and should be dropped.

⚠ **That is a licence to drop the species label, NOT to restate the result at genus
rank.** An earlier draft of this file said otherwise and was wrong. Genus *Segatella*
holds **494** ASVs — the tested feature is 21.9% of it; *Prevotella* **splits eight
ways** rather than renaming (Segatella 493, Leyella 193, Prevotella 47,
Massiliprevotella 32, Incertae_Sedis 25, Hoylesella 9, unassigned 9, Xylanibacter 4);
and the genus-rank test is **null** (*Prevotella*/`any_ige`: LinDA q = 0.379/0.345/0.453,
ALDEx2 q ≈ 0.80–0.93, ANCOM-BC2-only). Name the tested aggregate instead — see
[`TAXONOMY_COMPARISON.md`](TAXONOMY_COMPARISON.md) and `REVISION_PLAN.md` §6.4.

⚠ SILVA 138.2 also adopts the ICNP phylum names — *Firmicutes* → **Bacillota**,
*Proteobacteria* → **Pseudomonadota**, *Actinobacteriota* → **Actinomycetota**. If the
newer reference is cited anywhere in the revision, the manuscript's phylum nomenclature
must be updated to match or the choice stated explicitly.

Other frequent species labels: *Bacteroides_vulgatus* (76 ASVs), *B. plebeius* (51),
*B. massiliensis* (49), *Ruminococcus_bicirculans* (47), *Alistipes_shahii* (41).

## Phylogeny — it was built, and the command is recorded

```
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences representative_sequences.qza \
  --output-dir tree
```

That pipeline is MAFFT alignment → conservation/gap masking → FastTree → midpoint
rooting, producing `rooted_tree.qza`, exported to
`data/q2/tree/rooted_tree_15170/tree.nwk` — the path notebook 02 reads for Faith's PD.

**REBUILT AND VERIFIED.** Job 39166540 (1 h 49 m, 32 cores) re-ran DADA2 from
`paired-end-demux.qza` and rebuilt the tree. The result reproduces the published data
**exactly**:

| Check | Result |
|---|---|
| ASV identifiers | **15,170 / 15,170 identical** |
| Tree tips | 15,170, rooted, with branch lengths |
| Count matrix | **2 cells differ of 30,855,780** — one read swapped between two ASVs in sample 85661 (185↔186); all sample and grand totals unchanged |
| Grand total reads | 28,526,232 in both |

Because an ASV id is the MD5 of its sequence, 15,170 identical ids means 15,170 identical
sequences — which **independently proves the recovered parameters above are correct**.
Wrong truncation, trimming or filtering could not produce a matching set.

Environment recovered in passing: **R 4.0.5, DADA2 1.18.0, Rcpp 1.0.7**.

Artifacts: `KORA_data/q2_reprocessing/` (rooted_tree.nwk, rep_seqs.fasta,
representative_sequences.qza, table.qza, alignment/masked_alignment, denoising_stats).

The rebuild script is `revision/slurm/00d_rebuild_phylogeny.sh`. It pins
**qiime2-2021.11** on purpose and verifies the ASV overlap before endorsing the tree,
refusing to proceed if coverage is poor.

Blocker B2 is closed: **UniFrac and Faith's PD are computed** in
`revision/analysis/05b_phylo_diversity.R`. Both are null — see `REVISION_PLAN.md` §5.6.

### A gotcha for anyone re-running this

`qiime2-2021.11` lives at `/lustre/groups/shared/apps/conda/envs/` after an env
relocation. The `R` launcher is a shell script and was prefix-fixed during that move,
but **`Rscript` is a compiled binary with the original prefix baked in**
(`/home/icb/oleg.vlasovetc/miniconda3/envs/qiime2-2021.11/lib/R`). DADA2 therefore dies
with `Rscript execution error: No such file or directory`, and exporting `R_HOME` does
not help because Rscript fails before reading it. The fix is to repopulate the legacy
prefix with symlinks to the relocated env. **`qiime2-2021.8`, `qiime2-2022.2` and
`q2_ll` carry the same damage.**

## Export to plain formats

```
qiime tools export --input-path table.qza --output-path phyloseq
biom convert -i feature-table.biom -o asv_table.tsv --to-tsv
```

`taxonomy_clean.csv` is a cleaned derivative of the classifier output; the cleaning
script is still not in any repository.

## What remains unrecoverable from these notebooks

- **The entire RESCRIPt reference-construction workflow.** Confirmed as the source of
  `taxonomy_clean.csv`, but no RESCRIPt command appears in either notebook. Needed for
  the Methods: SILVA release and point version, culling/dereplication settings, whether
  reads were extracted to the V3–V4 region, and the classifier training call.
- The `--p-confidence` used at classification (not recorded; QIIME 2 default is 0.7).
- Why Table S2 reports 96% unknown at species level where `taxonomy_clean.csv` gives
  87.6%.
- The cleaning step that turned the raw classifier output into `taxonomy_clean.csv`.
