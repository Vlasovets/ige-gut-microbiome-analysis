# KEY RESOURCES TABLE (generated from the live environment)

Repository: `revision-r1` @ `5f5c5a562c`  **uncommitted changes present**

| Category | Reagent or resource | Source | Version / identifier | Used for |
|---|---|---|---|---|
| Software and algorithms | QIIME 2 | Bolyen et al. 2019 | see manuscript / to be confirmed | 16S read processing (upstream; version from manuscript: 2021.11) |
| Software and algorithms | DADA2 (q2-dada2) | Callahan et al. 2016 | see manuscript / to be confirmed | denoising - ABSENT from the submitted KRT |
| Software and algorithms | SILVA SSU reference | Quast et al. 2012 | see manuscript / to be confirmed | taxonomy reference (v138; release/point version to be confirmed) |
| Software and algorithms | RESCRIPt | Robeson et al. 2021 | see manuscript / to be confirmed | built the SILVA reference actually used for taxonomy - the pre-trained full-length classifier in the notebook is commented out and was NOT used |
| Software and algorithms | R | R Core Team | 4.4.3 | statistical analysis |
| Software and algorithms | Python | Python Software Foundation | 3.11.15 | analysis and machine learning |
| Software and algorithms | LinDA | Zhou et al. 2022 | 0.2.0 (HEAD @ af0f62fad8) | differential abundance |
| Software and algorithms | ANCOM-BC2 | Lin & Peddada 2020 | 2.8.0 | differential abundance (consensus) - ABSENT from the submitted KRT |
| Software and algorithms | ALDEx2 | Fernandes et al. 2014 | 1.38.0 | differential abundance (consensus) - ABSENT from the submitted KRT |
| Software and algorithms | SPIEC-EASI | Kurtz et al. 2015 | 1.1.1 (lowrank @ 0d347382aa) | network inference - ABSENT from the submitted KRT |
| Software and algorithms | NetCoMi | Peschel et al. 2021 | 1.2.0 (develop @ 4598486d53) | network construction and comparison - ABSENT from the submitted KRT |
| Software and algorithms | pulsar | Muller et al. 2016 | 0.3.13 | StARS model selection - ABSENT from the submitted KRT |
| Software and algorithms | vegan | Oksanen et al. | 2.7.5 | PERMANOVA, betadisper, distances |
| Software and algorithms | phyloseq | McMurdie & Holmes 2013 | 1.50.0 | microbiome data handling - ABSENT from the submitted KRT |
| Software and algorithms | TaxSEA | Pham et al. 2025 | 1.2.0 | taxon set enrichment |
| Software and algorithms | c-lasso | Simpson et al. 2021 | 1.0.11 | log-contrast classification - THE classifier used; ABSENT from the KRT |
| Software and algorithms | scikit-learn | Pedregosa et al. 2011 | 1.8.0 | predictive benchmark, cross-validation |
| Software and algorithms | statsmodels | Seabold & Perktold 2010 | 0.14.6 | multiple-testing correction |
| Software and algorithms | GUniFrac | Chen et al. 2012 | 1.9 | unweighted, generalised and weighted UniFrac (module 05b) |
| Software and algorithms | ape | Paradis & Schliep 2019 | see manuscript / to be confirmed | phylogeny handling and Faith's PD |
| Software | ZicoSeq (permutation-based DA) | GUniFrac (CRAN) | 1.9 | extended DA consensus (module 02d) |
| Software | limma-voom | Bioconductor | 3.62.2 | extended DA consensus (module 02d) |
| Software | XGBoost | revision/env/xgb-venv (pip) | 2.1.4 | ML benchmark (modules 04/04b/04c) |
| Software | QIIME 2 + RESCRIPt (independent reference) | conda q2-2026.7-slr | 2026.7; SILVA 138.2 SSU NR99, V3-V4-extracted (SLURM job 39237074) | taxonomy comparison (module 09) |
| Software | legacy Python (published-notebook semantics) | revision/env/legacy-venv (pip) | pandas==2.2.3, numpy==1.26.4 | published-design reproduction (module 07a) |

## Removed from the submitted Key Resources Table

- **GGLasso** — listed in the submitted KRT but used nowhere in the study

## Still required before submission

- **DADA2 parameters: RECOVERED AND VERIFIED.** trim-left 17/21 (primer
  lengths; primers removed by trimming, not cutadapt), trunc-len 245/245,
  max-EE 2/2, trunc-q 2, chimera consensus, under QIIME 2 2021.11. The re-run
  (SLURM job 39166540) reproduces the published feature table: 15,170/15,170
  ASV identifiers identical, 2 count cells of 30,855,780 differ.
- **UniFrac and Faith's PD are reported.** The phylogeny was rebuilt from the
  source reads and reproduces the published ASV set exactly (15,170 / 15,170
  identifiers; 2 count cells differ of 30,855,780). Both metrics are null.
- **The ORIGINAL RESCRIPt reference construction remains undocumented** (no
  RESCRIPt call in the preprocessing notebooks). Handled without guessing:
  module 09 built an independent, fully documented SILVA 138.2 V3-V4
  reference (see revision/TAXONOMY_COMPARISON.md) and compared; informative
  species assignment rises only 12.3% -> 16.8%, so species sparsity is
  inherent to the amplicon, not the reference choice.
- **No dietary variable exists.** `u3su_d3*` is a self-reported intolerance
  questionnaire, not intake, so diet cannot be adjusted for.
- **Bioconductor 3.20** pins the main environment; TaxSEA requires R >= 4.5 and
  therefore runs in a separate `taxsea-env`.
