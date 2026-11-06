# Differential abundance — six methods

The direct answer to *“use more than one DA method, and show whether they agree.”*

All six ran on **identical input** — the same matched, prevalence-filtered tables — with the same formula `~ ige_status + sex + age + bmi`, except the paired Wilcoxon, which uses the 1:1 matching directly instead of adjusting.

| Family | Methods |
|---|---|
| Log-ratio | LinDA, ANCOM-BC2, ALDEx2 |
| Permutation / reference-frame | ZicoSeq |
| Non-parametric, matched | Wilcoxon on CLR (paired) |
| Count model | limma-voom |

## Per-method significance rate

| method | n_tests | n_sig | pct_sig | min_q |
|---|---|---|---|---|
| ALDEx2 | 41405 | 0 | 0 | 0.08224 |
| ANCOMBC2 | 21836 | 6880 | 31.51 | 5.805e-19 |
| LinDA | 46080 | 9 | 0.01953 | 0.02102 |
| Wilcoxon_CLR_paired | 46080 | 0 | 0 | 0.05382 |
| ZicoSeq | 29762 | 4 | 0.01344 | 0.0272 |
| limma_voom | 46080 | 1119 | 2.428 | 1.748e-48 |


A method calling **18%** of taxa significant while two others call 0.0–0.1% is a statement about calibration, not biology.


## Pairwise overlap of significant calls

Restricted to taxon-analyses **both** methods actually evaluated. ZicoSeq drops features that are constant within a partition and ANCOM-BC2 applies internal filtering, so a taxon a method never tested must not be scored as that method finding it non-significant.

| method_a | method_b | n_common_evaluated | n_sig_a_in_common | n_sig_b_in_common | n_shared | n_union | jaccard |
|---|---|---|---|---|---|---|---|
| ANCOMBC2 | limma_voom | 21836 | 6880 | 905 | 447 | 7338 | 0.06092 |
| LinDA | limma_voom | 46080 | 9 | 1119 | 3 | 1125 | 0.002667 |
| ANCOMBC2 | LinDA | 21836 | 6880 | 9 | 3 | 6886 | 0.0004357 |
| ANCOMBC2 | ZicoSeq | 21372 | 6506 | 4 | 2 | 6508 | 0.0003073 |
| ALDEx2 | Wilcoxon_CLR_paired | 41405 | 0 | 0 | 0 | 0 |  |
| ALDEx2 | LinDA | 41405 | 0 | 9 | 0 | 9 | 0 |
| ALDEx2 | ANCOMBC2 | 21836 | 0 | 6880 | 0 | 6880 | 0 |
| ANCOMBC2 | Wilcoxon_CLR_paired | 21836 | 6880 | 0 | 0 | 6880 | 0 |
| ALDEx2 | limma_voom | 41405 | 0 | 1119 | 0 | 1119 | 0 |
| ALDEx2 | ZicoSeq | 29762 | 0 | 4 | 0 | 4 | 0 |
| LinDA | Wilcoxon_CLR_paired | 46080 | 9 | 0 | 0 | 9 | 0 |
| LinDA | ZicoSeq | 29762 | 9 | 4 | 0 | 13 | 0 |
| Wilcoxon_CLR_paired | ZicoSeq | 29762 | 0 | 4 | 0 | 4 | 0 |
| Wilcoxon_CLR_paired | limma_voom | 46080 | 0 | 1119 | 0 | 1119 | 0 |
| ZicoSeq | limma_voom | 29762 | 4 | 1055 | 0 | 1059 | 0 |


## Effect-size agreement

The fairer measure: two methods can rank taxa near-identically yet share no significance calls if one is simply more conservative. Reporting only shared calls would overstate genuine disagreement.

> ⚠ **Not every pair is an independent comparison.** ALDEx2, ZicoSeq and the paired Wilcoxon all report a **CLR difference**. For a balanced matched design the mean of paired differences is *algebraically equal* to the difference of means, so ZicoSeq's and Wilcoxon's effect columns are the **same number** — they correlate at ρ = 1.000 by construction, not by agreement. ZicoSeq's column is not even a ZicoSeq statistic: its own R² is unsigned, so the CLR mean difference is substituted purely to recover direction. Pairs sharing an estimand are flagged below and must never be cited as cross-method validation.

### Independent comparisons — different estimands

| method_a | method_b | spearman | n | pct_same_sign | estimand_a | estimand_b | same_estimand | independent_comparison |
|---|---|---|---|---|---|---|---|---|
| LinDA | ZicoSeq | 0.9773 | 29762 | 94.68 | log2FC relative abundance | CLR difference | False | True |
| ALDEx2 | LinDA | 0.9621 | 41405 | 92.67 | CLR difference | log2FC relative abundance | False | True |
| LinDA | Wilcoxon_CLR_paired | 0.9595 | 46080 | 91.58 | log2FC relative abundance | CLR difference | False | True |
| ZicoSeq | limma_voom | 0.6253 | 29762 | 75 | CLR difference | log2FC logCPM | False | True |
| LinDA | limma_voom | 0.6095 | 46080 | 73.94 | log2FC relative abundance | log2FC logCPM | False | True |
| ALDEx2 | limma_voom | 0.5955 | 41405 | 76.33 | CLR difference | log2FC logCPM | False | True |
| Wilcoxon_CLR_paired | limma_voom | 0.5935 | 46080 | 71.79 | CLR difference | log2FC logCPM | False | True |
| ANCOMBC2 | ZicoSeq | 0.2374 | 21314 | 58.89 | log2FC absolute abundance | CLR difference | False | True |
| ANCOMBC2 | Wilcoxon_CLR_paired | 0.2357 | 21762 | 58.9 | log2FC absolute abundance | CLR difference | False | True |
| ANCOMBC2 | LinDA | 0.2307 | 21762 | 58.85 | log2FC absolute abundance | log2FC relative abundance | False | True |
| ALDEx2 | ANCOMBC2 | 0.1898 | 21762 | 57.12 | CLR difference | log2FC absolute abundance | False | True |
| ANCOMBC2 | limma_voom | 0.1736 | 21762 | 57.34 | log2FC absolute abundance | log2FC logCPM | False | True |


### Same-estimand pairs — agreement is partly definitional

| method_a | method_b | spearman | n | pct_same_sign | estimand_a | estimand_b | same_estimand | independent_comparison |
|---|---|---|---|---|---|---|---|---|
| Wilcoxon_CLR_paired | ZicoSeq | 1 | 29762 | 100 | CLR difference | CLR difference | True | False |
| ALDEx2 | ZicoSeq | 0.9907 | 29762 | 96.2 | CLR difference | CLR difference | True | False |
| ALDEx2 | Wilcoxon_CLR_paired | 0.9871 | 41405 | 95.29 | CLR difference | CLR difference | True | False |


Among the **independent** comparisons, five methods agree closely on direction (91–95% same-sign). **ANCOM-BC2 is the sole outlier at 57–59% — coin-flip.** That is the finding: the method calling 31.5% of taxa significant is also the one whose effect estimates do not reproduce in any other method.


## Consensus

Taxon-analyses called significant by ≥1 method: **7,558**

| n_methods | n_taxon_analyses |
|---|---|
| 3 | 1 |
| 2 | 452 |
| 1 | 7105 |


**Called by ≥2 methods: 453**

| level | prevalence | partition | taxon | n_methods | n_tested | methods | methods_not_tested |
|---|---|---|---|---|---|---|---|
| species | prev10 | any_ige | Prevotella_copri | 3 | 6 | ANCOMBC2,LinDA,limma_voom |  |
| asv | prev01 | any_ige | 05d757dba88ddd1e1e4e9974857a2b09 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 0d89c7f97685c1a34a36f86096c43246 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 155ace0850fc88498e3c881bc5252fcc | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 15ceb00ebfbfe7d69e79f08436f255f9 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 18f557ed1fb6026bb8a6398984d19dd4 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 19ff2a801156e1930f697799d400fa95 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 1f2dda5d2f361fc8d646e55a515ae1e4 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 26639d9f66c6c2106996b1c162605ff8 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 2dbb2410d57cf07b9448d6dc0d5cf236 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 2ee585c245a593473acb21bdf1d7a460 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 317c294f3f66c0443f743165841fce56 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 37020b47521b3745084f867b8c1f7bc3 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 38b7580c13cf081c9f279006c2dc70b9 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 3a893286e7a8c6965a9be85b7ec43c60 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 3fd0c1aebeadfd9c4f8a66a6680b310f | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 409f711b59152d57926cf444c5577087 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 414ea75fb1462d10153c7bed625b4b99 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 433d1758d4b53bad06cf0a8796ce3263 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 44cae38ace30f54ee81cc5c5b2ce47eb | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 463ccc4044935acfa08235ade550557e | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 47d74904a190fe918c6eea86295320d2 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 494e65560668c776b1b9e745b349c50d | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 4a5d2a184f229de74e17565e4c13a7d9 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 4f9106b58ba967a459b11a9d6f9ca4e0 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 53f976e60f93b527ffa62fbe48254d53 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 544bdbfbed9a11965f1b4457b4c00d1c | 2 | 5 | ANCOMBC2,limma_voom | ZicoSeq |
| asv | prev01 | any_ige | 557156f08d4eab8818ea1914d5f83ed0 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 582fae36f33acd1efdcaf7cacf00ef0a | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 60bed339f40d58c6de65c2b894d2132e | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 62c336937d9c98f93521f014e1f21a63 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 64cd0cec44c7015ce37e278c8c8a558c | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 6b6d13f1332283e053c18c3b7774b85c | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 6bea7adf118ca77f31985d2caecd0cff | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 6d88e732578bfda8c3afd4050222d185 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 6e20979e3b8626cbb61b64fac8dc2f0c | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 707940842caa2afe60491008e04a8173 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 747a284c95e1af0ec02d6582f6565774 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 75a7dd04040e23328468b763836841ac | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 77087c966bd1f46d43d2c4913475bc98 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 7842a433c5c0271371dba64a011c2d4b | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 9037a8807f265f6926427153c9f6f6f8 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 9243bc95f3bb4eba8fabc522fa335a7f | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 9591a959c4c592a119485d2f38ae27eb | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 97d5ed404ba83dcbaa164b1b8aa01d57 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 9af2911f6df6d008cf67880633e1d491 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 9b7a8a443b34fb607fae445c0b149466 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 9cb986fb9a5d33f51977148b33c14fcb | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | 9eaff8c3aca7b74962decea682803d9c | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | a1b6935c2608151d83edd47c8a72c1a6 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | a29990fffb08af05632b35e9d3d8fd77 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | a3d6369999378cdeb1c05d757e2ef00d | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | a79603bb90824ed6e253837fe96c6c49 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | b0cce80cd5e1c78e70e026edb63a65f9 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | b24ef5fe8c1ad555a6db6fb7800398d1 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | b74ab6944f5a7e98a99e29bcbea455c0 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | b8349e5459469537ac814c90cb3b99fa | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | b8f980ec6ad0c86efb83bd1f64c3a252 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | ba968750389f02d4d38a04a5bbc736cc | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | bc0730096682c48820678c59700fe66f | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | be1fe9bfe424ad27f568d8e2d7b42c19 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | c18dae5503801a39d3adaf67ecee7a9a | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | c1a89d5c446ea0400780188ff84729e3 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | c2abea83c2d018c7d888a1ecf9f5785e | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | c39cea2193e0d56af9669a1b2d5c21be | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | caf9cf329b5a33fce30625214429681b | 2 | 5 | ANCOMBC2,limma_voom | ZicoSeq |
| asv | prev01 | any_ige | ce0837f50ab0a7abc73bddc82e8b55e6 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | d0377d9209ed8843c1b85df7bc622a2c | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | d331c54836ce2164d5d51154f369d9a4 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | d347cad164508355e12daadf9ccea765 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | d4d386b7d01baf327f6b43aa568cb9e7 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | d5d977abb40d21dd360b1712ae9b839a | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | d6bff94349f93df9c98368bfbc740ed0 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | da08cefbe17c776cf78b0634a8377efa | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | ddf8c140c0dc12be5d69696cdddb91e0 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | df5c1ef25ab497ab767dcbb6be082f83 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | e31c7e269125382d9bd6ee8fcc1fde21 | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | e4cd8e09dbde32eab75a5a428425260d | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | ed7747f5b0284274c95a92cc42e8657d | 2 | 6 | ANCOMBC2,limma_voom |  |
| asv | prev01 | any_ige | f4fd4d5136d00b6ba8b67030e932e2be | 2 | 6 | ANCOMBC2,limma_voom |  |

_Showing 80 of 453 rows._


### By taxonomic rank

| rank | any_method | ge_2_methods | ge_3_methods |
|---|---|---|---|
| family | 221 | 3 | 0 |
| genus | 832 | 11 | 0 |
| species | 719 | 62 | 1 |
| asv | 5786 | 377 | 0 |


## Two cautions that must travel with these numbers

1. **Effect sizes are never averaged across methods.** They are on different
   scales — log2 fold change of relative abundance, natural-log LFC, CLR
   difference, log2 CPM fold change. An average would mean nothing.

2. **ZicoSeq at ASV rank must be read via `q_native`, not `q_within`.** A
   permutation p-value has a hard floor of 1/(nperm+1), so BH over *m* features
   cannot fall below *m*/(nperm+1). At nperm = 4999 that floor is 0.012 for 60
   features but **0.29 for 1,459** — structurally incapable of reaching 0.05.
   `q_native` is ZicoSeq's own permutation FDR, which is not floor-limited the
   same way. The scripts print this floor for every analysis.


---
_Generated 2026-09-26 by `revision/analysis/10_export_markdown.py`._
