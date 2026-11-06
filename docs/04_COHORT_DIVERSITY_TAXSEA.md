# Cohort, diversity, networks and classification


## Table 1 — full cohort (n = 508), formally tested

Tests the submitted Table 1 never reported.

| table | variable | column | type | is_collider | n_sensitized | n_non_sensitized | sensitized | non_sensitized | test | statistic | p_value | effect_name | effect_size | effect_ci_low | effect_ci_high | note | q_BH |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Table1_full508 | Age (years) | u3talter | cont | False | 275 | 233 | 51.09 (6.66) | 51.97 (7.01) | Wilcoxon rank-sum | 2.959e+04 | 0.1373 | Cliff's delta | -0.07641 |  |  | Shapiro p=6.11e-10 | 0.3295 |
| Table1_full508 | BMI (kg/m2) | u3tbmi | cont | False | 275 | 233 | 26.67 (5.17) | 26.85 (4.38) | Wilcoxon rank-sum | 3.002e+04 | 0.2201 | Cliff's delta | -0.06311 |  |  | Shapiro p=2.63e-17 | 0.3302 |
| Table1_full508 | Sex | u3csex | cat | False | 275 | 233 | 116 (42.2%) / 159 (57.8%) | 69 (29.6%) / 164 (70.4%) | Chi-squared | 8.604 | 0.003354 | Cramer's V | 0.1301 | -0.2397 | -0.007238 | levels: Male,Female; risk diff=-0.126 | 0.01341 |
| Table1_full508 | Education | u3tschul | cat | False | 275 | 233 | 86 (31.3%) / 100 (36.4%) / 89 (32.4%) | 90 (38.6%) / 71 (30.5%) / 72 (30.9%) | Chi-squared | 3.355 | 0.1869 | Cramer's V | 0.08126 |  |  | levels: 1,2,3 | 0.3302 |
| Table1_full508 | Smoking | u3tcigreg_sf | cat | False | 275 | 233 | 30 (10.9%) / 115 (41.8%) / 130 (47.3%) | 38 (16.3%) / 92 (39.5%) / 103 (44.2%) | Chi-squared | 3.175 | 0.2045 | Cramer's V | 0.07905 |  |  | levels: Regular,Ex-smoker,Never | 0.3302 |
| Table1_full508 | Pack-years | u3tpyrs | cont | False | 275 | 233 | 7.03 (12.27) | 8.39 (12.64) | Wilcoxon rank-sum | 3.029e+04 | 0.2651 | Cliff's delta | -0.05453 |  |  | Shapiro p=6.99e-30 | 0.3534 |
| Table1_full508 | Alcohol (g/day) | u3talkkon | cont | False | 275 | 233 | 9.34 (9.95) | 8.90 (10.01) | Wilcoxon rank-sum | 3.299e+04 | 0.5601 | Cliff's delta | 0.02975 |  |  | Shapiro p=6.61e-22 | 0.6138 |
| Table1_full508 | Physical activity | u3tphact | cat | False | 275 | 233 | 179 (65.1%) / 44 (16.0%) / 52 (18.9%) | 154 (66.1%) / 30 (12.9%) / 49 (21.0%) | Chi-squared | 1.15 | 0.5627 | Cramer's V | 0.04758 |  |  | levels: Regular,Irregular,None | 0.6138 |
| Table1_full508 | PPI use | u3tmppi | cat | False | 275 | 233 | 256 (93.1%) / 19 (6.9%) | 217 (93.1%) / 16 (6.9%) | Chi-squared | 0.0003491 | 0.9851 | Cramer's V | 0.000829 | -0.06395 | 0.06267 | levels: No,Yes; risk diff=0.000 | 0.9851 |
| Table1_full508 | Asthma | u3da13aa | cat | True | 275 | 233 | 236 (85.8%) / 39 (14.2%) | 227 (97.4%) / 6 (2.6%) | Chi-squared | 21.05 | 4.483e-06 | Cramer's V | 0.2035 | 0.05048 | 0.1761 | levels: No,Yes; risk diff=0.116 | 2.69e-05 |
| Table1_full508 | Dermatitis | u3da4aa | cat | True | 275 | 233 | 246 (89.5%) / 29 (10.5%) | 221 (94.8%) / 12 (5.2%) | Chi-squared | 4.948 | 0.02611 | Cramer's V | 0.0987 | -0.01342 | 0.1176 | levels: No,Yes; risk diff=0.054 | 0.07834 |
| Table1_full508 | Hay fever | u3da8aa | cat | True | 275 | 233 | 146 (53.1%) / 129 (46.9%) | 225 (96.6%) / 8 (3.4%) | Chi-squared | 121 | 3.734e-28 | Cramer's V | 0.4881 | 0.3447 | 0.5106 | levels: No,Yes; risk diff=0.435 | 4.481e-27 |


## Table 1 — matched subset QC

Confirms the 1:1 matching balanced the covariates it was built on.

| table | variable | column | type | is_collider | n_sensitized | n_non_sensitized | sensitized | non_sensitized | test | statistic | p_value | effect_name | effect_size | effect_ci_low | effect_ci_high | note | q_BH |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Table1_matched442 | Age (years) | u3talter | cont | False | 221 | 221 | 51.53 (6.83) | 51.70 (6.93) | Wilcoxon rank-sum | 2.403e+04 | 0.7687 | Cliff's delta | -0.01617 |  |  | Shapiro p=4.99e-09 | 0.8478 |
| Table1_matched442 | BMI (kg/m2) | u3tbmi | cont | False | 221 | 221 | 26.15 (4.11) | 26.48 (4.16) | Wilcoxon rank-sum | 2.318e+04 | 0.3567 | Cliff's delta | -0.0507 |  |  | Shapiro p=3.26e-11 | 0.6115 |
| Table1_matched442 | Sex | u3csex | cat | False | 221 | 221 | 68 (30.8%) / 153 (69.2%) | 68 (30.8%) / 153 (69.2%) | Chi-squared | 0 | 1 | Cramer's V | 0 | -0.1208 | 0.1208 | levels: Male,Female; risk diff=0.000 | 1 |
| Table1_matched442 | Education | u3tschul | cat | False | 221 | 221 | 63 (28.5%) / 86 (38.9%) / 72 (32.6%) | 84 (38.0%) / 67 (30.3%) / 70 (31.7%) | Chi-squared | 5.388 | 0.06762 | Cramer's V | 0.1104 |  |  | levels: 1,2,3 | 0.2029 |
| Table1_matched442 | Smoking | u3tcigreg_sf | cat | False | 221 | 221 | 22 (10.0%) / 90 (40.7%) / 109 (49.3%) | 35 (15.8%) / 86 (38.9%) / 100 (45.2%) | Chi-squared | 3.443 | 0.1788 | Cramer's V | 0.08826 |  |  | levels: Regular,Ex-smoker,Never | 0.429 |
| Table1_matched442 | Pack-years | u3tpyrs | cont | False | 221 | 221 | 6.44 (10.95) | 8.30 (12.82) | Wilcoxon rank-sum | 2.287e+04 | 0.2222 | Cliff's delta | -0.06349 |  |  | Shapiro p=5.58e-28 | 0.4444 |
| Table1_matched442 | Alcohol (g/day) | u3talkkon | cont | False | 221 | 221 | 8.81 (9.77) | 9.08 (10.04) | Wilcoxon rank-sum | 2.404e+04 | 0.7772 | Cliff's delta | -0.01544 |  |  | Shapiro p=1e-20 | 0.8478 |
| Table1_matched442 | Physical activity | u3tphact | cat | False | 221 | 221 | 145 (65.6%) / 36 (16.3%) / 40 (18.1%) | 148 (67.0%) / 28 (12.7%) / 45 (20.4%) | Chi-squared | 1.325 | 0.5156 | Cramer's V | 0.05475 |  |  | levels: Regular,Irregular,None | 0.6875 |
| Table1_matched442 | PPI use | u3tmppi | cat | False | 221 | 221 | 203 (91.9%) / 18 (8.1%) | 207 (93.7%) / 14 (6.3%) | Chi-squared | 0.539 | 0.4628 | Cramer's V | 0.03492 | -0.05137 | 0.08695 | levels: No,Yes; risk diff=0.018 | 0.6875 |
| Table1_matched442 | Asthma | u3da13aa | cat | True | 221 | 221 | 189 (85.5%) / 32 (14.5%) | 216 (97.7%) / 5 (2.3%) | Chi-squared | 21.5 | 3.533e-06 | Cramer's V | 0.2206 | 0.05261 | 0.1876 | levels: No,Yes; risk diff=0.122 | 2.12e-05 |
| Table1_matched442 | Dermatitis | u3da4aa | cat | True | 221 | 221 | 195 (88.2%) / 26 (11.8%) | 209 (94.6%) / 12 (5.4%) | Chi-squared | 5.643 | 0.01752 | Cramer's V | 0.113 | -0.01093 | 0.1355 | levels: No,Yes; risk diff=0.063 | 0.0701 |
| Table1_matched442 | Hay fever | u3da8aa | cat | True | 221 | 221 | 114 (51.6%) / 107 (48.4%) | 213 (96.4%) / 8 (3.6%) | Chi-squared | 115.2 | 7.12e-27 | Cramer's V | 0.5105 | 0.3493 | 0.5313 | levels: No,Yes; risk diff=0.448 | 8.544e-26 |


## Table 2 — IgE cluster descriptives

| cluster | n | age_mean | age_sd | bmi_mean | bmi_sd | pct_female |
|---|---|---|---|---|---|---|
| A | 12 | 51.17 | 6.132 | 25.85 | 2.765 | 50 |
| B | 24 | 53.88 | 6.688 | 26.86 | 5.812 | 62.5 |
| C | 26 | 51.81 | 7.531 | 25.81 | 4.942 | 76.92 |
| D | 45 | 48.8 | 5.98 | 25.76 | 3.991 | 51.11 |
| E | 22 | 50.14 | 7.213 | 26.08 | 3.201 | 40.91 |
| F | 10 | 45.2 | 4.872 | 25.94 | 2.208 | 50 |
| G | 22 | 48.36 | 5.141 | 26.17 | 4.385 | 54.55 |
| H | 40 | 51.52 | 7.197 | 28.2 | 5.198 | 50 |
| I | 74 | 52.97 | 5.979 | 27.2 | 6.655 | 66.22 |


## Table 2 — IgE clusters, formally tested

Clusters F (n = 10) and A (n = 12) are underpowered; read with that in mind.

| variable | column | type | test | statistic | df | p_value | effect_name | effect_size | q_BH | underpowered_note |
|---|---|---|---|---|---|---|---|---|---|---|
| Age (years) | u3talter | cont | Kruskal-Wallis | 29.48 | 8 | 0.0002614 | eta-squared | 0.08074 | 0.001569 | clusters F (n=10) and A (n=12) are underpowered |
| BMI (kg/m2) | u3tbmi | cont | Kruskal-Wallis | 7.187 | 8 | 0.5166 | eta-squared | -0.003057 | 0.6092 | clusters F (n=10) and A (n=12) are underpowered |
| Sex | u3csex | cat | Fisher's exact (simulated) |  |  | 0.1717 | Cramer's V | 0.2028 | 0.4372 | clusters F (n=10) and A (n=12) are underpowered |
| Education | u3tschul | cat | Fisher's exact (simulated) |  |  | 0.0033 | Cramer's V | 0.2547 | 0.0132 | clusters F (n=10) and A (n=12) are underpowered |
| Smoking | u3tcigreg_sf | cat | Fisher's exact (simulated) |  |  | 0.3522 | Cramer's V | 0.1767 | 0.4696 | clusters F (n=10) and A (n=12) are underpowered |
| Pack-years | u3tpyrs | cont | Kruskal-Wallis | 11.19 | 8 | 0.1911 | eta-squared | 0.012 | 0.4372 | clusters F (n=10) and A (n=12) are underpowered |
| Alcohol (g/day) | u3talkkon | cont | Kruskal-Wallis | 10.14 | 8 | 0.255 | eta-squared | 0.008061 | 0.4372 | clusters F (n=10) and A (n=12) are underpowered |
| Physical activity | u3tphact | cat | Fisher's exact (simulated) |  |  | 0.9459 | Cramer's V | 0.1243 | 0.9459 | clusters F (n=10) and A (n=12) are underpowered |
| PPI use | u3tmppi | cat | Fisher's exact (simulated) |  |  | 0.339 | Cramer's V | 0.1844 | 0.4696 | clusters F (n=10) and A (n=12) are underpowered |
| Asthma | u3da13aa | cat | Fisher's exact (simulated) |  |  | 0.2422 | Cramer's V | 0.1844 | 0.4372 | clusters F (n=10) and A (n=12) are underpowered |
| Dermatitis | u3da4aa | cat | Fisher's exact (simulated) |  |  | 0.5584 | Cramer's V | 0.1492 | 0.6092 | clusters F (n=10) and A (n=12) are underpowered |
| Hay fever | u3da8aa | cat | Fisher's exact (simulated) |  |  | 1e-05 | Cramer's V | 0.5798 | 0.00012 | clusters F (n=10) and A (n=12) are underpowered |


## TaxSEA, recomputed

The re-run that settled the impossible FDR < p rows in the submitted supplement.

| sub_database | taxonSetName | median_rank_of_set_members | Test_statistic | PValue | FDR | prevalence | rank_by | n_input_species | n_mapped |
|---|---|---|---|---|---|---|---|---|---|
| All_databases | BacDive_Utilizes_cellobiose | -0.323 | 0.4094 | 0.1275 | 1 | prev10 | stat | 44 | 38 |
| All_databases | GutMGene_producers_of_Equol | -0.2651 | 0.4737 | 0.1419 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_glucose | -0.5343 | 0.4079 | 0.1666 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_rhamnose | -0.5343 | 0.4079 | 0.1666 | 1 | prev10 | stat | 44 | 38 |
| All_databases | Valles-Colomer2019_MGB040.Menaquinone.synthesis.(vitamin.K2).I | -0.4291 | 0.3772 | 0.1854 | 1 | prev10 | stat | 44 | 38 |
| All_databases | GutMGene_producers_of_Acetate | 0.8092 | 0.4737 | 0.1935 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_maltose | -0.2105 | 0.2982 | 0.2233 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_arabinose | -0.323 | 0.3722 | 0.2945 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Enzyme_glutamyl-glutamate arylamidase | -0.1273 | 0.3211 | 0.2979 | 1 | prev10 | stat | 44 | 38 |
| All_databases | GMRepoV2_Decreased_in_Diabetes_Mellitus,_Type_2 | 0.4493 | 0.4211 | 0.3074 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_lactose | -0.182 | 0.2664 | 0.3123 | 1 | prev10 | stat | 44 | 38 |
| All_databases | GMRepoV2_Increased_in_COVID-19 | 0.1728 | 0.4211 | 0.3151 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Enzyme_beta-glucosidase | -0.182 | 0.2744 | 0.3247 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_salicin | -0.323 | 0.3459 | 0.3797 | 1 | prev10 | stat | 44 | 38 |
| All_databases | GMRepoV2_Decreased_in_Colitis,_Ulcerative | -0.0455 | 0.3421 | 0.3894 | 1 | prev10 | stat | 44 | 38 |
| All_databases | GMRepoV2_Increased_in_Adenoma | -0.7456 | 0.3895 | 0.4014 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_D-xylose | -0.1535 | 0.2679 | 0.4592 | 1 | prev10 | stat | 44 | 38 |
| All_databases | GMRepoV2_Decreased_in_Adenoma | 0.1269 | 0.3684 | 0.4706 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_trehalose | 0.6245 | 0.3421 | 0.4727 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_D-mannose | -0.1535 | 0.2316 | 0.5099 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_raffinose | -0.1535 | 0.2198 | 0.5099 | 1 | prev10 | stat | 44 | 38 |
| All_databases | GMRepoV2_Decreased_in_Cardiovascular_Diseases | 0.03581 | 0.2895 | 0.5278 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_mannose | -0.2383 | 0.2895 | 0.5278 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_xylose | -0.4496 | 0.3246 | 0.5368 | 1 | prev10 | stat | 44 | 38 |
| All_databases | MiMeDB_producers_of_3IAA | 0.3177 | 0.307 | 0.6019 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_sucrose | -0.1273 | 0.2039 | 0.6051 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Enzyme_acid phosphatase | -0.1282 | 0.2368 | 0.6738 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_esculin | -0.1273 | 0.203 | 0.6755 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_L-arabinose | -0.1012 | 0.2632 | 0.6872 | 1 | prev10 | stat | 44 | 38 |
| All_databases | GMRepoV2_Decreased_in_Obesity | -0.07333 | 0.1992 | 0.7104 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Enzyme_tryptophan deaminase | 0.1245 | 0.2434 | 0.7145 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_tryptophan | 0.1245 | 0.2434 | 0.7145 | 1 | prev10 | stat | 44 | 38 |
| All_databases | GMRepoV2_Decreased_in_Scleroderma,_Systemic | -0.06783 | 0.2061 | 0.7263 | 1 | prev10 | stat | 44 | 38 |
| All_databases | GMRepoV2_Increased_in_Spondylitis,_Ankylosing | 0.8092 | 0.2556 | 0.7378 | 1 | prev10 | stat | 44 | 38 |
| All_databases | MiMeDB_producers_of_Dimethylamine_ | 0.4559 | 0.2519 | 0.7503 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Enzyme_esterase (C 4) | 0.8092 | 0.2842 | 0.7574 | 1 | prev10 | stat | 44 | 38 |
| All_databases | MiMeDB_producers_of_GABA | 0.4559 | 0.1472 | 0.7631 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_D-glucose | -0.1012 | 0.1943 | 0.7633 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Enzyme_naphthol-AS-BI-phosphohydrolase | -0.1282 | 0.2368 | 0.7703 | 1 | prev10 | stat | 44 | 38 |
| All_databases | GMRepoV2_Decreased_in_COVID-19 | -0.1008 | 0.2406 | 0.7923 | 1 | prev10 | stat | 44 | 38 |
| All_databases | GMRepoV2_Increased_in_Autism_Spectrum_Disorder | -0.03914 | 0.2544 | 0.8049 | 1 | prev10 | stat | 44 | 38 |
| All_databases | MiMeDB_producers_of_TMAO | 0.1728 | 0.2737 | 0.808 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Gram_negative | -0.1012 | 0.161 | 0.8089 | 1 | prev10 | stat | 44 | 38 |
| All_databases | GMRepoV2_Decreased_in_Spondylitis,_Ankylosing | 0.09766 | 0.1842 | 0.8181 | 1 | prev10 | stat | 44 | 38 |
| All_databases | Valles-Colomer2019_MGB035.Isovaleric.acid.synthesis.II.(KADC.pathway) | 0.7859 | 0.2293 | 0.8305 | 1 | prev10 | stat | 44 | 38 |
| All_databases | GMRepoV2_Decreased_in_Tuberculosis | 0.1269 | 0.1866 | 0.8527 | 1 | prev10 | stat | 44 | 38 |
| All_databases | MiMeDB_producers_of_Vitamin_K | -0.1012 | 0.1579 | 0.8555 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Enzyme_leucine arylamidase | 0.3818 | 0.2105 | 0.8611 | 1 | prev10 | stat | 44 | 38 |
| All_databases | Valles-Colomer2019_MGB037.Inositol.synthesis | -0.1012 | 0.17 | 0.8625 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Enzyme_alpha-arabinosidase | -0.1012 | 0.1988 | 0.8693 | 1 | prev10 | stat | 44 | 38 |
| All_databases | Mucin_utilisers | 0.06364 | 0.2368 | 0.8753 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Enzyme_glutamate decarboxylase | -0.1012 | 0.177 | 0.8833 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Utilizes_L-glutamate | -0.1012 | 0.177 | 0.8833 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_rod-shaped | -0.1008 | 0.1579 | 0.8863 | 1 | prev10 | stat | 44 | 38 |
| All_databases | MiMeDB_producers_of_Isovaleric_acid | 0.03581 | 0.1546 | 0.8871 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Enzyme_alanine arylamidase | 0.03581 | 0.1316 | 0.9021 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Enzyme_leucyl glycin arylamidase | 0.03581 | 0.1316 | 0.9021 | 1 | prev10 | stat | 44 | 38 |
| All_databases | MiMeDB_producers_of_Vitamin_A | -0.06783 | 0.2193 | 0.9042 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Enzyme_beta-galactosidase | -0.07317 | 0.1206 | 0.9237 | 1 | prev10 | stat | 44 | 38 |
| All_databases | BacDive_Enzyme_N-acetyl-beta-glucosaminidase | -0.07333 | 0.1316 | 0.9278 | 1 | prev10 | stat | 44 | 38 |

_Showing 60 of 1,308 rows._


---
_Generated 2026-09-26 by `revision/analysis/10_export_markdown.py`._
