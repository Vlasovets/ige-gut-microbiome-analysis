# Provenance of `utils/`

Vendored on 2026-08-05 from:

    /home/itg/oleg.vlasovets/slr_example/Causal_Sparse_Low_Rank_Microbiome_Tutorial/utils/
    source commit: dc2b583d8f1ea19a1bfe19894293bd0ba9e7cda9

## Why this directory exists

The submitted version of this repository imported `utils` but did not contain it.
That made every notebook non-runnable and, more importantly, put two
methodologically load-bearing implementations outside the published record:

- `helper.transform_features` — the CLR and mCLR transforms
- `pair_matching.discrepancyMatrix` / `construct_network` — the 1:1 matching
  on sex, age (±5 yr) and BMI (±4) that the matched analyses depend on

A reviewer could not have checked either. They are now part of the repository.

## Files

| File | Provides |
|---|---|
| `helper.py` | `transform_features` (CLR/mCLR), `generate_taxa_dict`, `calculate_log_ratios`, `perform_bh_correction_and_filter`, `calculate_obs_stat_and_pvalues` |
| `pair_matching.py` | `discrepancyMatrix`, `construct_network`, `process_matched_pairs` |
| `preprocessing.py` | `filter_and_process_asv_table` and the feature-table filters |
| `plotting.py`, `heatmap.py` | figure helpers |
| `solver.py`, `stability_selection.py` | ADMM solver and StARS subsampling (used by the SLR analyses) |
| `simulation.py` | simulation helpers (not used by this manuscript) |

## Known issue, carried over deliberately

`construct_network` loops `range(0, N_treated + 1)`, which lets the first control
node be emitted as a pair when `i == N_treated`. The pipeline compensates by
dropping pairs whose treated index is out of bounds
(`00a_match_pairs.py`, "valid_keys"). The upstream code is vendored **unmodified**
so it matches what produced the published results; do not silently fix it here.
