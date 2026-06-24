---
name: SPIRE
description: >
  SPIRE (Synovial Pathway Immune-subpopulation RA-resemblance Evaluation) is a
  cross-species single-cell pathway evaluation framework for comparing mouse
  arthritis models against human RA at the pathway-activity level. Use this
  skill whenever the task involves Hallmark AUCell scoring, cross-species
  pathway-pattern comparison, mouse model evaluation against hsRA, Spearman
  correlation across matched synovial cell types, SPIRE-style figures, or
  methods/results writing for a disease-state-only model evaluation workflow.
---

# SPIRE — Synovial Pathway Immune-subpopulation RA-resemblance Evaluation

## Overview

SPIRE compares mouse arthritis models against human RA at the pathway-activity
level using single-cell RNA-seq data. It does not compare absolute AUCell
scores across species or datasets. Instead, it compares relative pathway
activity patterns across matched synovial cell subpopulations.

Core logic: if a mouse model and human RA show similar cell-type-specific
Hallmark pathway patterns, the mouse model is considered to better recapitulate
human RA at the pathway level.

Key design principle: SPIRE uses disease-state samples only. Healthy, naive, and
OA control samples are not required for the cross-species comparison. This
avoids asymmetric comparisons such as human RA versus OA but mouse disease
versus naive.

## Pipeline overview

```text
AUCell scoring
  -> Step 1: prepare inputs and equal-sample AUCell means
  -> Step 2: build model profiles and within-model column-wise z-score
  -> Step 3: hsRA-vs-mouse Spearman correlation matrices
  -> Step 4: pairwise and global permutation tests
  -> Step 5: core figures
  -> Step 6: pathway driver analysis
  -> Step 7: sensitivity analysis
  -> Step 8: final summary and protocol/cell-type comparison plots
```

## Pre-pipeline: AUCell scoring

AUCell scoring is performed before entering the SPIRE comparison workflow.

Requirements:

- QC, clustering, and cell-type annotation are completed independently for each dataset.
- AUCell scores are computed per cell for all 50 MSigDB Hallmark pathways.
- Ortholog mapping is applied before scoring when cross-species gene sets are used.
- Only disease-state samples are passed into the SPIRE comparison.
- No minimum post-AUCell cell-count threshold is applied inside SPIRE.

AUCell is used because it is rank-based within each cell and is therefore more
appropriate for cross-dataset and cross-species pathway activity comparison than
absolute expression module scores.

## Step 1 — Prepare inputs

Input: per-dataset CSV/RDS/Seurat-derived tables containing cell-level AUCell
scores.

Required columns:

- `cell_id`
- `model_id`
- `cell_type`
- `sample_id`, recommended when available
- 50 `HALLMARK_*` AUCell columns

Operations:

1. Standardize metadata columns and canonical cell-type labels.
2. Assign model labels, for example `hsRA`, `STIA`, `GPI`, `CIA`, and `C&AIA`.
3. Retain disease-state samples only.
4. Retain all cell subgroups present after upstream QC and AUCell scoring; do
   not remove subgroups because of low cell number.
5. Compute sample-level AUCell means first:

```text
cell-level AUCell -> sample_id × cell_type mean
```

6. Compute dataset/source-level means using equal sample weight:

```text
sample_id × cell_type means -> source × cell_type mean
```

This avoids giving larger-cell-count samples larger influence.

If `sample_id` is unavailable, the source is treated as `n = 1` and interpreted
accordingly.

Output:

- `harmonized_cell_level_AUCell.csv`
- `sample_level_mean_AUCell.csv`
- `dataset_equal_sample_mean_AUCell.csv`
- `model_equal_source_mean_AUCell.csv`
- QC tables documenting sample and cell-type coverage

## Step 2 — Z-score profiles

This is the most critical step.

Z-score rule: for each `model_id` independently, each Hallmark pathway column is
z-scored across the included cell-type rows.

```r
for (mid in unique(profiles$model_id)) {
  idx <- which(profiles$model_id == mid)
  mat <- as.matrix(profiles[idx, ..pathways])
  z_mat <- scale(mat)  # column-wise: each pathway across cell types
  profiles[idx, (pathways) := as.data.table(z_mat)]
}
```

This answers:

```text
For this model, in which cell types is this pathway relatively high or low?
```

This is not row-wise z-scoring across pathways within each cell type.

Wrong direction:

```r
z_mat <- t(apply(mat, 1, scale))  # do not use for SPIRE
```

Why column-wise z-score is used: it removes pathway-specific baselines shared
across cell types and focuses the comparison on cell-type-specific pathway
preferences.

Verification: after z-scoring, for each `model_id` and each Hallmark pathway,
the mean across available cell types should be approximately 0 and the standard
deviation should be approximately 1, unless too few or invariant cell types are
available.

Missing cell types remain missing and are not imputed.

Outputs:

- `zscore_model_profiles_wide.csv`
- `zscore_model_profiles_long.csv`
- `zscore_verification.csv`, recommended

## Step 3 — Correlation matrices

For each mouse model, compute the Spearman correlation between hsRA and mouse
cell-type profiles.

Matrix definition:

```text
Rows    = hsRA cell types
Columns = mouse model cell types
Values  = Spearman rho between 50-Hallmark z-score vectors
```

The main protocol comparison is restricted to five included dimensions:

```text
Macrophage, Monocyte, DC, T cell, FLS
```

If a model lacks one dimension, that dimension remains missing rather than being
replaced by another cell type.

Self-match rank: for each hsRA cell type, rank all available mouse cell types by
Spearman rho. Rank 1 means the same-name mouse cell type is the top match.

Outputs:

- `all_hsRA_vs_mouse_spearman_long.csv`
- `hsRA_vs_[model]_spearman_matrix.csv`
- `matched_celltype_correlations.csv`
- `self_match_ranks.csv`
- `model_level_summary.csv`

## Step 4 — Permutation tests

Step 4 evaluates statistical support for the correlation structure.

### 4A. Pairwise pathway-label permutation

For each hsRA cell type × mouse cell type pair:

1. Keep the hsRA vector fixed.
2. Shuffle the mouse vector's Hallmark pathway labels.
3. Recompute Spearman rho.
4. Repeat, for example 1000 times.

This tests whether a pairwise pathway-pattern correlation is stronger than
expected from random pathway alignment.

### 4B. Global cell-label permutation

For each mouse model:

1. Compute the observed statistic:

```text
mean matched rho - mean mismatched rho
```

2. Shuffle mouse cell-type labels within the included comparable dimensions.
3. Recompute the statistic.
4. Repeat, for example 10000 times.

This tests whether the same-name cell-type matching structure is stronger than
expected by chance.

Outputs:

- `pairwise_pathway_permutation_results.csv`
- `matched_celltype_significance.csv`
- `global_matched_vs_mismatched_permutation.csv`
- `global_null_distribution.csv`

## Step 5 — Core figures

Recommended core figures:

| Figure | Content | Role |
|---|---|---|
| F1 | 5 × 5 Spearman rho heatmaps per model | Core cross-species evidence |
| F2 | Matched cell-type correlation dot plot with BH labels | Pairwise statistical summary |
| F3 | Matched-minus-mismatched advantage plot | Global structure summary |
| F4 | Self-match rank heatmap | Cell-type identity validation |

These figures should be restricted to the five included dimensions unless a
specific supplementary analysis states otherwise.

## Step 6 — Pathway driver analysis

For each mouse model × matched cell type, classify all 50 Hallmark pathways:

- Concordantly UP: `hsRA_z > 0` and `mouse_z > 0`
- Concordantly DOWN: `hsRA_z < 0` and `mouse_z < 0`
- Discordant: opposite signs

Contribution score:

```text
product_contribution = hsRA_z × mouse_z
```

Interpretation:

- Positive product: directional concordance.
- Negative product: directional discordance.
- Larger magnitude: stronger contribution to similarity or divergence.

Outputs:

- `pathway_driver_all_pairs.csv`
- `top_concordant_up_pathways.csv`
- `top_concordant_down_pathways.csv`
- `top_discordant_pathways.csv`
- `pathway_concordance_summary.csv`
- `model_pathway_summary.csv`

## Step 7 — Sensitivity analysis

Step 7 evaluates whether the model comparison is robust to pathway and cell-type
composition.

Analyses:

| Analysis | Pathways or dimensions | Purpose |
|---|---|---|
| `full_50_pathways` | All 50 Hallmark pathways | Baseline |
| `without_broad_pathways` | Remove broad/cell-cycle/metabolic pathways | Test whether broad biological state drives the result |
| `leave_one_celltype_out` | Drop one of the five included cell types at a time | Test whether results depend on a single cell type |

The immune-only pathway subset is not part of the current SPIRE default
workflow, because it can introduce subjective pathway-classification decisions
and reduce pathway coverage.

Outputs:

- `sensitivity_model_summary.csv`
- `sensitivity_leave_one_celltype_out_summary.csv`
- `model_rank_by_each_sensitivity_analysis.csv`
- `model_rank_stability.csv`
- `target_dimension_coverage.csv`

## Step 8 — Final summary and comparison plots

Step 8 integrates the results without collapsing them into a single composite
score. SPIRE should not announce a single "best model" based on an arbitrary
weighted score.

Recommended outputs:

- `final_model_summary.csv`
- `final_matched_celltype_summary.csv`
- `final_pathway_driver_summary.csv`
- `final_sensitivity_summary.csv`
- `figure_caption_statistics.csv`

### F5 — Protocol-dimension dot plot

F5 summarizes the five protocol dimensions:

1. Mean matched Spearman rho.
2. Rank 1 self-match fraction.
3. Global structure support, `-log10(p_value_greater)`.
4. Sensitivity retention ratio.
5. Discordant pathway burden, shown as a reverse-oriented metric when scaled.

The plot should use scaled positions only for visualization and should label raw
values directly. This avoids misinterpreting a scaled 0 as an original value of
0.

### F6 — Cell-type matched rho radar

F6 compares the matched Spearman rho for the five included cell types:

```text
Macrophage, Monocyte, DC, T cell, FLS
```

Each model is shown as one polygon. The axis values are original matched
Spearman rho values. Missing dimensions, such as STIA missing Monocyte, should
remain missing and should not be imputed.

## Composite figure layout

Recommended manuscript figure structure:

```text
Panel A: UMAP / cell-type context for each dataset or model
Panel B: Hallmark z-score overview heatmap
Panel C: hsRA-vs-mouse 5 × 5 correlation heatmaps
Panel D: Self-match rank heatmap
Panel E: Cell-type matched rho radar
Panel F: Matched cell-type correlation dot plot
Panel G: Matched-minus-mismatched global structure plot
Supplementary: protocol-dimension dot plot, sensitivity tables, pathway drivers
```

## Key methodological decisions

1. AUCell is used because it is rank-based within each cell.
2. All 50 Hallmark pathways are used by default to avoid data-driven pathway selection.
3. Z-score is column-wise within each model, across cell types for each pathway.
4. Disease-state-only comparison avoids asymmetric control definitions across species.
5. No data-driven pathway filtering is performed before correlation.
6. Partial cell-type coverage is allowed, but missing dimensions are documented and not imputed.
7. Spearman correlation is used because it is rank-based and robust to scale differences.
8. No single composite score is used to declare one universal best model.

## Data requirements

| Field | Required | Notes |
|---|---:|---|
| `model_id` | Yes | Model or dataset label |
| `cell_type` | Yes | Canonical cell-type name |
| `sample_id` | Recommended | Used for equal-sample-weight means; if absent, source is treated as n = 1 |
| 50 `HALLMARK_*` columns | Yes | AUCell scores |

Minimum requirement for z-score: at least two comparable cell types per model.
No minimum post-AUCell cell-count threshold is applied inside SPIRE.

## Scripts

Run scripts in order:

```text
01_step1_prepare_inputs.R
02_step2_make_profiles.R
03_step3_correlation_matrices.R
04_step4_permutation_tests.R
05_step5_make_figures.R
06_step6_pathway_drivers.R
07_step7_sensitivity_analysis.R
08_step8_final_summary_and_radar.R
```

Dependencies:

```r
data.table
ggplot2
```

