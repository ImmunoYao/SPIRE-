---
name: SPIRE
version: 2.0 (code-verified)
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

> **Document status.** This version has been reconciled line-by-line against the
> eight production R scripts (`01_step1` … `08_step8`) and against the C&AIA
> manuscript. Every parameter below reflects what the code actually does, not
> what the original protocol intended. Divergences between code, documentation
> and manuscript are catalogued in **§ Verification findings**.

---

## Overview

SPIRE compares mouse arthritis models against human RA at the pathway-activity
level using single-cell RNA-seq data. It does not compare absolute AUCell
scores across species or datasets. Instead, it compares relative pathway
activity patterns across matched synovial cell subpopulations.

**Core logic.** If a mouse model and human RA show similar cell-type-specific
Hallmark pathway patterns, the mouse model is considered to better recapitulate
human RA at the pathway level.

**Key design principle.** SPIRE uses disease-state samples only. Healthy, naive
and OA control samples are not required for the cross-species comparison. This
avoids asymmetric comparisons such as human RA versus OA against mouse disease
versus naive. The cost is that SPIRE cannot assess whether mouse and human
disease deviate from their respective healthy baselines in the same direction.

**What SPIRE produces.** A cell-type × pathway concordance/discordance
decomposition. Model-level summary metrics are a compression of that
decomposition, not a substitute for it. SPIRE does not emit a single composite
score and does not declare a universally best model.

---

## Model label mapping (must be resolved before release)

The production scripts use the internal label `GPI` for the
glucose-6-phosphate-isomerase-induced arthritis dataset (GSE184609). The
manuscript uses `GIA` throughout. The label appears in:

| Location | Occurrence |
|---|---|
| `01_step1_prepare_inputs.R` | `model_id = "GPI"`; `metadata_assumptions` table |
| `05_step5_make_figures.R` | `PREFERRED_MODEL_ORDER` |
| `06_step6_pathway_drivers.R` | `MODEL_ORDER` |
| `07_step7_sensitivity_analysis.R` | `MODEL_ORDER` |
| `08_step8_final_summary_and_radar.R` | `MODEL_ORDER`, `MODEL_COLORS` key |

Renaming must be applied in all five files simultaneously. `MODEL_COLORS` is
keyed by name, so a partial rename silently produces uncoloured or dropped
polygons rather than an error.

---

## Pipeline overview

```text
AUCell scoring (external to SPIRE)
  -> Step 1: harmonize inputs; three-level equal-weight AUCell means
  -> Step 2: within-model column-wise z-score across cell types   [CRITICAL]
  -> Step 3: hsRA-vs-mouse Spearman correlation matrices, full cell-type space
  -> Step 4: pairwise pathway-label and global cell-label permutation
  -> Step 5: Figure X panels, restricted to the five target cell types
  -> Step 6: pathway driver decomposition, five target cell types
  -> Step 7: sensitivity analyses, five target cell types
  -> Step 8: final summaries and Figure Y panels
```

`04_step4_permutation_tests.R` ends with an explicit `source()` call to
`05_step5_make_figures.R`. Steps 4 and 5 therefore execute as one unit. Remove
that line if the steps are to be run independently.

---

## Analytical scope by step

This table is the single most important reference in this document. The scripts
do **not** apply a uniform cell-type scope, and several reported metrics are
computed on different cell-type spaces.

| Step | Cell types entering the computation | Consequence |
|---|---|---|
| 1 | All annotated cell types | No filtering by cell count or target list |
| 2 (z-score) | **All annotated cell types per model** | The z-score background is the full annotated compartment set, not the five target types |
| 3 (correlations) | All × all | Full rectangular matrix exported |
| 3 (self-match rank) | All available mouse cell types | Superseded by Steps 5/8 for reporting |
| 4A (pairwise permutation) | All pairs in the full matrix | BH within model computed over all pairs |
| 4B (global permutation) | **All mouse cell types** | Different scope from every other reported metric |
| 5 (figures) | Five target types only | Self-match rank and BH q recomputed on this subset |
| 6 (pathway drivers) | Five target types only | z-scores retain the full-set background |
| 7 (sensitivity) | Five target types only | z-scores are **not** recomputed |
| 8 (final metrics) | Five target types, **except** global p from Step 4 | Mixed scope within one figure |

Target cell types, defined identically in Steps 5–8:

```r
TARGET_CELLTYPES <- c("Macrophage", "Monocyte", "DC", "T cell", "FLS")
```

---

## Pre-pipeline: AUCell scoring

Performed before entering SPIRE.

- QC, clustering and cell-type annotation are completed independently for each
  dataset.
- AUCell scores are computed per cell for all 50 MSigDB Hallmark pathways.
- Ortholog mapping is applied before scoring when cross-species gene sets are
  used.
- Only disease-state samples are passed into the SPIRE comparison.
- No minimum post-AUCell cell-count threshold is applied inside SPIRE. Cell
  filtering is entirely upstream.

AUCell is used because it is rank-based within each cell and is therefore more
appropriate for cross-dataset and cross-species comparison than absolute
expression module scores such as AddModuleScore or GSVA.

---

## Step 1 — Prepare inputs

**Script:** `01_step1_prepare_inputs.R`

### Input contract

Per-dataset CSV files containing cell-level AUCell scores.

| Column | Required | Notes |
|---|---|---|
| `cell_id` | Yes | Rewritten as `source_id::cell_id` on output |
| cell-type column | Yes | Column name varies per dataset and is passed explicitly |
| sample column | Yes, or derivable | Either a column, or parsed from `cell_id` via `sample_fun` |
| `HALLMARK_*` | Exactly 50 | The script errors if the count is not exactly 50 |

### Dataset registry as implemented

| model_id | source_id | Accession | Sample filter | Nominal n |
|---|---|---|---|---:|
| `hsRA` | `GSE109449_FLS` | GSE109449 | `cell_id` prefix `RA8`/`RA9` | 2 |
| `hsRA` | `GAS003842_infiltrating` | GAS003842 | `sample ∈ RA1…RA6` | 6 |
| `STIA` | `GSE129087_FLS` | GSE129087 | `sample ∈ inflamed_R1…R3` | 3 |
| `STIA` | `GSE254560_infiltrating` | GSE254560 | `sample == "PBS"` | 1 |
| `GPI` → `GIA` | `GSE184609` | GSE184609 | `sample ∈ D6, D14, D25` | 3 |
| `CIA` | `GSE192504` | GSE192504 | `sample == "CIA"` | 1 |
| `C&AIA` | `C&AIA` | in-house | `sample_id ∈ B1, B2, B3` | 3 |

For GIA, the three retained samples are disease time points (D6/D14/D25) rather
than biological replicates of a single stage. They are averaged with equal
weight, so the GIA profile is a time-point-averaged disease state.

### Cell-type harmonization

A single `common_map` is applied to all datasets:

```r
Fibroblast, Fibroblast0, Fibroblast1, Fibroblast2  -> FLS
B                                                   -> B cell
Macrophage, Monocyte, T cell, DC, Neutrophil,
NK, Endothelial, Pericyte, Osteoblast, Chondrocyte,
Muscle, RBC, Schwann                                -> identity
```

Labels absent from the map are retained verbatim. This is permissive by design;
a typo in an upstream annotation will propagate as a new cell type rather than
raising an error. Inspect `model_celltype_availability.csv` before Step 2.

### Aggregation — three levels, equal weight at each

```text
cell-level AUCell
  -> mean by (model_id, source_id, sample_id, cell_type)   [sample level]
  -> mean by (model_id, source_id, cell_type)              [equal sample weight]
  -> mean by (model_id, cell_type)                         [equal source weight]
```

Equal weighting at each level prevents samples or sources with more cells from
dominating. For `hsRA` and `STIA`, whose FLS and infiltrating compartments come
from different accessions, the final level averages across sources; because the
two sources contribute disjoint cell types, no actual cross-source averaging
occurs for any single cell type.

### Outputs

```text
harmonized_cell_level_AUCell.csv / .rds
sample_celltype_counts.csv
source_celltype_qc.csv
model_celltype_availability.csv
sample_level_mean_AUCell.csv
dataset_equal_sample_mean_AUCell.csv
model_equal_source_mean_AUCell.csv     <- input to Step 2
metadata_assumptions.csv
```

---

## Step 2 — Z-score profiles

**Script:** `02_step2_make_profiles.R`. This is the step where an error changes
every downstream result.

### Rule

For each `model_id` independently, each Hallmark pathway **column** is z-scored
across the **cell-type rows** available for that model.

```r
for (mid in unique(profiles$model_id)) {
  idx <- which(profiles$model_id == mid)
  mat <- as.matrix(profiles[idx, ..pathways])
  z_mat <- scale(mat)   # column-wise: each pathway across cell types
  profiles[idx, (pathways) := as.data.table(z_mat)]
}
```

This answers: *for this model, in which cell types is this pathway relatively
high or low?*

**Wrong direction — do not use:**

```r
z_mat <- t(apply(mat, 1, scale))   # row-wise across pathways within a cell type
```

Row-wise z-scoring does not remove the pathway baseline shared across cell
types and inflates all pairwise Spearman correlations into roughly 0.7–0.95,
including mismatched pairs.

### Background definition

The z-score background is **all cell types annotated for that model**, taken
from `model_equal_source_mean_AUCell.csv`. It is *not* restricted to the five
target cell types. A model annotated with Neutrophil, B cell, NK, Endothelial
and Pericyte in addition to the five target types is normalized against all of
them.

Consequences that must be stated in any Methods section:

1. The reported correlations describe relative pathway preference within each
   model's own annotated compartment set, and those sets differ between models.
2. Statements of the form "z-scoring used the five matched cell types as the
   reference background" are **incorrect** for this implementation.

### Edge-case handling

| Condition | Behaviour |
|---|---|
| Model has < 2 cell types | All pathways set to `NA`, warning emitted |
| Pathway has zero variance across cell types | Set to `NA`, message listing the pathways |
| Missing cell types | Remain missing; never imputed |

### Verification

`zscore_verification.csv` reports, per model, the mean of per-pathway column
means and the mean of per-pathway column SDs. Acceptance thresholds hard-coded
in the script:

```r
abs(mean_of_col_means) <= 0.01
abs(mean_of_col_sds - 1) <= 0.05
```

A violation raises a warning, not an error. Check this file before proceeding.

### Outputs

```text
raw_model_profiles_wide.csv / _long.csv
zscore_model_profiles_wide.csv / _long.csv     <- input to Steps 3, 4, 6, 7
celltype_overlap_summary.csv
profile_qc.csv
zscore_verification.csv
step2_profiles.rds
```

---

## Step 3 — Correlation matrices

**Script:** `03_step3_correlation_matrices.R`

### Definition

```text
Rows    = hsRA cell types (all annotated)
Columns = mouse model cell types (all annotated)
Values  = Spearman rho between 50-dimensional Hallmark z-score vectors
```

A correlation is computed only when at least **10** pathways are finite in both
vectors (`sum(ok) >= 10L`); otherwise the cell is `NA`. This threshold matters
only when many pathways were nulled for zero variance.

### Self-match rank

For each hsRA cell type, all available mouse cell types are ranked by Spearman
rho using `rank(-vals, ties.method = "min")`. Rank 1 means the same-name mouse
cell type is the top match.

**At this step the candidate set is all annotated mouse cell types.** Steps 5
and 8 recompute the rank on the five target cell types only, and it is the
recomputed version that appears in the figures. See § Verification findings.

### Outputs

```text
hsRA_vs_[model]_spearman_matrix.csv
all_hsRA_vs_mouse_spearman_long.csv
matched_celltype_correlations.csv
matched_vs_mismatched_summary.csv
self_match_ranks.csv                 <- full-space version, not plotted
model_level_summary.csv
step3_correlation_results.rds
```

---

## Step 4 — Permutation tests

**Script:** `04_step4_permutation_tests.R`

```r
set.seed(20260623L)
N_PERM_PAIRWISE    <- 1000L
N_PERM_GLOBAL      <- 10000L
MIN_PATHWAYS_FOR_COR <- 10L
```

### 4A — Pairwise pathway-label permutation

For each hsRA × mouse cell-type pair:

1. Hold the hsRA vector fixed.
2. Shuffle the mouse vector's pathway labels.
3. Recompute Spearman rho. Repeat 1,000 times.

The empirical p value is **two-sided**:

```r
p = (1 + sum(abs(null) >= abs(observed))) / (n_perm + 1)
```

Because of the `1 +` correction, the smallest attainable p is 1/1001 ≈ 9.99e-4.
A pairwise q value can therefore never fall far below 1e-3 regardless of how
strong the correlation is.

Two BH corrections are computed and both are exported:

| Column | Family |
|---|---|
| `q_value_BH_all_pairs` | All pairs across all models |
| `q_value_BH_within_model` | All pairs within each model |

Step 5 computes a **third** correction, `q_value_BH_target_matched_within_model`,
over only the five matched pairs within each model. This is the version that
carries the asterisks in the published dot plot.

### 4B — Global cell-label permutation

Per model:

1. Observed statistic = mean matched rho − mean mismatched rho.
2. Permute mouse cell-type labels across the correlation matrix.
3. Recompute. Repeat 10,000 times.

**Scope:** `mouse_types <- sort(unique(d$mouse_cell_type))` — the permutation
operates over the **full annotated cell-type space**, not the five target types.
Both `p_value_greater` (one-sided) and `p_value_two_sided` are exported; the
`-log10 global p` dimension in Figure Y uses `p_value_greater`.

**Resolution ceiling.** With 10,000 permutations the smallest attainable p is
1/10001, so `-log10 p` is capped at exactly **4.0000**. Any model reported at
4.00 is censored at the resolution limit, not measured at that value. Report as
`P < 1e-4` rather than as a point estimate.

### Outputs

```text
pairwise_pathway_permutation_results.csv
matched_celltype_significance.csv
global_matched_vs_mismatched_permutation.csv
global_null_distribution.csv           <- 10,000 rows per model, feeds Figure Y-C
step4_permutation_results.rds
```

---

## Step 5 — Figure X panels

**Script:** `05_step5_make_figures.R`

```r
TARGET_CELLTYPES        <- c("Macrophage", "Monocyte", "DC", "T cell", "FLS")
PREFERRED_MODEL_ORDER   <- c("C&AIA", "STIA", "GPI", "CIA")
```

The correlation table is filtered to the five target types on both axes, then
completed against a full 5 × 5 grid per model so that unavailable combinations
render as grey rather than disappearing.

Three quantities are **recomputed** at this step rather than inherited:

| Quantity | Recomputed as |
|---|---|
| BH q value | Over the five matched pairs within each model |
| Self-match rank | Among the five target cell types only |
| Matched-minus-mismatched | Over the 5 × 5 restricted matrix (descriptive only, no permutation) |

### Rendering conventions

| Element | Value |
|---|---|
| Diverging scale | low `#2166AC`, mid white, high `#B2182B`, limits −1 to 1 |
| Missing values | `grey85` |
| In-tile label colour | white when `abs(rho) >= 0.55`, otherwise black |
| Significance | `*` q ≤ 0.05, `**` q ≤ 0.01, `***` q ≤ 0.001 |
| Output | Vector PDF, `useDingbats = FALSE`, RGB |

### Panels

| Panel | File stem | Content |
|---|---|---|
| X-C | `F1_all_models_correlation_heatmap` | 5 × 5 rho heatmaps, faceted by model, diagonal outlined |
| — | `F1_[model]_correlation_heatmap` | Per-model enlarged version, for supplementary use |
| X-D | `F2_matched_celltype_correlation_dotplot` | Size = abs(rho), colour = rho, asterisks = BH q |
| — | `F3_restricted_matched_vs_mismatched_summary` | Descriptive bar chart, no test attached |
| X-E | `F4_self_match_rank_heatmap` | Rank 1 dark red, NA grey |

Panel X-A (per-dataset UMAPs) and Panel X-B (the 50-pathway z-score landscape
heatmap) are **not generated by any script in this pipeline**. They are produced
separately from the Seurat objects and the Step 2 z-score table respectively.

---

## Step 6 — Pathway driver analysis

**Script:** `06_step6_pathway_drivers.R`. Restricted to the five target cell
types. `TOP_N <- 10L` per model × cell type.

### Classification

```text
product_contribution = hsRA_z * mouse_z
```

| Class | Condition |
|---|---|
| Concordantly UP | `hsRA_z > 0` and `mouse_z > 0` |
| Concordantly DOWN | `hsRA_z < 0` and `mouse_z < 0` |
| Discordant | Opposite signs |
| Boundary/zero | Either z-score exactly 0 |
| Missing | Either z-score `NA` |

The last two classes exist in the implementation but are absent from the
original documentation and from the manuscript figure legend. They are normally
empty; confirm this in `pathway_concordance_summary.csv` rather than assuming
it.

Auxiliary quantities: `abs_difference = abs(hsRA_z - mouse_z)` and
`discordance_score`, defined as `abs(product_contribution)` for discordant
pathways only.

The model-level discordance metric consumed by Figure Y-B is the **mean of
`n_discordant` across the five cell types**, i.e. a mean count out of 50, not a
proportion.

### Outputs

```text
pathway_driver_all_pairs.csv
top_concordant_up_pathways.csv
top_concordant_down_pathways.csv
top_discordant_pathways.csv
pathway_concordance_summary.csv
model_pathway_summary.csv
pathway_driver_wide_for_heatmap.csv
target_dimension_coverage.csv
step6_pathway_driver_results.rds
```

---

## Step 7 — Sensitivity analysis

**Script:** `07_step7_sensitivity_analysis.R`. The script deletes its own prior
outputs via `unlink()` before running, so a failed partial run leaves no stale
files.

### Analyses

| Analysis name | Definition |
|---|---|
| `full_50_pathways` | All 50 Hallmark pathways, five cell types — the reference |
| `without_broad_pathways` | 40 pathways after removing the ten listed below |
| `leave_out_[cell type]` | All 50 pathways, four cell types (five variants) |

```r
BROAD_PATHWAYS <- c(
  "MYC_TARGETS_V1", "MYC_TARGETS_V2", "E2F_TARGETS", "G2M_CHECKPOINT",
  "OXIDATIVE_PHOSPHORYLATION", "MTORC1_SIGNALING", "DNA_REPAIR",
  "GLYCOLYSIS", "FATTY_ACID_METABOLISM", "PROTEIN_SECRETION"
)
```

### Critical implementation detail

**Z-scores are not recomputed for any sensitivity analysis.** Step 7 reads the
Step 2 matrix and subsets it.

- For `without_broad_pathways` this is correct: z-scoring is per-pathway across
  cell types, so dropping pathway columns leaves the remaining columns valid.
- For `leave_out_*` this is a deliberate but consequential choice: the omitted
  cell type is removed from the set of correlated pairs but **remains in the
  z-score normalization background**. The analysis therefore tests whether the
  mean matched correlation is driven by a single pair, not whether the
  normalization is robust to cell-type composition.

Any claim that "results do not depend on any single cellular compartment" must
be scoped accordingly. A stricter version of this test, re-running Step 2 on the
reduced cell-type set, is not implemented.

An immune-only pathway subset is deliberately excluded from the default
workflow, because classifying Hallmark pathways as immune or non-immune
introduces a subjective decision and reduces pathway coverage.

### Outputs

```text
pathway_set_definitions.csv
sensitivity_correlation_long.csv
sensitivity_model_summary.csv
sensitivity_leave_one_celltype_out_correlations.csv
sensitivity_leave_one_celltype_out_summary.csv
model_rank_by_each_sensitivity_analysis.csv
model_rank_stability.csv
target_dimension_coverage.csv
step7_sensitivity_results.rds
```

---

## Step 8 — Final summary and Figure Y

**Script:** `08_step8_final_summary_and_radar.R`

### The five evaluation dimensions, exactly as computed

| # | Dimension | Source | Definition |
|---|---|---|---|
| 1 | Mean matched rho | Step 3, restricted to 5 | Mean of available matched rho |
| 2 | Rank-1 fraction | Step 3 rho, rank recomputed on 5 | `n_rank1 / n_available_matched_celltypes` |
| 3 | `-log10` global p | **Step 4, full cell-type space** | `-log10(p_value_greater)` |
| 4 | Sensitivity retention | Step 7 | Mean of `pmin(ratio, 1)` over 6 checks |
| 5 | Discordant burden | Step 6, restricted to 5 | Mean `n_discordant` across cell types |

**Dimension 2 denominator.** The denominator is the number of *available*
matched cell types, not 5. A model with four available compartments and three
rank-1 matches scores 0.75, and a model with four available and four rank-1
matches would score 1.00 — identical to a model with five of five. Coverage and
rank quality are therefore separate facts and must be reported together.

**Dimension 3 scope mismatch.** This is the only dimension computed on the full
annotated cell-type space. Dimensions 1, 2, 4 and 5 use the five target types.
Either re-run Step 4B on the restricted matrix, or state the scope difference
explicitly in the figure legend.

**Dimension 4 capping.** `retention_ratio_capped = pmin(retention_ratio, 1)`.
A sensitivity analysis that *increases* the mean matched correlation is recorded
as 1.00, not as its true ratio. The reported retention is therefore bounded
above by 1 by construction and is a conservative summary, not a symmetric one.
The reference denominator is always `full_50_pathways` for that model, and the
mean is taken over six checks: one pathway-set analysis plus five
leave-one-cell-type-out analyses.

### Scaling for display

```r
minmax_scale(x)          # (x - min) / (max - min), across the four models
minmax_scale_reverse(x)  # applied to discordant burden
```

Scaling is computed **across the four models present**, so on every axis one
model is pinned at 1.0 and another at 0.0 regardless of the absolute spread. A
0.02 difference and a 0.5 difference produce identical radar geometry. Raw
values must be printed alongside, and the axes must never be summed or averaged
into a composite score.

### Fixed model colours

```text
C&AIA = #B2182B
STIA  = #2166AC
GPI   = #1B7837    (rename to GIA together with MODEL_ORDER)
CIA   = #762A83
```

### Panels

| Panel | Content | Underlying evidence |
|---|---|---|
| Y-A | Matched rho radar across the five cell types, one polygon per model, raw rho on a −1 to 1 scale | Step 3 |
| Y-B | Five-dimension dot plot and radar, min–max scaled positions with raw value labels | Steps 3, 4, 6, 7 |
| Y-C | Global permutation null distributions, one sub-panel per model, observed value as a vertical line | Step 4B |
| Y-D | Sensitivity heatmap, models × conditions, colour = mean matched rho | Step 7 |
| Y-E | Pathway concordance/discordance heatmap, colour = `product_contribution` | Step 6 |
| — | `FigureY_composite_A_to_E.pdf` | Assembled composite |

Missing dimensions, such as the STIA monocyte axis, remain missing and are never
imputed.

### Outputs

```text
final_model_summary.csv
final_matched_celltype_summary.csv
figure_y_panel_manifest.csv
figure_y_protocol_dimension_table.csv
figure_y_sensitivity_heatmap_table.csv
figure_y_pathway_heatmap_table.csv
figures/FigureY_A … FigureY_E, FigureY_composite_A_to_E.pdf
```

---

## Verification findings

Reconciliation of code, this document and the manuscript. Ordered by severity.

### Blocking — must be resolved before submission

**V1. Self-match rank candidate set is misreported in the manuscript.**
The Methods state that self-match ranks were determined by ranking *all
available mouse cell types*, and the Figure 6e legend repeats this. The panel is
generated by `05_step5_make_figures.R`, which filters to the five target cell
types **before** ranking. The published ranks are therefore rank-of-five, not
rank-of-all. Rank 1 among five candidates is a substantially weaker claim than
rank 1 among the full annotated set. Either correct both text locations, or
re-plot from `self_match_ranks.csv` (Step 3), which contains the full-space
ranks. The two versions may not agree.

**V2. Global permutation operates on a different cell-type space from every
other metric.** Step 4B permutes labels across all annotated mouse cell types,
while dimensions 1, 2, 4 and 5 of Figure 7b are computed on the five target
types. Figure 7d inherits the full-space statistic. This is defensible — a
larger label pool is a more stringent test — but it is currently undocumented in
both Methods and legend.

**V3. Model label `GPI` versus `GIA`.** Five scripts emit `GPI`; the manuscript
uses `GIA`. `MODEL_COLORS` is keyed by name, so a partial rename fails silently.

### Material — should be corrected in Methods

**V4. Z-score background is the full annotated cell-type set, not the five
target types.** The manuscript wording ("each pathway was then z-scored across
the available cell types") happens to be compatible with the code, but the
Discussion framework's phrase *five-cell z-score background* is incorrect and
must not be used. The correct scope statement is that normalization was
performed within each dataset across all annotated cell types, and that the
background composition therefore differs between models.

**V5. `-log10 P = 4.00` is a censoring artefact.** With 10,000 permutations and
the `1 +` correction the value is bounded at 4.0000. Reporting a range of
"3.22 to 4.00" implies a measured upper value. Report the ceiling as
`P < 1 × 10⁻⁴`.

**V6. Sensitivity retention is capped at 1.** Values above 1 are truncated
before averaging. This should appear in the Methods sentence describing the
metric.

**V7. Leave-one-cell-type-out does not re-normalize.** The dropped cell type
remains in the z-score background. Scope the robustness claim to the averaging,
not to the normalization.

**V8. Three different BH families exist in the outputs.** The asterisks in the
published dot plot come from a correction over the five matched pairs within
each model, computed in Step 5 — not from `q_value_BH_within_model` in Step 4.
The Methods sentence "Benjamini–Hochberg correction within each model" is
ambiguous between them and should name the family explicitly.

### Documentation gaps — now closed in this version

| Item | Value |
|---|---|
| Random seed | `20260623L`, set in Step 4 only |
| Minimum pathways for a correlation | 10 finite pairs |
| Pairwise permutation sidedness | Two-sided |
| Global permutation p used for Figure Y | One-sided `p_value_greater` |
| Zero-variance pathways | Set to `NA` in Step 2 |
| Top-N pathways per model × cell type | 10 |
| Direction classes | Five, not three |
| Step 4 → Step 5 chaining | Explicit `source()` at the end of Step 4 |
| Dependencies | `data.table`, `ggplot2` only |

### Documented but not implemented

| Claimed artefact | Status |
|---|---|
| Panel X-A per-dataset UMAPs | Generated outside this pipeline |
| Panel X-B 50-pathway landscape heatmap | Generated outside this pipeline |
| S1 raw correlation heatmaps without z-score | Not implemented; raw profiles are exported by Step 2 and the panel could be built from them |
| S3 KEGG cross-validation | Not implemented in any script |

---

## Key methodological decisions

1. AUCell is used because it is rank-based within each cell.
2. All 50 Hallmark pathways are used by default; no data-driven pathway
   selection is performed at any point.
3. Z-scoring is column-wise within each model, per pathway across cell types,
   using all annotated cell types as the background.
4. Disease-state-only comparison avoids asymmetric control definitions across
   species, at the cost of being unable to assess directional change relative to
   healthy tissue.
5. Partial cell-type coverage is allowed; missing dimensions are documented and
   never imputed.
6. Spearman correlation is used because it is rank-based and robust to scale
   differences.
7. No composite score is computed and no universally best model is declared.
8. Equal weighting is applied at the sample and source aggregation levels so
   that cell-count differences do not drive the profiles.

---

## Data requirements

| Field | Required | Notes |
|---|---|---|
| `model_id` | Yes | Model or dataset label |
| `cell_type` | Yes | Canonical name after harmonization |
| `sample_id` | Yes, or derivable | Used for equal-sample-weight means; a source with one sample is treated as n = 1 |
| 50 `HALLMARK_*` columns | Yes | Exactly 50; the script errors otherwise |

Minimum for z-scoring: at least two comparable cell types per model. No minimum
post-AUCell cell-count threshold is applied inside SPIRE.

---

## Scripts and execution

```text
01_step1_prepare_inputs.R
02_step2_make_profiles.R
03_step3_correlation_matrices.R
04_step4_permutation_tests.R      -> sources Step 5 on completion
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

All paths are absolute Windows paths defined at the top of each script. Step 8
probes two candidate directories (`Figure 6&7`, then `Figure 6`) and takes the
first that exists; the other scripts do not, so a directory rename breaks
Steps 1–7 while Step 8 continues to run against stale inputs.

### Checkpoints

Inspect these before advancing:

| After | File | Check |
|---|---|---|
| Step 1 | `model_celltype_availability.csv` | Unexpected cell-type labels from unmapped annotations |
| Step 2 | `zscore_verification.csv` | Column means ≈ 0, column SDs ≈ 1 |
| Step 2 | Console messages | Pathways nulled for zero variance |
| Step 3 | `matched_celltype_correlations.csv` | Matched rho plausible and non-missing |
| Step 4 | `global_matched_vs_mismatched_permutation.csv` | Whether any p sits at the 1/10001 floor |
| Step 7 | `model_rank_stability.csv` | Rank range across analyses |
| Step 8 | `figure_y_protocol_dimension_table.csv` | Raw values match the printed labels |
