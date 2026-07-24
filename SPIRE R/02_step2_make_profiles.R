options(stringsAsFactors = FALSE)

# Step 2 only (CORRECTED)
# Purpose:
# 1. Read Step 1 model-level AUCell profiles.
# 2. Build raw model x cell_type x Hallmark pathway profiles.
# 3. Compute z-score within each model_id, for each pathway across cell types.
# 4. Export cell-type overlap summaries for Step 3 correlation analysis.
#
# Z-score rule (CORRECTED):
# For each model_id separately, for each PATHWAY column:
#   z(cell_type) = (AUCell_mean(cell_type) - mean(all cell_types)) /
#                  sd(all cell_types)
#
# This normalizes each pathway across cell types within each model,
# removing the shared baseline that inflates all cross-cell-type correlations.
# After z-scoring, each pathway column has mean=0 and sd=1 across cell types.
#
# IMPORTANT: This is NOT row-wise z-score across pathways within a cell type.
# Row-wise z-score only tells you "which pathways are relatively active within
# this cell type" but does NOT remove the shared pattern across cell types,
# leading to inflated Spearman correlations for all pairs (~0.7-0.95).
#
# NORMALIZATION BACKGROUND (documented for Methods):
# The z-score background is ALL cell types annotated for that model, as
# present in model_equal_source_mean_AUCell.csv. It is NOT restricted to the
# five benchmarked compartments. A model annotated with Neutrophil, B cell,
# NK, Endothelial or Pericyte in addition to the five target types is
# normalized against all of them. Consequently the background composition
# differs between models, and the resulting correlations describe relative
# pathway preference within each model's own annotated compartment set.
# Do not describe this as a "five-cell-type z-score background".

required_pkgs <- c("data.table")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing R package(s): ", paste(missing_pkgs, collapse = ", "),
    "\nInstall first in RStudio: install.packages('data.table')",
    call. = FALSE
  )
}

library(data.table)

FIG6_DIR <- "E:/OneDrive/01.Academic/Papers/In Writing/C&AIA/Manuscript/8nd Nauture Communications/Figure 6"
STEP1_DIR <- file.path(FIG6_DIR, "RA-Scores", "01_step1_prepare_inputs")
OUT_DIR <- file.path(FIG6_DIR, "RA-Scores", "02_step2_make_profiles")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

step1_profile_file <- file.path(STEP1_DIR, "model_equal_source_mean_AUCell.csv")
if (!file.exists(step1_profile_file)) {
  stop(
    "Cannot find Step 1 profile file:\n",
    step1_profile_file,
    "\nRun 01_step1_prepare_inputs.R first.",
    call. = FALSE
  )
}

raw_profiles <- fread(step1_profile_file, showProgress = FALSE)

required_meta <- c("model_id", "cell_type")
missing_meta <- setdiff(required_meta, names(raw_profiles))
if (length(missing_meta) > 0) {
  stop(
    "Step 1 profile file is missing column(s): ",
    paste(missing_meta, collapse = ", "),
    call. = FALSE
  )
}

pathways <- grep("^HALLMARK_", names(raw_profiles), value = TRUE)
if (length(pathways) != 50L) {
  stop(
    "Expected exactly 50 HALLMARK_ columns; found ",
    length(pathways),
    call. = FALSE
  )
}

setorder(raw_profiles, model_id, cell_type)

# ──────────────────────────────────────────────────────────────
# CORRECTED z-score: column-wise (per pathway) across cell types
# within each model_id.
#
# scale() defaults to column-wise: for each column, subtract mean
# and divide by sd computed across all rows.
# We apply this separately per model_id so that human and each
# mouse model are z-scored independently.
# ──────────────────────────────────────────────────────────────

zscore_profiles <- copy(raw_profiles)

models <- sort(unique(zscore_profiles$model_id))

for (mid in models) {
  idx <- which(zscore_profiles$model_id == mid)
  
  if (length(idx) < 2L) {
    warning(
      "model_id '", mid, "' has only ", length(idx),
      " cell type(s); z-score requires >= 2. Setting to NA.",
      call. = FALSE
    )
    zscore_profiles[idx, (pathways) := NA_real_]
    next
  }
  
  mat <- as.matrix(zscore_profiles[idx, ..pathways])
  
  # scale(): for each COLUMN (pathway), z-score across ROWS (cell types)
  z_mat <- scale(mat)
  
  # Handle pathways with zero variance (constant across cell types)
  zero_var_cols <- which(attr(z_mat, "scaled:scale") == 0 | !is.finite(attr(z_mat, "scaled:scale")))
  if (length(zero_var_cols) > 0) {
    z_mat[, zero_var_cols] <- NA_real_
    message(
      "  model_id '", mid, "': ",
      length(zero_var_cols), " pathway(s) with zero variance across cell types set to NA: ",
      paste(pathways[zero_var_cols], collapse = ", ")
    )
  }
  
  zscore_profiles[idx, (pathways) := as.data.table(z_mat)]
}

# ──────────────────────────────────────────────────────────────
# Verification: after z-scoring, each pathway column within each
# model should have mean ≈ 0 and sd ≈ 1 across cell types.
# ──────────────────────────────────────────────────────────────

zscore_verification <- zscore_profiles[, {
  mat <- as.matrix(.SD)
  data.table(
    n_cell_types = nrow(mat),
    n_pathways = ncol(mat),
    mean_of_col_means = mean(colMeans(mat, na.rm = TRUE), na.rm = TRUE),
    mean_of_col_sds = mean(apply(mat, 2, sd, na.rm = TRUE), na.rm = TRUE)
  )
}, by = model_id, .SDcols = pathways]

message("\n--- Z-score verification (expect col_means ~ 0, col_sds ~ 1) ---")
print(zscore_verification)

# Check for any model where z-score looks wrong
bad_models <- zscore_verification[abs(mean_of_col_means) > 0.01 | abs(mean_of_col_sds - 1) > 0.05]
if (nrow(bad_models) > 0) {
  warning(
    "Z-score verification failed for model(s): ",
    paste(bad_models$model_id, collapse = ", "),
    "\nExpected column means ~ 0 and column sds ~ 1.",
    call. = FALSE
  )
}

# ──────────────────────────────────────────────────────────────
# Long format exports
# ──────────────────────────────────────────────────────────────

zscore_long <- melt(
  zscore_profiles,
  id.vars = c("model_id", "cell_type"),
  measure.vars = pathways,
  variable.name = "pathway",
  value.name = "zscore"
)

raw_long <- melt(
  raw_profiles,
  id.vars = c("model_id", "cell_type"),
  measure.vars = pathways,
  variable.name = "pathway",
  value.name = "AUCell_mean"
)

# ──────────────────────────────────────────────────────────────
# Cell-type overlap summary
# ──────────────────────────────────────────────────────────────

mouse_models <- setdiff(models, "hsRA")
hsra_celltypes <- sort(unique(raw_profiles[model_id == "hsRA", cell_type]))

celltype_overlap_summary <- rbindlist(lapply(mouse_models, function(m) {
  mouse_celltypes <- sort(unique(raw_profiles[model_id == m, cell_type]))
  shared <- intersect(hsra_celltypes, mouse_celltypes)
  
  data.table(
    mouse_model = m,
    n_hsRA_celltypes = length(hsra_celltypes),
    n_mouse_celltypes = length(mouse_celltypes),
    n_shared_celltypes = length(shared),
    shared_celltypes = paste(shared, collapse = "; "),
    hsRA_only_celltypes = paste(setdiff(hsra_celltypes, mouse_celltypes), collapse = "; "),
    mouse_only_celltypes = paste(setdiff(mouse_celltypes, hsra_celltypes), collapse = "; ")
  )
}))

# ──────────────────────────────────────────────────────────────
# Profile QC
# ──────────────────────────────────────────────────────────────

profile_qc <- raw_profiles[, {
  x <- as.numeric(.SD)
  data.table(
    n_pathways = length(x),
    n_missing_pathways = sum(is.na(x)),
    raw_mean_across_pathways = mean(x, na.rm = TRUE),
    raw_sd_across_pathways = stats::sd(x, na.rm = TRUE),
    zscore_possible = is.finite(stats::sd(x, na.rm = TRUE)) && stats::sd(x, na.rm = TRUE) > 0
  )
}, by = .(model_id, cell_type), .SDcols = pathways]

# ──────────────────────────────────────────────────────────────
# Export
# ──────────────────────────────────────────────────────────────

fwrite(raw_profiles, file.path(OUT_DIR, "raw_model_profiles_wide.csv"))
fwrite(raw_long, file.path(OUT_DIR, "raw_model_profiles_long.csv"))
fwrite(zscore_profiles, file.path(OUT_DIR, "zscore_model_profiles_wide.csv"))
fwrite(zscore_long, file.path(OUT_DIR, "zscore_model_profiles_long.csv"))
fwrite(celltype_overlap_summary, file.path(OUT_DIR, "celltype_overlap_summary.csv"))
fwrite(profile_qc, file.path(OUT_DIR, "profile_qc.csv"))
fwrite(zscore_verification, file.path(OUT_DIR, "zscore_verification.csv"))

saveRDS(
  list(
    raw_profiles = raw_profiles,
    zscore_profiles = zscore_profiles,
    pathways = pathways,
    celltype_overlap_summary = celltype_overlap_summary,
    profile_qc = profile_qc,
    zscore_verification = zscore_verification,
    zscore_rule = "within each model_id, for each pathway column, z-score across cell types (rows)"
  ),
  file.path(OUT_DIR, "step2_profiles.rds")
)

message("\nStep 2 complete.")
message("Input file: ", normalizePath(step1_profile_file, winslash = "/", mustWork = FALSE))
message("Output directory: ", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE))
message("Models: ", paste(models, collapse = ", "))
message("Hallmark pathways: ", length(pathways))
message("Z-score rule: within each model_id, for each pathway (column), z-score across cell types (rows).")
message("Check these files before Step 3:")
message("  - zscore_model_profiles_wide.csv")
message("  - zscore_verification.csv  <-- confirm col_means ~ 0, col_sds ~ 1")
message("  - celltype_overlap_summary.csv")
message("  - profile_qc.csv")