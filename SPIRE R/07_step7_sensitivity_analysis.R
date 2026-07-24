options(stringsAsFactors = FALSE)

# Step 7 only
# Purpose:
# Sensitivity analysis for hsRA vs mouse model similarity.
#
# Scope:
# Only the five protocol-included cell types are used:
#   Macrophage, Monocyte, DC, T cell, FLS
#
# Sensitivity analyses:
# 1. full_50_pathways
# 2. without_broad_pathways
# 3. leave-one-celltype-out using all 50 Hallmark pathways
#
# SCOPE NOTE (documented for Methods):
# Z-scores are NOT recomputed for any sensitivity analysis. This script reads
# the Step 2 matrix and subsets it.
#   - For the pathway-set analysis this is correct: z-scoring is per-pathway
#     across cell types, so dropping pathway columns leaves the remaining
#     columns valid.
#   - For leave-one-celltype-out the dropped cell type is removed from the
#     correlated pairs but REMAINS in the z-score normalization background.
#     These analyses therefore test whether the mean matched correlation is
#     driven by a single compartment, not whether the normalization itself is
#     robust to cell-type composition. Scope any robustness claim accordingly.
#
# This script reads the corrected Step 2 z-score matrix:
#   02_step2_make_profiles/zscore_model_profiles_wide.csv

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
STEP2_DIR <- file.path(FIG6_DIR, "RA-Scores", "02_step2_make_profiles")
OUT_DIR <- file.path(FIG6_DIR, "RA-Scores", "07_step7_sensitivity_analysis")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

old_outputs <- file.path(OUT_DIR, c(
  "pathway_set_definitions.csv",
  "sensitivity_correlation_long.csv",
  "sensitivity_model_summary.csv",
  "sensitivity_leave_one_celltype_out_correlations.csv",
  "sensitivity_leave_one_celltype_out_summary.csv",
  "model_rank_by_each_sensitivity_analysis.csv",
  "model_rank_stability.csv",
  "target_dimension_coverage.csv",
  "step7_sensitivity_results.rds"
))
unlink(old_outputs[file.exists(old_outputs)])

TARGET_CELLTYPES <- c("Macrophage", "Monocyte", "DC", "T cell", "FLS")
MODEL_ORDER <- c("C&AIA", "STIA", "GIA", "CIA")
MIN_PATHWAYS_FOR_COR <- 10L

BROAD_PATHWAYS <- paste0("HALLMARK_", c(
  "MYC_TARGETS_V1",
  "MYC_TARGETS_V2",
  "E2F_TARGETS",
  "G2M_CHECKPOINT",
  "OXIDATIVE_PHOSPHORYLATION",
  "MTORC1_SIGNALING",
  "DNA_REPAIR",
  "GLYCOLYSIS",
  "FATTY_ACID_METABOLISM",
  "PROTEIN_SECRETION"
))

profile_file <- file.path(STEP2_DIR, "zscore_model_profiles_wide.csv")
if (!file.exists(profile_file)) {
  stop("Cannot find Step 2 z-score profile file:\n", profile_file, call. = FALSE)
}

profiles <- fread(profile_file, showProgress = FALSE)

required_profile_cols <- c("model_id", "cell_type")
missing_profile_cols <- setdiff(required_profile_cols, names(profiles))
if (length(missing_profile_cols) > 0) {
  stop(
    "Step 2 profile file is missing column(s): ",
    paste(missing_profile_cols, collapse = ", "),
    call. = FALSE
  )
}

all_pathways <- grep("^HALLMARK_", names(profiles), value = TRUE)
if (length(all_pathways) != 50L) {
  stop("Expected exactly 50 HALLMARK_ columns; found ", length(all_pathways), call. = FALSE)
}

if (!"hsRA" %in% profiles$model_id) {
  stop("No hsRA rows found in Step 2 profile file.", call. = FALSE)
}

available_models <- intersect(MODEL_ORDER, sort(unique(profiles$model_id)))
available_models <- setdiff(available_models, "hsRA")

if (length(available_models) == 0L) {
  stop("No mouse model rows found in Step 2 profile file.", call. = FALSE)
}

pathway_sets <- list(
  full_50_pathways = all_pathways,
  without_broad_pathways = setdiff(all_pathways, intersect(BROAD_PATHWAYS, all_pathways))
)

pathway_set_definitions <- rbindlist(lapply(names(pathway_sets), function(nm) {
  data.table(
    analysis_name = nm,
    n_pathways = length(pathway_sets[[nm]]),
    pathway = pathway_sets[[nm]]
  )
}), use.names = TRUE)

spearman_one <- function(x, y, min_pathways = MIN_PATHWAYS_FOR_COR) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < min_pathways) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}

get_vector <- function(d, model_value, cell_value, pathways) {
  x <- d[model_id == model_value & cell_type == cell_value]
  if (nrow(x) == 0L) {
    return(rep(NA_real_, length(pathways)))
  }
  if (nrow(x) > 1L) {
    stop(
      "Duplicated profile rows for model_id = ", model_value,
      ", cell_type = ", cell_value,
      call. = FALSE
    )
  }
  as.numeric(x[, ..pathways])
}

compute_correlation_long <- function(d, pathways, analysis_name, celltypes = TARGET_CELLTYPES) {
  out <- vector("list", length(available_models) * length(celltypes) * length(celltypes))
  k <- 0L

  for (mouse_model in available_models) {
    for (human_ct in celltypes) {
      h <- get_vector(d, "hsRA", human_ct, pathways)

      for (mouse_ct in celltypes) {
        m <- get_vector(d, mouse_model, mouse_ct, pathways)
        k <- k + 1L

        out[[k]] <- data.table(
          analysis_name = analysis_name,
          n_pathways_in_analysis = length(pathways),
          mouse_model = mouse_model,
          human_cell_type = human_ct,
          mouse_cell_type = mouse_ct,
          matched_cell_type = human_ct == mouse_ct,
          spearman_rho = spearman_one(h, m),
          n_finite_pathways = sum(is.finite(h) & is.finite(m))
        )
      }
    }
  }

  rbindlist(out, use.names = TRUE)
}

summarize_model <- function(cor_dt, analysis_label) {
  cor_dt[, {
    matched <- spearman_rho[matched_cell_type == TRUE & is.finite(spearman_rho)]
    mismatched <- spearman_rho[matched_cell_type == FALSE & is.finite(spearman_rho)]

    self_rank_dt <- .SD[, {
      vals <- spearman_rho
      names(vals) <- mouse_cell_type
      finite_vals <- vals[is.finite(vals)]
      self_val <- vals[human_cell_type[1]]

      if (length(finite_vals) == 0L || length(self_val) == 0L || !is.finite(self_val)) {
        data.table(
          self_available = FALSE,
          self_rho = NA_real_,
          self_rank = NA_integer_
        )
      } else {
        data.table(
          self_available = TRUE,
          self_rho = as.numeric(self_val),
          self_rank = as.integer(rank(-finite_vals, ties.method = "min")[human_cell_type[1]])
        )
      }
    }, by = human_cell_type]

    data.table(
      analysis_name = analysis_label,
      n_matched_pairs = length(matched),
      n_mismatched_pairs = length(mismatched),
      mean_matched_rho = mean(matched, na.rm = TRUE),
      median_matched_rho = stats::median(matched, na.rm = TRUE),
      mean_mismatched_rho = mean(mismatched, na.rm = TRUE),
      matched_minus_mismatched = mean(matched, na.rm = TRUE) - mean(mismatched, na.rm = TRUE),
      n_self_match_available = sum(self_rank_dt$self_available),
      n_rank1_self_matches = sum(self_rank_dt$self_rank == 1L, na.rm = TRUE),
      mean_self_rho = mean(self_rank_dt$self_rho, na.rm = TRUE),
      median_self_rho = stats::median(self_rank_dt$self_rho, na.rm = TRUE)
    )
  }, by = mouse_model]
}

main_correlations <- rbindlist(lapply(names(pathway_sets), function(nm) {
  compute_correlation_long(
    d = profiles,
    pathways = pathway_sets[[nm]],
    analysis_name = nm,
    celltypes = TARGET_CELLTYPES
  )
}), use.names = TRUE)

main_model_summary <- rbindlist(lapply(names(pathway_sets), function(nm) {
  summarize_model(main_correlations[analysis_name == nm], nm)
}), use.names = TRUE)

main_model_summary[, rank_by_mean_matched_rho := frank(-mean_matched_rho, ties.method = "min"), by = analysis_name]
main_model_summary[, rank_by_matched_advantage := frank(-matched_minus_mismatched, ties.method = "min"), by = analysis_name]
setorder(main_model_summary, analysis_name, rank_by_mean_matched_rho, mouse_model)

leave_one_results <- rbindlist(lapply(TARGET_CELLTYPES, function(drop_ct) {
  kept <- setdiff(TARGET_CELLTYPES, drop_ct)
  compute_correlation_long(
    d = profiles,
    pathways = all_pathways,
    analysis_name = paste0("leave_out_", gsub(" ", "_", drop_ct)),
    celltypes = kept
  )[, left_out_cell_type := drop_ct][]
}), use.names = TRUE)

leave_one_summary <- rbindlist(lapply(unique(leave_one_results$analysis_name), function(nm) {
  left_out <- unique(leave_one_results[analysis_name == nm, left_out_cell_type])
  summarize_model(leave_one_results[analysis_name == nm], nm)[, left_out_cell_type := left_out][]
}), use.names = TRUE)

leave_one_summary[, rank_by_mean_matched_rho := frank(-mean_matched_rho, ties.method = "min"), by = analysis_name]
leave_one_summary[, rank_by_matched_advantage := frank(-matched_minus_mismatched, ties.method = "min"), by = analysis_name]
setcolorder(leave_one_summary, c("analysis_name", "left_out_cell_type", setdiff(names(leave_one_summary), c("analysis_name", "left_out_cell_type"))))
setorder(leave_one_summary, left_out_cell_type, rank_by_mean_matched_rho, mouse_model)

rank_stability_input <- rbind(
  main_model_summary[, .(
    analysis_group = "pathway_set",
    analysis_name,
    mouse_model,
    rank_by_mean_matched_rho,
    rank_by_matched_advantage,
    mean_matched_rho,
    matched_minus_mismatched
  )],
  leave_one_summary[, .(
    analysis_group = "leave_one_celltype_out",
    analysis_name,
    mouse_model,
    rank_by_mean_matched_rho,
    rank_by_matched_advantage,
    mean_matched_rho,
    matched_minus_mismatched
  )],
  use.names = TRUE
)

model_rank_stability <- rank_stability_input[, .(
  n_analyses = .N,
  median_rank_mean_matched_rho = stats::median(rank_by_mean_matched_rho, na.rm = TRUE),
  min_rank_mean_matched_rho = min(rank_by_mean_matched_rho, na.rm = TRUE),
  max_rank_mean_matched_rho = max(rank_by_mean_matched_rho, na.rm = TRUE),
  median_rank_matched_advantage = stats::median(rank_by_matched_advantage, na.rm = TRUE),
  min_rank_matched_advantage = min(rank_by_matched_advantage, na.rm = TRUE),
  max_rank_matched_advantage = max(rank_by_matched_advantage, na.rm = TRUE),
  mean_of_mean_matched_rho = mean(mean_matched_rho, na.rm = TRUE),
  mean_of_matched_advantage = mean(matched_minus_mismatched, na.rm = TRUE)
), by = mouse_model]
setorder(model_rank_stability, median_rank_mean_matched_rho, mouse_model)

target_dimension_coverage <- CJ(
  mouse_model = available_models,
  cell_type = TARGET_CELLTYPES,
  unique = TRUE
)
coverage_observed <- profiles[
  model_id %in% c("hsRA", available_models) &
    cell_type %in% TARGET_CELLTYPES,
  .(profile_available = TRUE),
  by = .(model_id, cell_type)
]

mouse_coverage <- coverage_observed[model_id != "hsRA", .(
  mouse_model = model_id,
  cell_type,
  mouse_profile_available = profile_available
)]

hsra_coverage <- coverage_observed[model_id == "hsRA", .(
  cell_type,
  hsRA_profile_available = profile_available
)]

target_dimension_coverage <- merge(target_dimension_coverage, mouse_coverage, by = c("mouse_model", "cell_type"), all.x = TRUE)
target_dimension_coverage <- merge(target_dimension_coverage, hsra_coverage, by = "cell_type", all.x = TRUE)
target_dimension_coverage[is.na(mouse_profile_available), mouse_profile_available := FALSE]
target_dimension_coverage[is.na(hsRA_profile_available), hsRA_profile_available := FALSE]
target_dimension_coverage[, comparable := mouse_profile_available & hsRA_profile_available]
setorder(target_dimension_coverage, mouse_model, cell_type)

fwrite(pathway_set_definitions, file.path(OUT_DIR, "pathway_set_definitions.csv"))
fwrite(main_correlations, file.path(OUT_DIR, "sensitivity_correlation_long.csv"))
fwrite(main_model_summary, file.path(OUT_DIR, "sensitivity_model_summary.csv"))
fwrite(leave_one_results, file.path(OUT_DIR, "sensitivity_leave_one_celltype_out_correlations.csv"))
fwrite(leave_one_summary, file.path(OUT_DIR, "sensitivity_leave_one_celltype_out_summary.csv"))
fwrite(rank_stability_input, file.path(OUT_DIR, "model_rank_by_each_sensitivity_analysis.csv"))
fwrite(model_rank_stability, file.path(OUT_DIR, "model_rank_stability.csv"))
fwrite(target_dimension_coverage, file.path(OUT_DIR, "target_dimension_coverage.csv"))

saveRDS(
  list(
    pathway_set_definitions = pathway_set_definitions,
    sensitivity_correlation_long = main_correlations,
    sensitivity_model_summary = main_model_summary,
    sensitivity_leave_one_celltype_out_correlations = leave_one_results,
    sensitivity_leave_one_celltype_out_summary = leave_one_summary,
    model_rank_by_each_sensitivity_analysis = rank_stability_input,
    model_rank_stability = model_rank_stability,
    target_dimension_coverage = target_dimension_coverage,
    target_celltypes = TARGET_CELLTYPES,
    pathway_sets = pathway_sets
  ),
  file.path(OUT_DIR, "step7_sensitivity_results.rds")
)

message("Step 7 complete.")
message("Input z-score file: ", normalizePath(profile_file, winslash = "/", mustWork = FALSE))
message("Output directory: ", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE))
message("Target cell types: ", paste(TARGET_CELLTYPES, collapse = ", "))
message("Pathway-set analyses: ", paste(names(pathway_sets), collapse = ", "))
message("Check these files:")
message("  - sensitivity_model_summary.csv")
message("  - sensitivity_leave_one_celltype_out_summary.csv")
message("  - model_rank_stability.csv")
message("  - model_rank_by_each_sensitivity_analysis.csv")
message("  - target_dimension_coverage.csv")
