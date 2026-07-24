options(stringsAsFactors = FALSE)

# Step 6 only
# Purpose:
# Pathway driver analysis for matched hsRA vs mouse cell types.
#
# Scope:
# Only the five protocol-included cell types are analyzed:
#   Macrophage, Monocyte, DC, T cell, FLS
#
# For each mouse model x matched cell type:
#   compare the 50 Hallmark z-scores in hsRA vs the mouse model.
#
# Outputs:
#   - pathway_driver_all_pairs.csv
#   - top_concordant_up_pathways.csv
#   - top_concordant_down_pathways.csv
#   - top_discordant_pathways.csv
#   - pathway_concordance_summary.csv
#   - pathway_driver_wide_for_heatmap.csv

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
STEP4_DIR <- file.path(FIG6_DIR, "RA-Scores", "04_step4_permutation_tests")
OUT_DIR <- file.path(FIG6_DIR, "RA-Scores", "06_step6_pathway_drivers")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TARGET_CELLTYPES <- c("Macrophage", "Monocyte", "DC", "T cell", "FLS")
MODEL_ORDER <- c("C&AIA", "STIA", "GIA", "CIA")
TOP_N <- 10L

profile_file <- file.path(STEP2_DIR, "zscore_model_profiles_wide.csv")
sig_file <- file.path(STEP4_DIR, "matched_celltype_significance.csv")

if (!file.exists(profile_file)) {
  stop("Cannot find Step 2 z-score profile file:\n", profile_file, call. = FALSE)
}
if (!file.exists(sig_file)) {
  stop("Cannot find Step 4 matched significance file:\n", sig_file, call. = FALSE)
}

profiles <- fread(profile_file, showProgress = FALSE)
matched_sig <- fread(sig_file, showProgress = FALSE)

required_profile_cols <- c("model_id", "cell_type")
missing_profile_cols <- setdiff(required_profile_cols, names(profiles))
if (length(missing_profile_cols) > 0) {
  stop(
    "Step 2 profile file is missing column(s): ",
    paste(missing_profile_cols, collapse = ", "),
    call. = FALSE
  )
}

pathways <- grep("^HALLMARK_", names(profiles), value = TRUE)
if (length(pathways) != 50L) {
  stop("Expected exactly 50 HALLMARK_ columns; found ", length(pathways), call. = FALSE)
}

if (!"hsRA" %in% profiles$model_id) {
  stop("No hsRA rows found in Step 2 profile file.", call. = FALSE)
}

profile_long <- melt(
  profiles,
  id.vars = c("model_id", "cell_type"),
  measure.vars = pathways,
  variable.name = "pathway",
  value.name = "zscore"
)

profile_long <- profile_long[cell_type %in% TARGET_CELLTYPES]

hsra_long <- profile_long[model_id == "hsRA", .(
  cell_type,
  pathway,
  hsRA_z = zscore
)]

mouse_long <- profile_long[model_id != "hsRA", .(
  mouse_model = model_id,
  cell_type,
  pathway,
  mouse_z = zscore
)]

mouse_long <- mouse_long[mouse_model %in% MODEL_ORDER]

drivers <- merge(
  hsra_long,
  mouse_long,
  by = c("cell_type", "pathway"),
  all.y = TRUE,
  sort = FALSE
)

drivers[, matched_cell_type := TRUE]
drivers[, product_contribution := hsRA_z * mouse_z]
drivers[, abs_difference := abs(hsRA_z - mouse_z)]
drivers[, direction_class := fifelse(
  is.na(hsRA_z) | is.na(mouse_z), "Missing",
  fifelse(hsRA_z > 0 & mouse_z > 0, "Concordantly UP",
    fifelse(hsRA_z < 0 & mouse_z < 0, "Concordantly DOWN",
      fifelse(hsRA_z == 0 | mouse_z == 0, "Boundary/zero", "Discordant")
    )
  )
)]

drivers[, concordant := direction_class %in% c("Concordantly UP", "Concordantly DOWN")]
drivers[, discordance_score := fifelse(direction_class == "Discordant", abs(product_contribution), NA_real_)]

sig_cols <- intersect(
  c(
    "mouse_model", "human_cell_type", "mouse_cell_type",
    "observed_spearman_rho", "p_value_two_sided",
    "q_value_BH_within_model", "q_value_BH_all_pairs"
  ),
  names(matched_sig)
)

if (all(c("mouse_model", "human_cell_type", "mouse_cell_type") %in% sig_cols)) {
  sig_info <- matched_sig[
    human_cell_type %in% TARGET_CELLTYPES &
      mouse_cell_type %in% TARGET_CELLTYPES &
      human_cell_type == mouse_cell_type,
    ..sig_cols
  ]
  setnames(sig_info, c("human_cell_type", "mouse_cell_type"), c("cell_type", "mouse_cell_type_sig"))
  sig_info[, mouse_cell_type_sig := NULL]
  drivers <- merge(drivers, sig_info, by = c("mouse_model", "cell_type"), all.x = TRUE, sort = FALSE)
}

setcolorder(
  drivers,
  c(
    "mouse_model", "cell_type", "pathway",
    "hsRA_z", "mouse_z",
    "product_contribution", "abs_difference",
    "direction_class", "concordant",
    setdiff(names(drivers), c(
      "mouse_model", "cell_type", "pathway",
      "hsRA_z", "mouse_z",
      "product_contribution", "abs_difference",
      "direction_class", "concordant"
    ))
  )
)

setorder(drivers, mouse_model, cell_type, -product_contribution)

top_concordant_up <- drivers[
  direction_class == "Concordantly UP" & is.finite(product_contribution),
  head(.SD[order(-product_contribution)], TOP_N),
  by = .(mouse_model, cell_type)
]

top_concordant_down <- drivers[
  direction_class == "Concordantly DOWN" & is.finite(product_contribution),
  head(.SD[order(-product_contribution)], TOP_N),
  by = .(mouse_model, cell_type)
]

top_discordant <- drivers[
  direction_class == "Discordant" & is.finite(abs(product_contribution)),
  head(.SD[order(-abs(product_contribution))], TOP_N),
  by = .(mouse_model, cell_type)
]

pathway_concordance_summary <- drivers[, .(
  n_pathways = .N,
  n_concordant_up = sum(direction_class == "Concordantly UP", na.rm = TRUE),
  n_concordant_down = sum(direction_class == "Concordantly DOWN", na.rm = TRUE),
  n_discordant = sum(direction_class == "Discordant", na.rm = TRUE),
  concordant_fraction = mean(concordant, na.rm = TRUE),
  mean_product_contribution = mean(product_contribution, na.rm = TRUE),
  median_product_contribution = stats::median(product_contribution, na.rm = TRUE),
  mean_abs_difference = mean(abs_difference, na.rm = TRUE)
), by = .(mouse_model, cell_type)]

setorder(pathway_concordance_summary, mouse_model, cell_type)

model_pathway_summary <- drivers[, .(
  mean_product_contribution = mean(product_contribution, na.rm = TRUE),
  n_concordant_pairs = sum(concordant, na.rm = TRUE),
  n_discordant_pairs = sum(direction_class == "Discordant", na.rm = TRUE),
  concordant_fraction = mean(concordant, na.rm = TRUE)
), by = .(mouse_model, pathway)]

setorder(model_pathway_summary, mouse_model, -mean_product_contribution)

pathway_driver_wide_for_heatmap <- dcast(
  drivers,
  mouse_model + cell_type ~ pathway,
  value.var = "product_contribution"
)

target_dimension_coverage <- CJ(
  mouse_model = MODEL_ORDER,
  cell_type = TARGET_CELLTYPES,
  unique = TRUE
)
coverage_observed <- drivers[, .(
  n_pathways_available = sum(is.finite(hsRA_z) & is.finite(mouse_z))
), by = .(mouse_model, cell_type)]
target_dimension_coverage <- merge(
  target_dimension_coverage,
  coverage_observed,
  by = c("mouse_model", "cell_type"),
  all.x = TRUE,
  sort = FALSE
)
target_dimension_coverage[is.na(n_pathways_available), n_pathways_available := 0L]
target_dimension_coverage[, available := n_pathways_available > 0]

fwrite(drivers, file.path(OUT_DIR, "pathway_driver_all_pairs.csv"))
fwrite(top_concordant_up, file.path(OUT_DIR, "top_concordant_up_pathways.csv"))
fwrite(top_concordant_down, file.path(OUT_DIR, "top_concordant_down_pathways.csv"))
fwrite(top_discordant, file.path(OUT_DIR, "top_discordant_pathways.csv"))
fwrite(pathway_concordance_summary, file.path(OUT_DIR, "pathway_concordance_summary.csv"))
fwrite(model_pathway_summary, file.path(OUT_DIR, "model_pathway_summary.csv"))
fwrite(pathway_driver_wide_for_heatmap, file.path(OUT_DIR, "pathway_driver_wide_for_heatmap.csv"))
fwrite(target_dimension_coverage, file.path(OUT_DIR, "target_dimension_coverage.csv"))

saveRDS(
  list(
    pathway_driver_all_pairs = drivers,
    top_concordant_up_pathways = top_concordant_up,
    top_concordant_down_pathways = top_concordant_down,
    top_discordant_pathways = top_discordant,
    pathway_concordance_summary = pathway_concordance_summary,
    model_pathway_summary = model_pathway_summary,
    pathway_driver_wide_for_heatmap = pathway_driver_wide_for_heatmap,
    target_dimension_coverage = target_dimension_coverage,
    target_celltypes = TARGET_CELLTYPES,
    contribution_rule = "product_contribution = hsRA_z * mouse_z; positive means directional concordance, negative means discordance"
  ),
  file.path(OUT_DIR, "step6_pathway_driver_results.rds")
)

message("Step 6 complete.")
message("Input z-score file: ", normalizePath(profile_file, winslash = "/", mustWork = FALSE))
message("Input significance file: ", normalizePath(sig_file, winslash = "/", mustWork = FALSE))
message("Output directory: ", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE))
message("Target cell types: ", paste(TARGET_CELLTYPES, collapse = ", "))
message("Top N per model x cell type: ", TOP_N)
message("Check these files:")
message("  - pathway_driver_all_pairs.csv")
message("  - top_concordant_up_pathways.csv")
message("  - top_concordant_down_pathways.csv")
message("  - top_discordant_pathways.csv")
message("  - pathway_concordance_summary.csv")
message("  - model_pathway_summary.csv")
message("  - target_dimension_coverage.csv")
