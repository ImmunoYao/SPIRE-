options(stringsAsFactors = FALSE)

# Step 3 only
# Purpose:
# 1. Read Step 2 z-score profiles.
# 2. Compute Spearman correlation matrices:
#      hsRA cell_type x mouse model cell_type
#    using the 50 Hallmark z-score vector for each cell type.
# 3. Export:
#      - one wide matrix per mouse model
#      - one combined long-format correlation table
#      - matched vs mismatched summaries
#      - self-match rank table
#
# This step does NOT perform permutation tests or p-value correction.

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
OUT_DIR <- file.path(FIG6_DIR, "RA-Scores", "03_step3_correlation_matrices")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

profile_file <- file.path(STEP2_DIR, "zscore_model_profiles_wide.csv")
if (!file.exists(profile_file)) {
  stop(
    "Cannot find Step 2 z-score profile file:\n",
    profile_file,
    "\nRun 02_step2_make_profiles.R first.",
    call. = FALSE
  )
}

profiles <- fread(profile_file, showProgress = FALSE)

required_meta <- c("model_id", "cell_type")
missing_meta <- setdiff(required_meta, names(profiles))
if (length(missing_meta) > 0) {
  stop(
    "Step 2 profile file is missing column(s): ",
    paste(missing_meta, collapse = ", "),
    call. = FALSE
  )
}

pathways <- grep("^HALLMARK_", names(profiles), value = TRUE)
if (length(pathways) != 50L) {
  stop(
    "Expected exactly 50 HALLMARK_ columns; found ",
    length(pathways),
    call. = FALSE
  )
}

if (!"hsRA" %in% profiles$model_id) {
  stop("No hsRA profiles found in Step 2 file.", call. = FALSE)
}

make_profile_matrix <- function(d, model_value, pathways) {
  x <- d[d[["model_id"]] == model_value]
  if (nrow(x) == 0L) {
    stop("No profile rows for model_id = ", model_value, call. = FALSE)
  }

  if (anyDuplicated(x$cell_type)) {
    dup <- unique(x$cell_type[duplicated(x$cell_type)])
    stop(
      "Duplicated cell_type rows found for model_id = ", model_value, ": ",
      paste(dup, collapse = ", "),
      "\nStep 3 expects one row per model_id x cell_type.",
      call. = FALSE
    )
  }

  mat <- as.matrix(x[, ..pathways])
  rownames(mat) <- x$cell_type
  mat
}

spearman_matrix <- function(human_mat, mouse_mat) {
  ans <- matrix(
    NA_real_,
    nrow = nrow(human_mat),
    ncol = nrow(mouse_mat),
    dimnames = list(rownames(human_mat), rownames(mouse_mat))
  )

  for (i in seq_len(nrow(human_mat))) {
    for (j in seq_len(nrow(mouse_mat))) {
      h <- as.numeric(human_mat[i, ])
      m <- as.numeric(mouse_mat[j, ])
      ok <- is.finite(h) & is.finite(m)

      if (sum(ok) >= 10L) {
        ans[i, j] <- suppressWarnings(stats::cor(h[ok], m[ok], method = "spearman"))
      }
    }
  }

  ans
}

matrix_to_long <- function(mat, mouse_model) {
  out <- as.data.table(as.table(mat))
  setnames(out, c("human_cell_type", "mouse_cell_type", "spearman_rho"))
  out[, mouse_model := mouse_model]
  out[, matched_cell_type := human_cell_type == mouse_cell_type]
  setcolorder(out, c("mouse_model", "human_cell_type", "mouse_cell_type", "matched_cell_type", "spearman_rho"))
  out[]
}

human_mat <- make_profile_matrix(profiles, "hsRA", pathways)
mouse_models <- sort(setdiff(unique(profiles$model_id), "hsRA"))

if (length(mouse_models) == 0L) {
  stop("No mouse model profiles found in Step 2 file.", call. = FALSE)
}

correlation_long_list <- vector("list", length(mouse_models))
names(correlation_long_list) <- mouse_models

for (m in mouse_models) {
  mouse_mat <- make_profile_matrix(profiles, m, pathways)
  cor_mat <- spearman_matrix(human_mat, mouse_mat)

  matrix_out <- data.table(hsRA_cell_type = rownames(cor_mat), cor_mat, check.names = FALSE)
  fwrite(matrix_out, file.path(OUT_DIR, paste0("hsRA_vs_", m, "_spearman_matrix.csv")))

  correlation_long_list[[m]] <- matrix_to_long(cor_mat, m)
}

correlation_long <- rbindlist(correlation_long_list, use.names = TRUE)
setorder(correlation_long, mouse_model, human_cell_type, mouse_cell_type)

matched_correlations <- correlation_long[matched_cell_type == TRUE]
setorder(matched_correlations, mouse_model, human_cell_type)

matched_vs_mismatched_summary <- correlation_long[, .(
  n_pairs = sum(is.finite(spearman_rho)),
  mean_rho = mean(spearman_rho, na.rm = TRUE),
  median_rho = stats::median(spearman_rho, na.rm = TRUE),
  min_rho = min(spearman_rho, na.rm = TRUE),
  max_rho = max(spearman_rho, na.rm = TRUE)
), by = .(mouse_model, matched_cell_type)]
setorder(matched_vs_mismatched_summary, mouse_model, -matched_cell_type)

# SCOPE NOTE (documented for Methods):
# The self-match ranks computed here use ALL available mouse cell types as the
# candidate set. Steps 5 and 8 RECOMPUTE this rank restricted to
# TARGET_CELLTYPES (Macrophage, Monocyte, DC, T cell, FLS), and it is the
# restricted version that is plotted and reported. self_match_ranks.csv below
# is therefore the full-space reference, not the published figure source.
self_match_ranks <- correlation_long[, {
  vals <- spearman_rho
  names(vals) <- mouse_cell_type
  self_value <- vals[human_cell_type[1]]

  if (length(self_value) == 0L || is.na(self_value)) {
    data.table(
      self_match_available = FALSE,
      self_rho = NA_real_,
      self_rank = NA_integer_,
      n_candidate_mouse_celltypes = sum(is.finite(vals)),
      top_mouse_cell_type = names(vals)[which.max(vals)],
      top_rho = max(vals, na.rm = TRUE)
    )
  } else {
    data.table(
      self_match_available = TRUE,
      self_rho = as.numeric(self_value),
      self_rank = as.integer(rank(-vals, ties.method = "min")[human_cell_type[1]]),
      n_candidate_mouse_celltypes = sum(is.finite(vals)),
      top_mouse_cell_type = names(vals)[which.max(vals)],
      top_rho = max(vals, na.rm = TRUE)
    )
  }
}, by = .(mouse_model, human_cell_type)]
setorder(self_match_ranks, mouse_model, self_rank, human_cell_type)

model_level_summary <- self_match_ranks[, .(
  n_human_celltypes = .N,
  n_self_match_available = sum(self_match_available),
  mean_self_rho = mean(self_rho, na.rm = TRUE),
  median_self_rho = stats::median(self_rho, na.rm = TRUE),
  n_rank1_self_matches = sum(self_rank == 1L, na.rm = TRUE)
), by = mouse_model]
setorder(model_level_summary, -mean_self_rho)

fwrite(correlation_long, file.path(OUT_DIR, "all_hsRA_vs_mouse_spearman_long.csv"))
fwrite(matched_correlations, file.path(OUT_DIR, "matched_celltype_correlations.csv"))
fwrite(matched_vs_mismatched_summary, file.path(OUT_DIR, "matched_vs_mismatched_summary.csv"))
fwrite(self_match_ranks, file.path(OUT_DIR, "self_match_ranks.csv"))
fwrite(model_level_summary, file.path(OUT_DIR, "model_level_summary.csv"))

saveRDS(
  list(
    correlation_long = correlation_long,
    matched_correlations = matched_correlations,
    matched_vs_mismatched_summary = matched_vs_mismatched_summary,
    self_match_ranks = self_match_ranks,
    model_level_summary = model_level_summary,
    pathways = pathways,
    correlation_rule = "Spearman correlation between 50-Hallmark z-score vectors"
  ),
  file.path(OUT_DIR, "step3_correlation_results.rds")
)

message("Step 3 complete.")
message("Input file: ", normalizePath(profile_file, winslash = "/", mustWork = FALSE))
message("Output directory: ", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE))
message("Mouse models: ", paste(mouse_models, collapse = ", "))
message("Correlation rule: Spearman correlation between 50-Hallmark z-score vectors.")
message("Check these files before Step 4:")
message("  - hsRA_vs_[model]_spearman_matrix.csv")
message("  - all_hsRA_vs_mouse_spearman_long.csv")
message("  - matched_celltype_correlations.csv")
message("  - self_match_ranks.csv")
message("  - model_level_summary.csv")
