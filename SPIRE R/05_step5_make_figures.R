options(stringsAsFactors = FALSE)

# Step 5 only
# Purpose:
# Make publication-oriented figures from Step 3 and Step 4 outputs.
#
# Inputs:
#   Step 3:
#     - all_hsRA_vs_mouse_spearman_long.csv
#     - self_match_ranks.csv
#     - model_level_summary.csv
#   Step 4:
#     - matched_celltype_significance.csv
#     - global_matched_vs_mismatched_permutation.csv
#
# Outputs:
#   - correlation heatmaps
#   - matched cell-type correlation dot plot
#   - restricted matched-vs-mismatched summary plot
#   - self-match rank heatmap
#   - figure-ready CSV tables

required_pkgs <- c("data.table", "ggplot2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing R package(s): ", paste(missing_pkgs, collapse = ", "),
    "\nInstall first in RStudio: install.packages(c('data.table', 'ggplot2'))",
    call. = FALSE
  )
}

library(data.table)
library(ggplot2)

FIG6_DIR <- "E:/OneDrive/01.Academic/Papers/In Writing/C&AIA/Manuscript/8nd Nauture Communications/Figure 6"
STEP3_DIR <- file.path(FIG6_DIR, "RA-Scores", "03_step3_correlation_matrices")
STEP4_DIR <- file.path(FIG6_DIR, "RA-Scores", "04_step4_permutation_tests")
OUT_DIR <- file.path(FIG6_DIR, "RA-Scores", "05_step5_figures")
FIG_DIR <- file.path(OUT_DIR, "figures")
TABLE_DIR <- file.path(OUT_DIR, "figure_tables")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

read_required <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot find required file:\n", path, call. = FALSE)
  }
  fread(path, showProgress = FALSE)
}

sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

sig_label <- function(q) {
  fifelse(
    is.na(q), "",
    fifelse(q <= 0.001, "***",
      fifelse(q <= 0.01, "**",
        fifelse(q <= 0.05, "*", "")
      )
    )
  )
}

save_plot <- function(plot, filename, width, height) {
  pdf_path <- file.path(FIG_DIR, paste0(filename, ".pdf"))

  vector_pdf_device <- function(file, width, height, bg = "white", ...) {
    grDevices::pdf(
      file = file,
      width = width,
      height = height,
      bg = bg,
      onefile = TRUE,
      useDingbats = FALSE,
      colormodel = "rgb"
    )
  }

  ggsave(
    pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = vector_pdf_device
  )
}

cor_long <- read_required(file.path(STEP3_DIR, "all_hsRA_vs_mouse_spearman_long.csv"))
self_ranks <- read_required(file.path(STEP3_DIR, "self_match_ranks.csv"))
model_summary <- read_required(file.path(STEP3_DIR, "model_level_summary.csv"))
matched_sig <- read_required(file.path(STEP4_DIR, "matched_celltype_significance.csv"))
global_perm <- read_required(file.path(STEP4_DIR, "global_matched_vs_mismatched_permutation.csv"))

required_cor_cols <- c("mouse_model", "human_cell_type", "mouse_cell_type", "matched_cell_type", "spearman_rho")
required_self_cols <- c("mouse_model", "human_cell_type", "self_match_available", "self_rho", "self_rank", "top_mouse_cell_type", "top_rho")
required_sig_cols <- c("mouse_model", "human_cell_type", "mouse_cell_type", "observed_spearman_rho", "q_value_BH_within_model")
required_global_cols <- c("mouse_model", "observed_matched_minus_mismatched", "p_value_greater", "q_value_BH_greater")

stopifnot(all(required_cor_cols %in% names(cor_long)))
stopifnot(all(required_self_cols %in% names(self_ranks)))
stopifnot(all(required_sig_cols %in% names(matched_sig)))
stopifnot(all(required_global_cols %in% names(global_perm)))

# Only display and summarize the protocol-included dimensions.
# Missing combinations are kept as NA so they appear as blank/grey cells.
TARGET_CELLTYPES <- c("Macrophage", "Monocyte", "DC", "T cell", "FLS")
PREFERRED_MODEL_ORDER <- c("C&AIA", "STIA", "GIA", "CIA")

model_order <- intersect(PREFERRED_MODEL_ORDER, sort(unique(cor_long$mouse_model)))
model_order <- c(model_order, setdiff(sort(unique(cor_long$mouse_model)), model_order))

cor_long <- cor_long[
  human_cell_type %in% TARGET_CELLTYPES &
    mouse_cell_type %in% TARGET_CELLTYPES
]

complete_cor_grid <- CJ(
  mouse_model = model_order,
  human_cell_type = TARGET_CELLTYPES,
  mouse_cell_type = TARGET_CELLTYPES,
  unique = TRUE
)

cor_long <- merge(
  complete_cor_grid,
  cor_long,
  by = c("mouse_model", "human_cell_type", "mouse_cell_type"),
  all.x = TRUE,
  sort = FALSE
)
cor_long[, matched_cell_type := human_cell_type == mouse_cell_type]

restricted_global_summary <- cor_long[, {
  matched <- spearman_rho[matched_cell_type == TRUE & is.finite(spearman_rho)]
  mismatched <- spearman_rho[matched_cell_type == FALSE & is.finite(spearman_rho)]
  data.table(
    n_matched_pairs = length(matched),
    n_mismatched_pairs = length(mismatched),
    mean_matched_rho = mean(matched, na.rm = TRUE),
    mean_mismatched_rho = mean(mismatched, na.rm = TRUE),
    observed_matched_minus_mismatched = mean(matched, na.rm = TRUE) - mean(mismatched, na.rm = TRUE)
  )
}, by = mouse_model]

matched_sig <- matched_sig[
  human_cell_type %in% TARGET_CELLTYPES &
    mouse_cell_type %in% TARGET_CELLTYPES &
    human_cell_type == mouse_cell_type
]

if ("p_value_two_sided" %in% names(matched_sig)) {
  matched_sig[, q_value_BH_target_matched_within_model := stats::p.adjust(p_value_two_sided, method = "BH"),
    by = mouse_model
  ]
} else {
  matched_sig[, q_value_BH_target_matched_within_model := q_value_BH_within_model]
}

complete_matched_grid <- CJ(
  mouse_model = model_order,
  human_cell_type = TARGET_CELLTYPES,
  unique = TRUE
)
complete_matched_grid[, mouse_cell_type := human_cell_type]

matched_sig <- merge(
  complete_matched_grid,
  matched_sig,
  by = c("mouse_model", "human_cell_type", "mouse_cell_type"),
  all.x = TRUE,
  sort = FALSE
)
matched_sig[, q_label := sig_label(q_value_BH_target_matched_within_model)]
matched_sig[, abs_rho_for_size := abs(observed_spearman_rho)]

# SCOPE NOTE (documented for Methods and the Fig. 6e legend):
# cor_long has already been filtered to TARGET_CELLTYPES on both axes, so the
# rank computed below is the position of the same-name mouse population among
# the five prespecified compartments available for that model (four for STIA,
# which lacks monocytes) - NOT among all annotated mouse cell types. This is
# the version that appears in the published figure. The full-space ranks are
# in 03_step3.../self_match_ranks.csv and are not plotted.
self_ranks_plot <- cor_long[, {
  vals <- spearman_rho
  names(vals) <- mouse_cell_type
  self_value <- vals[human_cell_type[1]]
  finite_vals <- vals[is.finite(vals)]

  if (length(finite_vals) == 0L || length(self_value) == 0L || !is.finite(self_value)) {
    data.table(
      self_match_available = FALSE,
      self_rho = NA_real_,
      self_rank = NA_integer_,
      n_candidate_mouse_celltypes = length(finite_vals),
      top_mouse_cell_type = NA_character_,
      top_rho = NA_real_
    )
  } else {
    data.table(
      self_match_available = TRUE,
      self_rho = as.numeric(self_value),
      self_rank = as.integer(rank(-finite_vals, ties.method = "min")[human_cell_type[1]]),
      n_candidate_mouse_celltypes = length(finite_vals),
      top_mouse_cell_type = names(finite_vals)[which.max(finite_vals)],
      top_rho = max(finite_vals, na.rm = TRUE)
    )
  }
}, by = .(mouse_model, human_cell_type)]

restricted_model_summary <- self_ranks_plot[, .(
  n_target_human_celltypes = .N,
  n_self_match_available = sum(self_match_available),
  mean_self_rho = mean(self_rho, na.rm = TRUE),
  median_self_rho = stats::median(self_rho, na.rm = TRUE),
  n_rank1_self_matches = sum(self_rank == 1L, na.rm = TRUE)
), by = mouse_model]

cor_long[, mouse_model := factor(mouse_model, levels = model_order)]
cor_long[, human_cell_type := factor(human_cell_type, levels = rev(TARGET_CELLTYPES))]
cor_long[, mouse_cell_type := factor(mouse_cell_type, levels = TARGET_CELLTYPES)]
cor_long[, rho_label := ifelse(is.finite(spearman_rho), sprintf("%.2f", spearman_rho), "")]
cor_long[, rho_text_color := ifelse(is.finite(spearman_rho) & abs(spearman_rho) >= 0.55, "white", "black")]

matched_sig[, mouse_model := factor(mouse_model, levels = model_order)]
matched_sig[, human_cell_type := factor(human_cell_type, levels = rev(TARGET_CELLTYPES))]

restricted_global_summary[, mouse_model := factor(mouse_model, levels = model_order)]

self_ranks_plot[, mouse_model := factor(mouse_model, levels = model_order)]
self_ranks_plot[, human_cell_type := factor(human_cell_type, levels = rev(TARGET_CELLTYPES))]

# Figure 1: all model correlation heatmap
fig1 <- ggplot(cor_long, aes(x = mouse_cell_type, y = human_cell_type, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_tile(
    data = cor_long[matched_cell_type == TRUE],
    aes(x = mouse_cell_type, y = human_cell_type),
    fill = NA,
    color = "black",
    linewidth = 0.45
  ) +
  geom_text(
    aes(label = rho_label, color = rho_text_color),
    size = 2.4,
    na.rm = TRUE,
    show.legend = FALSE
  ) +
  facet_wrap(~ mouse_model, scales = "free_x") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey85",
    name = "Spearman rho",
    guide = guide_colourbar(raster = FALSE)
  ) +
  scale_color_identity() +
  labs(
    title = "hsRA vs mouse model cell-type correlation",
    x = "Mouse model cell type",
    y = "hsRA cell type"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_rect(fill = "grey92", color = "grey60"),
    plot.title = element_text(hjust = 0)
  )

save_plot(fig1, "F1_all_models_correlation_heatmap", width = 12, height = 8)

# Figure 1 per model
for (m in levels(cor_long$mouse_model)) {
  d <- cor_long[mouse_model == m]
  if (nrow(d) == 0L) next

  fig_m <- ggplot(d, aes(x = mouse_cell_type, y = human_cell_type, fill = spearman_rho)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_tile(
      data = d[matched_cell_type == TRUE],
      aes(x = mouse_cell_type, y = human_cell_type),
      fill = NA,
      color = "black",
      linewidth = 0.55
    ) +
    geom_text(
      aes(label = rho_label, color = rho_text_color),
      size = 3,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "grey85",
      name = "Spearman rho",
      guide = guide_colourbar(raster = FALSE)
    ) +
    scale_color_identity() +
    labs(
      title = paste0("hsRA vs ", m, " correlation matrix"),
      x = paste0(m, " cell type"),
      y = "hsRA cell type"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(hjust = 0)
    )

  save_plot(fig_m, paste0("F1_", sanitize_filename(as.character(m)), "_correlation_heatmap"), width = 7, height = 6)
}

# Figure 2: matched cell-type correlations with BH significance
fig2 <- ggplot(matched_sig, aes(x = mouse_model, y = human_cell_type)) +
  geom_point(aes(size = abs_rho_for_size, color = observed_spearman_rho), na.rm = TRUE) +
  geom_text(aes(label = q_label), nudge_x = 0.25, size = 4, na.rm = TRUE) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey70",
    name = "Matched rho",
    guide = guide_colourbar(raster = FALSE)
  ) +
  scale_size_continuous(range = c(2, 7), name = "|Matched rho|") +
  labs(
    title = "Matched cell-type correlations",
    subtitle = "* q<=0.05, ** q<=0.01, *** q<=0.001; BH adjusted within each model",
    x = "Mouse model",
    y = "Cell type"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0)
  )

save_plot(fig2, "F2_matched_celltype_correlation_dotplot", width = 8, height = 6)

# Figure 3: restricted matched-vs-mismatched descriptive summary
fig3 <- ggplot(restricted_global_summary, aes(x = mouse_model, y = observed_matched_minus_mismatched)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_col(fill = "#B2182B", width = 0.65) +
  labs(
    title = "Matched-vs-mismatched advantage in five included cell types",
    subtitle = "Restricted to Macrophage, Monocyte, DC, T cell and FLS; statistic = mean matched rho - mean mismatched rho",
    x = "Mouse model",
    y = "Matched - mismatched Spearman rho"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0)
  )

save_plot(fig3, "F3_restricted_matched_vs_mismatched_summary", width = 7, height = 5)

# Figure 4: self-match rank heatmap
self_ranks_plot[self_match_available == FALSE, self_rank := NA_integer_]

fig4 <- ggplot(self_ranks_plot, aes(x = mouse_model, y = human_cell_type, fill = self_rank)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(is.na(self_rank), "NA", as.character(self_rank))), size = 3) +
  scale_fill_gradient(
    low = "#B2182B",
    high = "#FEE8C8",
    na.value = "grey85",
    name = "Self-match rank",
    guide = guide_colourbar(raster = FALSE)
  ) +
  labs(
    title = "Self-match rank of each hsRA cell type",
    subtitle = "Rank 1 means the same-name mouse cell type is the top match",
    x = "Mouse model",
    y = "hsRA cell type"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0)
  )

save_plot(fig4, "F4_self_match_rank_heatmap", width = 7, height = 6)

# Export figure-ready tables
fwrite(cor_long, file.path(TABLE_DIR, "F1_correlation_heatmap_table.csv"))
fwrite(matched_sig, file.path(TABLE_DIR, "F2_matched_celltype_dotplot_table.csv"))
fwrite(restricted_global_summary, file.path(TABLE_DIR, "F3_restricted_matched_vs_mismatched_table.csv"))
fwrite(self_ranks_plot, file.path(TABLE_DIR, "F4_self_match_rank_table.csv"))
fwrite(restricted_model_summary, file.path(TABLE_DIR, "restricted_model_level_summary.csv"))

message("Step 5 complete.")
message("Output directory: ", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE))
message("Figures:")
message("  - figures/F1_all_models_correlation_heatmap.pdf")
message("  - figures/F1_[model]_correlation_heatmap.pdf")
message("  - figures/F2_matched_celltype_correlation_dotplot.pdf")
message("  - figures/F3_restricted_matched_vs_mismatched_summary.pdf")
message("  - figures/F4_self_match_rank_heatmap.pdf")
message("Figure-ready tables:")
message("  - figure_tables/F1_correlation_heatmap_table.csv")
message("  - figure_tables/F2_matched_celltype_dotplot_table.csv")
message("  - figure_tables/F3_restricted_matched_vs_mismatched_table.csv")
message("  - figure_tables/F4_self_match_rank_table.csv")
