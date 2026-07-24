options(stringsAsFactors = FALSE)

# Step 8 only
# Purpose:
# Generate protocol-aligned final summaries and Figure Y outputs.
#
# This script keeps existing Figure X panels as reusable evidence files and
# creates the Figure Y model-evaluation panels:
#   Y-B protocol-dimension radar plot
#   Y-C global permutation null distributions
#   Y-D sensitivity analysis heatmap
#   Y-E pathway concordance/discordance heatmap
#
# Main outputs:
#   - final_model_summary.csv
#   - final_matched_celltype_summary.csv
#   - figure_y_panel_manifest.csv
#   - figure_y_protocol_dimension_table.csv
#   - figure_y_sensitivity_heatmap_table.csv
#   - figure_y_pathway_heatmap_table.csv
#   - figures/FigureY_A_celltype_matched_rho_radar.pdf
#   - figures/FigureY_B_protocol_dimension_radar.pdf
#   - figures/FigureY_B_protocol_dimension_dotplot.pdf
#   - figures/FigureY_B_protocol_dimension_radar_[model].pdf
#   - figures/FigureY_B_protocol_dimension_individual_radars.pdf
#   - figures/FigureY_C_global_permutation_null_distributions.pdf
#   - figures/FigureY_D_sensitivity_analysis_heatmap.pdf
#   - figures/FigureY_E_pathway_concordance_discordance_heatmap.pdf
#   - figures/FigureY_composite_A_to_E.pdf

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

FIG6_CANDIDATE_DIRS <- c(
  "E:/OneDrive/01.Academic/Papers/In Writing/C&AIA/Manuscript/8nd Nauture Communications/Figure 6&7",
  "E:/OneDrive/01.Academic/Papers/In Writing/C&AIA/Manuscript/8nd Nauture Communications/Figure 6"
)
FIG6_DIR <- FIG6_CANDIDATE_DIRS[dir.exists(FIG6_CANDIDATE_DIRS)][1]
if (is.na(FIG6_DIR)) {
  stop(
    "Cannot find Figure 6/6&7 directory. Checked:\n",
    paste(FIG6_CANDIDATE_DIRS, collapse = "\n"),
    call. = FALSE
  )
}
RA_DIR <- file.path(FIG6_DIR, "RA-Scores")
STEP3_DIR <- file.path(RA_DIR, "03_step3_correlation_matrices")
STEP4_DIR <- file.path(RA_DIR, "04_step4_permutation_tests")
STEP5_FIG_DIR <- file.path(RA_DIR, "05_step5_figures", "figures")
STEP6_DIR <- file.path(RA_DIR, "06_step6_pathway_drivers")
STEP7_DIR <- file.path(RA_DIR, "07_step7_sensitivity_analysis")
STEP8_DIR <- file.path(RA_DIR, "08_step8_final_summary_and_radar")
STEP8_FIG_DIR <- file.path(STEP8_DIR, "figures")

OUT_DIR <- file.path(RA_DIR, "08_step8_final_summary_and_radar")
FIG_DIR <- file.path(OUT_DIR, "figures")
TABLE_DIR <- file.path(OUT_DIR, "tables")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

TARGET_CELLTYPES <- c("Macrophage", "Monocyte", "DC", "T cell", "FLS")
MODEL_ORDER <- c("C&AIA", "STIA", "GIA", "CIA")

MODEL_COLORS <- c(
  "C&AIA" = "#B2182B",
  "STIA" = "#2166AC",
  "GIA" = "#1B7837",
  "CIA" = "#762A83"
)

read_required <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot find required file:\n", path, call. = FALSE)
  }
  fread(path, showProgress = FALSE)
}

save_vector_pdf <- function(plot, filename, width, height) {
  pdf_path <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  grDevices::pdf(
    file = pdf_path,
    width = width,
    height = height,
    onefile = TRUE,
    useDingbats = FALSE,
    colormodel = "rgb",
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()
}

copy_existing_pdf <- function(source_path, target_name) {
  target_path <- file.path(FIG_DIR, target_name)
  if (file.exists(source_path)) {
    file.copy(source_path, target_path, overwrite = TRUE)
  }
  target_path
}

minmax_scale <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    return(rep(0.5, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

minmax_scale_reverse <- function(x) {
  y <- minmax_scale(x)
  ifelse(is.na(y), NA_real_, 1 - y)
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

raw_label <- function(metric, value) {
  ifelse(
    is.na(value), "NA",
    ifelse(
      metric %in% c("mean_matched_rho", "sensitivity_retention_ratio"),
      sprintf("%.3f", value),
      ifelse(
        metric == "rank1_self_match_fraction",
        sprintf("%.2f", value),
        ifelse(metric == "global_structure_neg_log10_p", sprintf("%.2f", value), sprintf("%.1f", value))
      )
    )
  )
}

make_protocol_dimension_radar <- function(dt, model_levels, model_colors,
                                          title, subtitle) {
  d <- copy(dt)
  d[, mouse_model := factor(as.character(mouse_model), levels = model_levels)]
  d[, axis_label := as.character(metric_label)]
  axis_levels <- rev(levels(dt$metric_label))
  d[, axis_label := factor(axis_label, levels = axis_levels)]

  axis_meta <- data.table(
    axis_label = factor(axis_levels, levels = axis_levels),
    axis_index = seq_along(axis_levels)
  )
  single_axis_labels <- c(
    "Mean Spearman rho" = "Mean rho",
    "Rank=1 fraction" = "Rank=1\nfraction",
    "-log10 global p" = "-log10 P",
    "Sensitivity retention" = "Sensitivity\nretention",
    "Low discordance" = "Low\ndiscordance"
  )
  axis_meta[, axis_display_label := single_axis_labels[as.character(axis_label)]]
  axis_meta[is.na(axis_display_label), axis_display_label := as.character(axis_label)]
  axis_meta[, angle := pi / 2 - 2 * pi * (axis_index - 1) / length(axis_levels)]
  axis_meta[, x_end := cos(angle)]
  axis_meta[, y_end := sin(angle)]
  axis_meta[, x_label := 1.22 * cos(angle)]
  axis_meta[, y_label := 1.22 * sin(angle)]
  axis_meta[, hjust := fifelse(x_label > 0.1, 0, fifelse(x_label < -0.1, 1, 0.5))]
  axis_meta[, vjust := fifelse(y_label > 0.1, 0, fifelse(y_label < -0.1, 1, 0.5))]

  d <- merge(d, axis_meta[, .(axis_label, axis_index, angle)], by = "axis_label", all.x = TRUE, sort = FALSE)
  d[, radius := pmax(0, pmin(1, scaled_value))]
  d[!is.finite(radius), radius := NA_real_]
  d[, x := radius * cos(angle)]
  d[, y := radius * sin(angle)]
  setorder(d, mouse_model, axis_index)

  closed <- d[is.finite(x) & is.finite(y), {
    if (.N < 3L) .SD else rbind(.SD, .SD[1])
  }, by = mouse_model]

  grid <- rbindlist(lapply(c(0.25, 0.5, 0.75, 1), function(r) {
    theta <- seq(0, 2 * pi, length.out = 241)
    data.table(ring = r, x = r * cos(theta), y = r * sin(theta))
  }))
  ring_labels <- data.table(
    ring = c(0.25, 0.5, 0.75, 1),
    label = c("0.25", "0.50", "0.75", "1.00"),
    x = 0.04,
    y = c(0.25, 0.5, 0.75, 1)
  )

  ggplot() +
    geom_path(data = grid, aes(x = x, y = y, group = ring), color = "grey84", linewidth = 0.3) +
    geom_segment(
      data = axis_meta,
      aes(x = 0, y = 0, xend = x_end, yend = y_end),
      color = "grey84",
      linewidth = 0.3
    ) +
    geom_text(
      data = axis_meta,
      aes(x = x_label, y = y_label, label = axis_label, hjust = hjust, vjust = vjust),
      size = 3,
      color = "grey25"
    ) +
    geom_text(
      data = ring_labels,
      aes(x = x, y = y, label = label),
      size = 2.4,
      color = "grey45",
      hjust = 0
    ) +
    geom_polygon(
      data = closed,
      aes(x = x, y = y, group = mouse_model, color = mouse_model, fill = mouse_model),
      alpha = 0.08,
      linewidth = 0.7,
      na.rm = TRUE
    ) +
    geom_path(
      data = closed,
      aes(x = x, y = y, group = mouse_model, color = mouse_model),
      linewidth = 0.9,
      na.rm = TRUE
    ) +
    geom_point(
      data = d,
      aes(x = x, y = y, color = mouse_model),
      size = 2,
      na.rm = TRUE
    ) +
    coord_equal(xlim = c(-1.45, 1.55), ylim = c(-1.35, 1.35), clip = "off") +
    scale_color_manual(values = model_colors, drop = FALSE) +
    scale_fill_manual(values = model_colors, drop = FALSE) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_void(base_size = 10) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0, size = 13),
      plot.subtitle = element_text(hjust = 0, size = 9),
      plot.margin = margin(12, 36, 12, 36)
    )
}

make_single_protocol_dimension_radar <- function(dt, model_name, model_color,
                                                 title = model_name) {
  d <- copy(dt[as.character(mouse_model) == model_name])
  d[, axis_label := as.character(metric_label)]
  axis_levels <- rev(levels(dt$metric_label))
  d[, axis_label := factor(axis_label, levels = axis_levels)]

  axis_meta <- data.table(
    axis_label = factor(axis_levels, levels = axis_levels),
    axis_index = seq_along(axis_levels)
  )
  single_axis_labels <- c(
    "Mean Spearman rho" = "Mean rho",
    "Rank=1 fraction" = "Rank=1\nfraction",
    "-log10 global p" = "-log10 P",
    "Sensitivity retention" = "Sensitivity\nretention",
    "Low discordance" = "Low\ndiscordance"
  )
  axis_meta[, axis_display_label := single_axis_labels[as.character(axis_label)]]
  axis_meta[is.na(axis_display_label), axis_display_label := as.character(axis_label)]
  axis_meta[, angle := pi / 2 - 2 * pi * (axis_index - 1) / length(axis_levels)]
  axis_meta[, x_end := cos(angle)]
  axis_meta[, y_end := sin(angle)]
  axis_meta[, x_label := 1.23 * cos(angle)]
  axis_meta[, y_label := 1.23 * sin(angle)]
  axis_meta[, hjust := fifelse(x_label > 0.1, 0, fifelse(x_label < -0.1, 1, 0.5))]
  axis_meta[, vjust := fifelse(y_label > 0.1, 0, fifelse(y_label < -0.1, 1, 0.5))]

  d <- merge(d, axis_meta[, .(axis_label, axis_index, angle)], by = "axis_label", all.x = TRUE, sort = FALSE)
  d[, radius := pmax(0, pmin(1, scaled_value))]
  d[!is.finite(radius), radius := NA_real_]
  d[, x := radius * cos(angle)]
  d[, y := radius * sin(angle)]
  d[, value_x := pmin(1.10, radius + 0.07) * cos(angle)]
  d[, value_y := pmin(1.10, radius + 0.07) * sin(angle)]
  d[, value_hjust := fifelse(value_x > 0.1, 0, fifelse(value_x < -0.1, 1, 0.5))]
  d[, value_vjust := fifelse(value_y > 0.1, 0, fifelse(value_y < -0.1, 1, 0.5))]
  setorder(d, axis_index)

  closed <- d[is.finite(x) & is.finite(y)]
  if (nrow(closed) >= 3L) {
    closed <- rbind(closed, closed[1])
  }

  grid <- rbindlist(lapply(c(0.25, 0.5, 0.75, 1), function(r) {
    theta <- seq(0, 2 * pi, length.out = 241)
    data.table(ring = r, x = r * cos(theta), y = r * sin(theta))
  }))

  ggplot() +
    geom_path(data = grid, aes(x = x, y = y, group = ring), color = "grey78", linewidth = 0.35) +
    geom_segment(
      data = axis_meta,
      aes(x = 0, y = 0, xend = x_end, yend = y_end),
      color = "grey78",
      linewidth = 0.35
    ) +
    geom_polygon(
      data = closed,
      aes(x = x, y = y),
      fill = model_color,
      color = model_color,
      alpha = 0.86,
      linewidth = 0.9,
      na.rm = TRUE
    ) +
    geom_path(
      data = closed,
      aes(x = x, y = y),
      color = model_color,
      linewidth = 1,
      na.rm = TRUE
    ) +
    geom_point(
      data = d,
      aes(x = x, y = y),
      color = "white",
      fill = "white",
      size = 2,
      na.rm = TRUE
    ) +
    geom_text(
      data = d,
      aes(x = value_x, y = value_y, label = raw_value_label, hjust = value_hjust, vjust = value_vjust),
      color = model_color,
      fontface = "bold",
      size = 3.8,
      na.rm = TRUE
    ) +
    geom_text(
      data = axis_meta,
      aes(x = x_label, y = y_label, label = axis_display_label, hjust = hjust, vjust = vjust),
      color = "grey25",
      fontface = "bold",
      size = 4.2
    ) +
    coord_equal(xlim = c(-1.45, 1.55), ylim = c(-1.35, 1.35), clip = "off") +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(color = model_color, hjust = 0.5, size = 16, face = "bold"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(14, 28, 14, 28)
    )
}

cor_long <- read_required(file.path(STEP3_DIR, "all_hsRA_vs_mouse_spearman_long.csv"))
matched_sig <- read_required(file.path(STEP4_DIR, "matched_celltype_significance.csv"))
global_perm <- read_required(file.path(STEP4_DIR, "global_matched_vs_mismatched_permutation.csv"))
global_null <- read_required(file.path(STEP4_DIR, "global_null_distribution.csv"))
pathway_summary <- read_required(file.path(STEP6_DIR, "pathway_concordance_summary.csv"))
model_pathway_summary <- read_required(file.path(STEP6_DIR, "model_pathway_summary.csv"))
pathway_driver <- read_required(file.path(STEP6_DIR, "pathway_driver_all_pairs.csv"))
sensitivity_summary <- read_required(file.path(STEP7_DIR, "sensitivity_model_summary.csv"))
leave_one_summary <- read_required(file.path(STEP7_DIR, "sensitivity_leave_one_celltype_out_summary.csv"))
rank_stability <- read_required(file.path(STEP7_DIR, "model_rank_stability.csv"))

cor_long <- cor_long[
  mouse_model %in% MODEL_ORDER &
    human_cell_type %in% TARGET_CELLTYPES &
    mouse_cell_type %in% TARGET_CELLTYPES
]

matched_cor <- cor_long[human_cell_type == mouse_cell_type, .(
  mouse_model,
  cell_type = human_cell_type,
  matched_rho = spearman_rho
)]

matched_sig_target <- matched_sig[
  mouse_model %in% MODEL_ORDER &
    human_cell_type %in% TARGET_CELLTYPES &
    mouse_cell_type %in% TARGET_CELLTYPES &
    human_cell_type == mouse_cell_type
]

sig_cols <- intersect(
  c(
    "mouse_model", "human_cell_type", "observed_spearman_rho",
    "p_value_two_sided", "q_value_BH_within_model", "q_value_BH_all_pairs"
  ),
  names(matched_sig_target)
)
matched_sig_target <- matched_sig_target[, ..sig_cols]
setnames(matched_sig_target, "human_cell_type", "cell_type")

final_matched_celltype_summary <- merge(
  CJ(mouse_model = MODEL_ORDER, cell_type = TARGET_CELLTYPES, unique = TRUE),
  matched_cor,
  by = c("mouse_model", "cell_type"),
  all.x = TRUE,
  sort = FALSE
)
final_matched_celltype_summary <- merge(
  final_matched_celltype_summary,
  matched_sig_target,
  by = c("mouse_model", "cell_type"),
  all.x = TRUE,
  sort = FALSE
)

if ("q_value_BH_within_model" %in% names(final_matched_celltype_summary)) {
  final_matched_celltype_summary[, significance := sig_label(q_value_BH_within_model)]
  final_matched_celltype_summary[, significant_q05 := !is.na(q_value_BH_within_model) & q_value_BH_within_model <= 0.05]
} else {
  final_matched_celltype_summary[, significance := ""]
  final_matched_celltype_summary[, significant_q05 := FALSE]
}

coverage_summary <- final_matched_celltype_summary[, .(
  n_target_celltypes = length(TARGET_CELLTYPES),
  n_available_matched_celltypes = sum(is.finite(matched_rho)),
  n_significant_matched_celltypes = sum(significant_q05, na.rm = TRUE)
), by = mouse_model]

restricted_cor_summary <- cor_long[, {
  matched <- spearman_rho[human_cell_type == mouse_cell_type & is.finite(spearman_rho)]
  mismatched <- spearman_rho[human_cell_type != mouse_cell_type & is.finite(spearman_rho)]
  data.table(
    mean_matched_rho_step3 = mean(matched, na.rm = TRUE),
    median_matched_rho_step3 = stats::median(matched, na.rm = TRUE),
    mean_mismatched_rho_step3 = mean(mismatched, na.rm = TRUE),
    matched_minus_mismatched_step3 = mean(matched, na.rm = TRUE) - mean(mismatched, na.rm = TRUE)
  )
}, by = mouse_model]

restricted_self_ranks <- cor_long[, {
  vals <- spearman_rho
  names(vals) <- mouse_cell_type
  finite_vals <- vals[is.finite(vals)]
  self_value <- vals[human_cell_type[1]]

  if (length(finite_vals) == 0L || length(self_value) == 0L || !is.finite(self_value)) {
    data.table(self_match_available = FALSE, self_rho = NA_real_, self_rank = NA_integer_)
  } else {
    data.table(
      self_match_available = TRUE,
      self_rho = as.numeric(self_value),
      self_rank = as.integer(rank(-finite_vals, ties.method = "min")[human_cell_type[1]])
    )
  }
}, by = .(mouse_model, human_cell_type)]

self_rank_summary <- restricted_self_ranks[, .(
  n_rank1_self_matches_step3 = sum(self_rank == 1L, na.rm = TRUE),
  mean_self_rho_step3 = mean(self_rho, na.rm = TRUE),
  median_self_rho_step3 = stats::median(self_rho, na.rm = TRUE)
), by = mouse_model]

global_cols <- intersect(
  c(
    "mouse_model", "observed_matched_minus_mismatched",
    "p_value_greater", "q_value_BH_greater",
    "p_value_two_sided", "q_value_BH_two_sided"
  ),
  names(global_perm)
)
global_perm_summary <- global_perm[mouse_model %in% MODEL_ORDER, ..global_cols]

pathway_concordance_model_summary <- pathway_summary[
  mouse_model %in% MODEL_ORDER &
    cell_type %in% TARGET_CELLTYPES,
  .(
    mean_pathway_concordant_fraction = mean(concordant_fraction, na.rm = TRUE),
    mean_pathway_product_contribution = mean(mean_product_contribution, na.rm = TRUE),
    mean_pathway_abs_difference = mean(mean_abs_difference, na.rm = TRUE),
    mean_discordant_pathway_count = mean(n_discordant, na.rm = TRUE),
    total_discordant_pathway_count = sum(n_discordant, na.rm = TRUE)
  ),
  by = mouse_model
]

sensitivity_retention_input <- rbind(
  sensitivity_summary[mouse_model %in% MODEL_ORDER, .(mouse_model, analysis_name, mean_matched_rho)],
  leave_one_summary[mouse_model %in% MODEL_ORDER, .(mouse_model, analysis_name, mean_matched_rho)],
  use.names = TRUE
)

full_reference <- sensitivity_retention_input[
  analysis_name == "full_50_pathways",
  .(mouse_model, full_50_mean_matched_rho = mean_matched_rho)
]

sensitivity_retention_summary <- merge(
  sensitivity_retention_input[analysis_name != "full_50_pathways"],
  full_reference,
  by = "mouse_model",
  all.x = TRUE,
  sort = FALSE
)
# METRIC DEFINITION NOTE (documented for Methods):
# Sensitivity retention is the mean of pmin(ratio, 1) over six checks: one
# pathway-set analysis plus five leave-one-cell-type-out analyses, each
# referenced to that model's full_50_pathways value. Ratios above 1 (a
# sensitivity analysis that INCREASES the mean matched correlation) are
# truncated to 1 before averaging, so the reported retention is bounded above
# by 1 by construction and is a conservative rather than symmetric summary.
sensitivity_retention_summary[, retention_ratio := mean_matched_rho / full_50_mean_matched_rho]
sensitivity_retention_summary[, retention_ratio_capped := pmin(retention_ratio, 1)]

sensitivity_retention_model <- sensitivity_retention_summary[, .(
  n_sensitivity_checks = .N,
  mean_sensitivity_retention_ratio = mean(retention_ratio, na.rm = TRUE),
  mean_sensitivity_retention_ratio_capped = mean(retention_ratio_capped, na.rm = TRUE)
), by = mouse_model]

final_sensitivity_summary <- merge(
  sensitivity_summary[
    mouse_model %in% MODEL_ORDER,
    .(
      mouse_model,
      analysis_name,
      n_matched_pairs,
      mean_matched_rho,
      matched_minus_mismatched,
      n_rank1_self_matches,
      rank_by_mean_matched_rho,
      rank_by_matched_advantage
    )
  ],
  rank_stability[mouse_model %in% MODEL_ORDER],
  by = "mouse_model",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("", "_stability")
)

final_model_summary <- Reduce(
  function(x, y) merge(x, y, by = "mouse_model", all = TRUE, sort = FALSE),
  list(
    coverage_summary,
    restricted_cor_summary,
    self_rank_summary,
    global_perm_summary,
    pathway_concordance_model_summary,
    rank_stability,
    sensitivity_retention_model
  )
)

final_model_summary[, rank_mean_matched_rho := frank(-mean_matched_rho_step3, ties.method = "min")]
final_model_summary[, rank_global_structure_p := frank(p_value_greater, ties.method = "min")]
final_model_summary[, rank_sensitivity_retention := frank(-mean_sensitivity_retention_ratio_capped, ties.method = "min")]
final_model_summary[, rank_discordant_pathway_count := frank(mean_discordant_pathway_count, ties.method = "min")]
final_model_summary[, model_order_index := match(mouse_model, MODEL_ORDER)]
setorder(final_model_summary, model_order_index)

protocol_metric_raw <- final_model_summary[, .(
  mouse_model,
  mean_matched_rho = mean_matched_rho_step3,
  rank1_self_match_fraction = n_rank1_self_matches_step3 / pmax(n_available_matched_celltypes, 1),
  global_structure_neg_log10_p = -log10(pmax(p_value_greater, .Machine$double.xmin)),
  sensitivity_retention_ratio = mean_sensitivity_retention_ratio_capped,
  discordant_pathway_count = mean_discordant_pathway_count
)]

protocol_metric_scaled <- copy(protocol_metric_raw)
# SCOPE AND SCALING NOTES (documented for the Fig. 7b legend):
# 1. Dimensions 1, 2, 4 and 5 are computed on the five prespecified
#    compartments. Dimension 3 (global_structure_neg_log10_p) is inherited
#    from Step 4 and was computed on the FULL annotated cell-type space.
#    This scope difference must be stated in the legend.
# 2. Dimension 3 uses the ONE-SIDED p_value_greater and is capped at 4.00 by
#    the number of permutations.
# 3. min-max scaling below is computed ACROSS THE MODELS PRESENT, so on every
#    axis one model is pinned at 1.0 and another at 0.0 regardless of the
#    absolute spread. Raw values must be printed alongside, and the axes must
#    never be summed or averaged into a composite score.
protocol_metric_scaled[, mean_matched_rho := minmax_scale(mean_matched_rho)]
protocol_metric_scaled[, rank1_self_match_fraction := minmax_scale(rank1_self_match_fraction)]
protocol_metric_scaled[, global_structure_neg_log10_p := minmax_scale(global_structure_neg_log10_p)]
protocol_metric_scaled[, sensitivity_retention_ratio := minmax_scale(sensitivity_retention_ratio)]
protocol_metric_scaled[, discordant_pathway_count := minmax_scale_reverse(discordant_pathway_count)]

metric_names <- setdiff(names(protocol_metric_raw), "mouse_model")
metric_labels <- c(
  mean_matched_rho = "Mean Spearman rho",
  rank1_self_match_fraction = "Rank=1 fraction",
  global_structure_neg_log10_p = "-log10 global p",
  sensitivity_retention_ratio = "Sensitivity retention",
  discordant_pathway_count = "Low discordance"
)

protocol_metric_long <- merge(
  melt(protocol_metric_scaled, id.vars = "mouse_model", measure.vars = metric_names, variable.name = "metric", value.name = "scaled_value"),
  melt(protocol_metric_raw, id.vars = "mouse_model", measure.vars = metric_names, variable.name = "metric", value.name = "raw_value"),
  by = c("mouse_model", "metric"),
  all = TRUE,
  sort = FALSE
)
protocol_metric_long[, metric_label := factor(metric_labels[metric], levels = rev(metric_labels[metric_names]))]
protocol_metric_long[, mouse_model := factor(mouse_model, levels = MODEL_ORDER)]
protocol_metric_long[, raw_value_label := raw_label(metric, raw_value)]
protocol_metric_long[, label_x := pmin(1.02, scaled_value + 0.035)]

figure_y_b <- make_protocol_dimension_radar(
  dt = protocol_metric_long,
  model_levels = MODEL_ORDER,
  model_colors = MODEL_COLORS,
  title = "Y-B. Protocol-dimension radar",
  subtitle = "Axes use min-max scaled values; raw values are exported in the companion table"
)

figure_y_b_dotplot <- ggplot(protocol_metric_long, aes(x = scaled_value, y = metric_label, color = mouse_model)) +
  geom_vline(xintercept = c(0, 0.5, 1), color = "grey88", linewidth = 0.35) +
  geom_line(aes(group = metric_label), color = "grey82", linewidth = 0.3) +
  geom_point(size = 3, na.rm = TRUE) +
  geom_text(aes(x = label_x, label = raw_value_label), size = 2.8, hjust = 0, show.legend = FALSE, na.rm = TRUE) +
  facet_wrap(~ mouse_model, nrow = 1) +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE) +
  scale_x_continuous(limits = c(0, 1.18), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(
    title = "Y-B. Protocol-dimension dot plot",
    subtitle = "Point position is min-max scaled within each dimension; labels show raw values",
    x = "Scaled score within dimension",
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = "grey65"),
    axis.text.y = element_text(color = "grey20"),
    plot.title = element_text(hjust = 0)
  )

figure_y_b_single_radars <- lapply(MODEL_ORDER, function(model_name) {
  make_single_protocol_dimension_radar(
    dt = protocol_metric_long,
    model_name = model_name,
    model_color = MODEL_COLORS[[model_name]],
    title = model_name
  )
})
names(figure_y_b_single_radars) <- MODEL_ORDER

global_null_plot_dt <- merge(
  global_null[mouse_model %in% MODEL_ORDER],
  global_perm_summary[, .(mouse_model, observed_matched_minus_mismatched, p_value_greater)],
  by = "mouse_model",
  all.x = TRUE,
  sort = FALSE
)
global_null_plot_dt[, mouse_model := factor(mouse_model, levels = MODEL_ORDER)]

figure_y_c <- ggplot(global_null_plot_dt, aes(x = null_matched_minus_mismatched)) +
  geom_histogram(bins = 45, fill = "grey78", color = "white", linewidth = 0.15) +
  geom_vline(aes(xintercept = observed_matched_minus_mismatched, color = mouse_model), linewidth = 0.85) +
  facet_wrap(~ mouse_model, nrow = 2, scales = "free_y") +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE) +
  labs(
    title = "Y-C. Global cell-label permutation null distributions",
    subtitle = "Vertical line marks observed mean matched rho - mean mismatched rho",
    x = "Null matched-minus-mismatched statistic",
    y = "Permutation count"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = "grey65"),
    plot.title = element_text(hjust = 0)
  )

sensitivity_heatmap_table <- rbind(
  sensitivity_summary[mouse_model %in% MODEL_ORDER, .(mouse_model, analysis_name, mean_matched_rho)],
  leave_one_summary[mouse_model %in% MODEL_ORDER, .(mouse_model, analysis_name, mean_matched_rho)],
  use.names = TRUE
)
sensitivity_levels <- c(
  "full_50_pathways",
  "without_broad_pathways",
  paste0("leave_out_", TARGET_CELLTYPES)
)
sensitivity_heatmap_table[, analysis_name := factor(analysis_name, levels = sensitivity_levels)]
sensitivity_heatmap_table[, mouse_model := factor(mouse_model, levels = rev(MODEL_ORDER))]
sensitivity_heatmap_table[, value_label := ifelse(is.finite(mean_matched_rho), sprintf("%.2f", mean_matched_rho), "NA")]

figure_y_d <- ggplot(sensitivity_heatmap_table, aes(x = analysis_name, y = mouse_model, fill = mean_matched_rho)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = value_label), size = 2.8) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey85",
    name = "Mean matched rho",
    guide = guide_colourbar(raster = FALSE)
  ) +
  labs(
    title = "Y-D. Sensitivity analysis",
    subtitle = "Mean matched Spearman rho across pathway and cell-type sensitivity conditions",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0)
  )

pathway_driver_target <- pathway_driver[
  mouse_model %in% MODEL_ORDER &
    cell_type %in% TARGET_CELLTYPES &
    matched_cell_type == TRUE
]

top_pathways_by_class <- pathway_driver_target[
  direction_class %in% c("Concordantly UP", "Concordantly DOWN", "Discordant"),
  .(
    mean_product_contribution = mean(product_contribution, na.rm = TRUE),
    mean_abs_product_contribution = mean(abs(product_contribution), na.rm = TRUE)
  ),
  by = .(direction_class, pathway)
][
  order(direction_class, -mean_abs_product_contribution),
  head(.SD, 8),
  by = direction_class
]

pathway_heatmap_table <- merge(
  pathway_driver_target,
  top_pathways_by_class[, .(direction_class, pathway)],
  by = c("direction_class", "pathway"),
  all = FALSE,
  sort = FALSE
)
pathway_heatmap_table[, model_celltype := paste(mouse_model, cell_type, sep = " | ")]
pathway_heatmap_table <- pathway_heatmap_table[, .(
  product_contribution = mean(product_contribution, na.rm = TRUE)
), by = .(direction_class, pathway, mouse_model, cell_type, model_celltype)]

pathway_order <- top_pathways_by_class[
  order(
    match(direction_class, c("Concordantly UP", "Concordantly DOWN", "Discordant")),
    -mean_abs_product_contribution
  ),
  paste(direction_class, pathway, sep = " :: ")
]
pathway_heatmap_table[, pathway_panel_label := paste(direction_class, pathway, sep = " :: ")]
pathway_heatmap_table[, pathway_panel_label := factor(pathway_panel_label, levels = rev(pathway_order))]

model_celltype_levels <- as.vector(outer(MODEL_ORDER, TARGET_CELLTYPES, paste, sep = " | "))
pathway_heatmap_table[, model_celltype := factor(model_celltype, levels = model_celltype_levels)]
pathway_heatmap_table[, product_label := ifelse(is.finite(product_contribution), sprintf("%.1f", product_contribution), "")]

figure_y_e <- ggplot(pathway_heatmap_table, aes(x = model_celltype, y = pathway_panel_label, fill = product_contribution)) +
  geom_tile(color = "white", linewidth = 0.18) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    na.value = "grey85",
    name = "Product contribution",
    guide = guide_colourbar(raster = FALSE)
  ) +
  labs(
    title = "Y-E. Pathway concordance and discordance drivers",
    subtitle = "Rows show top matched-pair pathways within concordant UP, concordant DOWN, and discordant classes",
    x = "Model | cell type",
    y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(hjust = 0)
  )

make_clean_radar <- function(dt, model_col, axis_col, value_col,
                             axis_levels, model_levels, model_colors,
                             title, subtitle,
                             ring_breaks = c(0, 0.25, 0.5, 0.75, 1),
                             ring_labels = c("0", "0.25", "0.50", "0.75", "1.00"),
                             value_min = 0, value_max = 1) {
  d <- copy(dt)
  setnames(d, c(model_col, axis_col, value_col), c("model_id_plot", "axis_id_plot", "value_plot"))
  d[, model_id_plot := factor(as.character(model_id_plot), levels = model_levels)]
  d[, axis_id_plot := factor(as.character(axis_id_plot), levels = axis_levels)]

  axis_meta <- data.table(
    axis_id_plot = factor(axis_levels, levels = axis_levels),
    axis_index = seq_along(axis_levels)
  )
  axis_meta[, angle := pi / 2 - 2 * pi * (axis_index - 1) / length(axis_levels)]
  axis_meta[, x_end := cos(angle)]
  axis_meta[, y_end := sin(angle)]
  axis_meta[, x_label := 1.20 * cos(angle)]
  axis_meta[, y_label := 1.20 * sin(angle)]
  axis_meta[, hjust := fifelse(x_label > 0.1, 0, fifelse(x_label < -0.1, 1, 0.5))]
  axis_meta[, vjust := fifelse(y_label > 0.1, 0, fifelse(y_label < -0.1, 1, 0.5))]

  d <- merge(d, axis_meta[, .(axis_id_plot, axis_index, angle)], by = "axis_id_plot", all.x = TRUE, sort = FALSE)
  d[, radius := (value_plot - value_min) / (value_max - value_min)]
  d[!is.finite(radius), radius := NA_real_]
  d[, radius := pmax(0, pmin(1, radius))]
  d[, x := radius * cos(angle)]
  d[, y := radius * sin(angle)]
  setorder(d, model_id_plot, axis_index)

  closed <- d[is.finite(x) & is.finite(y), {
    if (.N < 3L) .SD else rbind(.SD, .SD[1])
  }, by = model_id_plot]

  grid <- rbindlist(lapply(ring_breaks, function(r) {
    theta <- seq(0, 2 * pi, length.out = 241)
    data.table(ring = r, x = r * cos(theta), y = r * sin(theta))
  }))

  ring_label_dt <- data.table(
    ring = ring_breaks,
    label = ring_labels,
    x = ring_breaks * cos(pi / 2),
    y = ring_breaks * sin(pi / 2)
  )[ring > 0]

  ggplot() +
    geom_path(data = grid, aes(x = x, y = y, group = ring), color = "grey84", linewidth = 0.3) +
    geom_segment(data = axis_meta, aes(x = 0, y = 0, xend = x_end, yend = y_end), color = "grey84", linewidth = 0.3) +
    geom_text(data = axis_meta, aes(x = x_label, y = y_label, label = axis_id_plot, hjust = hjust, vjust = vjust), size = 3.3, color = "grey25") +
    geom_text(data = ring_label_dt, aes(x = x, y = y, label = label), nudge_x = 0.04, size = 2.4, color = "grey45") +
    geom_polygon(data = closed, aes(x = x, y = y, group = model_id_plot, color = model_id_plot, fill = model_id_plot), alpha = 0.09, linewidth = 0.75, na.rm = TRUE) +
    geom_path(data = closed, aes(x = x, y = y, group = model_id_plot, color = model_id_plot), linewidth = 0.9, na.rm = TRUE) +
    geom_point(data = d, aes(x = x, y = y, color = model_id_plot), size = 2.2, na.rm = TRUE) +
    coord_equal(xlim = c(-1.45, 1.55), ylim = c(-1.35, 1.35), clip = "off") +
    scale_color_manual(values = model_colors, drop = FALSE) +
    scale_fill_manual(values = model_colors, drop = FALSE) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_void(base_size = 10) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0, size = 14),
      plot.subtitle = element_text(hjust = 0, size = 9),
      plot.margin = margin(12, 36, 12, 36)
    )
}

celltype_radar_metric_table <- final_matched_celltype_summary[, .(
  mouse_model,
  cell_type,
  matched_rho
)]
celltype_radar_metric_table[, mouse_model := factor(mouse_model, levels = MODEL_ORDER)]
celltype_radar_metric_table[, cell_type := factor(cell_type, levels = TARGET_CELLTYPES)]
celltype_radar_metric_table <- celltype_radar_metric_table[order(mouse_model, cell_type)]

figure_y_a_generated <- make_clean_radar(
  dt = celltype_radar_metric_table,
  model_col = "mouse_model",
  axis_col = "cell_type",
  value_col = "matched_rho",
  axis_levels = TARGET_CELLTYPES,
  model_levels = MODEL_ORDER,
  model_colors = MODEL_COLORS,
  title = "Y-A. Matched cell-type Spearman rho",
  subtitle = "Axis values are original rho; radial mapping uses rho from -1 to 1",
  ring_breaks = c(0, 0.25, 0.5, 0.75, 1),
  ring_labels = c("-1.00", "-0.50", "0", "0.50", "1.00"),
  value_min = -1,
  value_max = 1
)

existing_radar_source <- file.path(STEP8_FIG_DIR, "F6_celltype_matched_rho_radar.pdf")
existing_radar_target <- copy_existing_pdf(existing_radar_source, "FigureY_A_celltype_matched_rho_radar.pdf")
if (!file.exists(existing_radar_target)) {
  save_vector_pdf(figure_y_a_generated, "FigureY_A_celltype_matched_rho_radar", width = 8.2, height = 6.2)
}

copy_existing_pdf(file.path(STEP5_FIG_DIR, "F1_all_models_correlation_heatmap.pdf"), "FigureX_C_existing_correlation_heatmap.pdf")
copy_existing_pdf(file.path(STEP5_FIG_DIR, "F2_matched_celltype_correlation_dotplot.pdf"), "FigureX_D_existing_matched_dotplot.pdf")
copy_existing_pdf(file.path(STEP5_FIG_DIR, "F4_self_match_rank_heatmap.pdf"), "FigureX_E_existing_self_match_rank_heatmap.pdf")

save_vector_pdf(figure_y_b, "FigureY_B_protocol_dimension_radar", width = 8.2, height = 6.2)
save_vector_pdf(figure_y_b_dotplot, "FigureY_B_protocol_dimension_dotplot", width = 11.5, height = 4.8)
individual_radar_files <- character()
for (model_name in MODEL_ORDER) {
  model_file_id <- gsub("[^A-Za-z0-9_.-]+", "_", model_name)
  output_name <- paste0("FigureY_B_protocol_dimension_radar_", model_file_id)
  save_vector_pdf(figure_y_b_single_radars[[model_name]], output_name, width = 6.4, height = 4.6)
  individual_radar_files <- c(individual_radar_files, file.path(FIG_DIR, paste0(output_name, ".pdf")))
}
grDevices::pdf(
  file = file.path(FIG_DIR, "FigureY_B_protocol_dimension_individual_radars.pdf"),
  width = 12.8,
  height = 9.2,
  onefile = TRUE,
  useDingbats = FALSE,
  colormodel = "rgb",
  bg = "white"
)
grid::grid.newpage()
single_layout <- grid::grid.layout(nrow = 2, ncol = 2)
grid::pushViewport(grid::viewport(layout = single_layout))
for (i in seq_along(MODEL_ORDER)) {
  print(
    figure_y_b_single_radars[[MODEL_ORDER[i]]],
    vp = grid::viewport(
      layout.pos.row = ceiling(i / 2),
      layout.pos.col = ifelse(i %% 2 == 1, 1, 2)
    )
  )
}
grid::popViewport()
grDevices::dev.off()
save_vector_pdf(figure_y_c, "FigureY_C_global_permutation_null_distributions", width = 8.5, height = 6.5)
save_vector_pdf(figure_y_d, "FigureY_D_sensitivity_analysis_heatmap", width = 8.8, height = 4.8)
save_vector_pdf(figure_y_e, "FigureY_E_pathway_concordance_discordance_heatmap", width = 12.5, height = 8.5)

grDevices::pdf(
  file = file.path(FIG_DIR, "FigureY_composite_A_to_E.pdf"),
  width = 15,
  height = 12,
  onefile = TRUE,
  useDingbats = FALSE,
  colormodel = "rgb",
  bg = "white"
)
grid::grid.newpage()
layout <- grid::grid.layout(
  nrow = 3,
  ncol = 2,
  widths = grid::unit(c(1.15, 1), "null"),
  heights = grid::unit(c(1, 1, 1.1), "null")
)
grid::pushViewport(grid::viewport(layout = layout))
print(figure_y_a_generated, vp = grid::viewport(layout.pos.row = 1:2, layout.pos.col = 1))
print(figure_y_b, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
print(figure_y_c, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
print(figure_y_d, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
print(figure_y_e, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 2))
grid::popViewport()
grDevices::dev.off()

figure_caption_statistics <- data.table(
  item = c(
    "target_celltypes",
    "top_model_by_mean_matched_rho",
    "top_model_by_global_structure_p",
    "top_model_by_sensitivity_retention",
    "lowest_model_by_discordant_pathway_count",
    "note_STIA_coverage"
  ),
  value = c(
    paste(TARGET_CELLTYPES, collapse = ", "),
    final_model_summary[order(-mean_matched_rho_step3), mouse_model][1],
    final_model_summary[order(p_value_greater), mouse_model][1],
    final_model_summary[order(-mean_sensitivity_retention_ratio_capped), mouse_model][1],
    final_model_summary[order(mean_discordant_pathway_count), mouse_model][1],
    paste0(
      "STIA available matched cell types: ",
      final_model_summary[mouse_model == "STIA", n_available_matched_celltypes],
      " of ", length(TARGET_CELLTYPES)
    )
  )
)

panel_manifest <- data.table(
  figure = c(
    "Figure X",
    "Figure X",
    "Figure X",
    "Figure Y",
    "Figure Y",
    "Figure Y",
    "Figure Y",
    "Figure Y",
    "Figure Y",
    "Figure Y",
    "Figure Y"
  ),
  panel = c("X-C", "X-D", "X-E", "Y-A", "Y-B radar", "Y-B dotplot", "Y-B individual radars", "Y-C", "Y-D", "Y-E", "Y-composite"),
  role = c(
    "Existing correlation heatmap evidence",
    "Existing matched cell-type dot plot evidence",
    "Existing self-match rank heatmap evidence",
    "Existing radar panel copied from Step 8 when available",
    "New protocol-dimension total radar plot",
    "New protocol-dimension dot plot",
    "New protocol-dimension single-model radar composite",
    "New permutation null distribution panel",
    "New sensitivity heatmap panel",
    "New pathway concordance/discordance heatmap panel",
    "Generated composite layout for Figure Y panels A-E"
  ),
  file = file.path(FIG_DIR, c(
    "FigureX_C_existing_correlation_heatmap.pdf",
    "FigureX_D_existing_matched_dotplot.pdf",
    "FigureX_E_existing_self_match_rank_heatmap.pdf",
    "FigureY_A_celltype_matched_rho_radar.pdf",
    "FigureY_B_protocol_dimension_radar.pdf",
    "FigureY_B_protocol_dimension_dotplot.pdf",
    "FigureY_B_protocol_dimension_individual_radars.pdf",
    "FigureY_C_global_permutation_null_distributions.pdf",
    "FigureY_D_sensitivity_analysis_heatmap.pdf",
    "FigureY_E_pathway_concordance_discordance_heatmap.pdf",
    "FigureY_composite_A_to_E.pdf"
  ))
)
panel_manifest <- rbind(
  panel_manifest,
  data.table(
    figure = "Figure Y",
    panel = paste0("Y-B individual radar: ", MODEL_ORDER),
    role = "New protocol-dimension single-model radar using the matched model color",
    file = individual_radar_files
  ),
  use.names = TRUE
)
panel_manifest[, file_exists := file.exists(file)]

fwrite(final_model_summary, file.path(TABLE_DIR, "final_model_summary.csv"))
fwrite(final_matched_celltype_summary, file.path(TABLE_DIR, "final_matched_celltype_summary.csv"))
fwrite(pathway_concordance_model_summary, file.path(TABLE_DIR, "final_pathway_driver_summary.csv"))
fwrite(model_pathway_summary[mouse_model %in% MODEL_ORDER], file.path(TABLE_DIR, "model_pathway_summary.csv"))
fwrite(final_sensitivity_summary, file.path(TABLE_DIR, "final_sensitivity_summary.csv"))
fwrite(figure_caption_statistics, file.path(TABLE_DIR, "figure_caption_statistics.csv"))
fwrite(protocol_metric_raw, file.path(TABLE_DIR, "figure_y_protocol_dimension_table_raw.csv"))
fwrite(protocol_metric_scaled, file.path(TABLE_DIR, "figure_y_protocol_dimension_table_scaled.csv"))
fwrite(protocol_metric_long, file.path(TABLE_DIR, "figure_y_protocol_dimension_radar_table.csv"))
fwrite(global_null_plot_dt, file.path(TABLE_DIR, "figure_y_global_null_distribution_table.csv"))
fwrite(sensitivity_heatmap_table, file.path(TABLE_DIR, "figure_y_sensitivity_heatmap_table.csv"))
fwrite(pathway_heatmap_table, file.path(TABLE_DIR, "figure_y_pathway_heatmap_table.csv"))
fwrite(celltype_radar_metric_table, file.path(TABLE_DIR, "figure_y_celltype_radar_metric_table.csv"))
fwrite(panel_manifest, file.path(TABLE_DIR, "figure_y_panel_manifest.csv"))

fwrite(final_model_summary, file.path(OUT_DIR, "final_model_summary.csv"))
fwrite(final_matched_celltype_summary, file.path(OUT_DIR, "final_matched_celltype_summary.csv"))
fwrite(pathway_concordance_model_summary, file.path(OUT_DIR, "final_pathway_driver_summary.csv"))
fwrite(model_pathway_summary[mouse_model %in% MODEL_ORDER], file.path(OUT_DIR, "model_pathway_summary.csv"))
fwrite(final_sensitivity_summary, file.path(OUT_DIR, "final_sensitivity_summary.csv"))
fwrite(figure_caption_statistics, file.path(OUT_DIR, "figure_caption_statistics.csv"))
fwrite(protocol_metric_raw, file.path(OUT_DIR, "protocol_radar_metric_table_raw.csv"))
fwrite(protocol_metric_scaled, file.path(OUT_DIR, "protocol_radar_metric_table_scaled.csv"))
fwrite(protocol_metric_long, file.path(OUT_DIR, "protocol_dimension_radar_table.csv"))
fwrite(celltype_radar_metric_table, file.path(OUT_DIR, "celltype_radar_metric_table.csv"))

saveRDS(
  list(
    final_model_summary = final_model_summary,
    final_matched_celltype_summary = final_matched_celltype_summary,
    pathway_concordance_model_summary = pathway_concordance_model_summary,
    final_sensitivity_summary = final_sensitivity_summary,
    figure_caption_statistics = figure_caption_statistics,
    protocol_metric_raw = protocol_metric_raw,
    protocol_metric_scaled = protocol_metric_scaled,
    protocol_metric_long = protocol_metric_long,
    global_null_plot_dt = global_null_plot_dt,
    sensitivity_heatmap_table = sensitivity_heatmap_table,
    pathway_heatmap_table = pathway_heatmap_table,
    celltype_radar_metric_table = celltype_radar_metric_table,
    panel_manifest = panel_manifest,
    target_celltypes = TARGET_CELLTYPES
  ),
  file.path(OUT_DIR, "step8_final_summary_and_figure_y_results.rds")
)

message("Step 8 complete.")
message("Output directory: ", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE))
message("Figures:")
message("  - figures/FigureY_A_celltype_matched_rho_radar.pdf")
message("  - figures/FigureY_B_protocol_dimension_radar.pdf")
message("  - figures/FigureY_B_protocol_dimension_dotplot.pdf")
message("  - figures/FigureY_B_protocol_dimension_radar_[model].pdf")
message("  - figures/FigureY_B_protocol_dimension_individual_radars.pdf")
message("  - figures/FigureY_C_global_permutation_null_distributions.pdf")
message("  - figures/FigureY_D_sensitivity_analysis_heatmap.pdf")
message("  - figures/FigureY_E_pathway_concordance_discordance_heatmap.pdf")
message("  - figures/FigureY_composite_A_to_E.pdf")
message("Tables:")
message("  - tables/final_model_summary.csv")
message("  - tables/final_matched_celltype_summary.csv")
message("  - tables/figure_y_panel_manifest.csv")
message("  - tables/figure_y_protocol_dimension_radar_table.csv")
message("  - tables/figure_y_global_null_distribution_table.csv")
message("  - tables/figure_y_sensitivity_heatmap_table.csv")
message("  - tables/figure_y_pathway_heatmap_table.csv")
