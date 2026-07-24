options(stringsAsFactors = FALSE)

# Step 4 only
# Purpose:
# 1. Pairwise pathway-label permutation:
#    For each hsRA cell_type x mouse cell_type pair, test whether the observed
#    Spearman correlation is stronger than random Hallmark pathway alignment.
#
# 2. Model-level matched-vs-mismatched permutation:
#    For each mouse model, test whether same-name cell-type pairs have higher
#    correlations than non-matched cell-type pairs.
#
# This step reads Step 2 and Step 3 outputs.
# It does not make final figures.

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

set.seed(20260623L)

N_PERM_PAIRWISE <- 1000L
N_PERM_GLOBAL <- 10000L
MIN_PATHWAYS_FOR_COR <- 10L

# SCOPE AND PARAMETER NOTES (documented for Methods):
# 1. Pairwise p values are TWO-SIDED: abs(null) >= abs(observed).
# 2. Empirical p uses the (1 + n_extreme) / (n_perm + 1) correction, so the
#    smallest attainable pairwise p is 1/1001 and the smallest attainable
#    global p is 1/10001. -log10(global p) is therefore capped at exactly
#    4.0000; report such values as P < 1e-4, not as a point estimate.
# 3. Two BH families are exported here (all pairs / within model). Step 5
#    computes a THIRD family over the five matched pairs within each model,
#    and that third family supplies the asterisks in the published dot plot.

FIG6_DIR <- "E:/OneDrive/01.Academic/Papers/In Writing/C&AIA/Manuscript/8nd Nauture Communications/Figure 6"
STEP2_DIR <- file.path(FIG6_DIR, "RA-Scores", "02_step2_make_profiles")
STEP3_DIR <- file.path(FIG6_DIR, "RA-Scores", "03_step3_correlation_matrices")
OUT_DIR <- file.path(FIG6_DIR, "RA-Scores", "04_step4_permutation_tests")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

profile_file <- file.path(STEP2_DIR, "zscore_model_profiles_wide.csv")
correlation_file <- file.path(STEP3_DIR, "all_hsRA_vs_mouse_spearman_long.csv")

if (!file.exists(profile_file)) {
  stop("Cannot find Step 2 z-score profile file:\n", profile_file, call. = FALSE)
}
if (!file.exists(correlation_file)) {
  stop("Cannot find Step 3 correlation file:\n", correlation_file, call. = FALSE)
}

profiles <- fread(profile_file, showProgress = FALSE)
correlation_long <- fread(correlation_file, showProgress = FALSE)

pathways <- grep("^HALLMARK_", names(profiles), value = TRUE)
if (length(pathways) != 50L) {
  stop("Expected exactly 50 HALLMARK_ columns; found ", length(pathways), call. = FALSE)
}

required_cor_cols <- c(
  "mouse_model", "human_cell_type", "mouse_cell_type",
  "matched_cell_type", "spearman_rho"
)
missing_cor_cols <- setdiff(required_cor_cols, names(correlation_long))
if (length(missing_cor_cols) > 0) {
  stop(
    "Step 3 correlation file is missing column(s): ",
    paste(missing_cor_cols, collapse = ", "),
    call. = FALSE
  )
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
      call. = FALSE
    )
  }
  mat <- as.matrix(x[, ..pathways])
  rownames(mat) <- x$cell_type
  mat
}

spearman_one <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < MIN_PATHWAYS_FOR_COR) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}

pairwise_permutation_one <- function(h, m, observed, n_perm) {
  ok <- is.finite(h) & is.finite(m)
  h <- h[ok]
  m <- m[ok]

  if (length(h) < MIN_PATHWAYS_FOR_COR || !is.finite(observed)) {
    return(list(
      n_pathways_used = length(h),
      p_two_sided = NA_real_,
      null_mean = NA_real_,
      null_sd = NA_real_
    ))
  }

  null <- replicate(
    n_perm,
    suppressWarnings(stats::cor(h, sample(m), method = "spearman"))
  )

  list(
    n_pathways_used = length(h),
    p_two_sided = (1 + sum(abs(null) >= abs(observed), na.rm = TRUE)) / (n_perm + 1),
    null_mean = mean(null, na.rm = TRUE),
    null_sd = stats::sd(null, na.rm = TRUE)
  )
}

human_mat <- make_profile_matrix(profiles, "hsRA", pathways)
mouse_models <- sort(unique(correlation_long$mouse_model))

pairwise_results <- vector("list", nrow(correlation_long))

for (i in seq_len(nrow(correlation_long))) {
  row_i <- correlation_long[i]
  mouse_mat <- make_profile_matrix(profiles, row_i$mouse_model, pathways)

  h <- as.numeric(human_mat[row_i$human_cell_type, pathways])
  m <- as.numeric(mouse_mat[row_i$mouse_cell_type, pathways])

  test <- pairwise_permutation_one(
    h = h,
    m = m,
    observed = row_i$spearman_rho,
    n_perm = N_PERM_PAIRWISE
  )

  pairwise_results[[i]] <- data.table(
    mouse_model = row_i$mouse_model,
    human_cell_type = row_i$human_cell_type,
    mouse_cell_type = row_i$mouse_cell_type,
    matched_cell_type = row_i$matched_cell_type,
    observed_spearman_rho = row_i$spearman_rho,
    n_pathways_used = test$n_pathways_used,
    permutation_type = "pathway-label permutation within cell-type pair",
    n_permutations = N_PERM_PAIRWISE,
    p_value_two_sided = test$p_two_sided,
    null_mean = test$null_mean,
    null_sd = test$null_sd
  )
}

pairwise_permutation_results <- rbindlist(pairwise_results, use.names = TRUE)
pairwise_permutation_results[, q_value_BH_all_pairs := stats::p.adjust(p_value_two_sided, method = "BH")]
pairwise_permutation_results[, q_value_BH_within_model := stats::p.adjust(p_value_two_sided, method = "BH"),
  by = mouse_model
]
setorder(pairwise_permutation_results, mouse_model, human_cell_type, mouse_cell_type)

matched_celltype_significance <- pairwise_permutation_results[matched_cell_type == TRUE]
setorder(matched_celltype_significance, mouse_model, q_value_BH_within_model, human_cell_type)

global_stat_from_long <- function(d) {
  matched <- d[matched_cell_type == TRUE & is.finite(spearman_rho), spearman_rho]
  mismatched <- d[matched_cell_type == FALSE & is.finite(spearman_rho), spearman_rho]

  if (length(matched) == 0L || length(mismatched) == 0L) {
    return(NA_real_)
  }

  mean(matched) - mean(mismatched)
}

global_results <- vector("list", length(mouse_models))
global_null_list <- vector("list", length(mouse_models))
names(global_results) <- mouse_models
names(global_null_list) <- mouse_models

# SCOPE NOTE (documented for Methods):
# The global cell-label permutation below operates on the FULL annotated
# cell-type space of each model (mouse_types = all unique mouse_cell_type in
# the correlation matrix), not on the five benchmarked compartments. This is
# more conservative than a permutation restricted to the matched space,
# because the label pool is larger. Note that dimensions 1, 2, 4 and 5 of the
# Step 8 summary ARE restricted to the five compartments, so this dimension
# has a different analytical scope and must be described as such.
for (m in mouse_models) {
  d <- correlation_long[mouse_model == m]
  observed <- global_stat_from_long(d)

  human_types <- sort(unique(d$human_cell_type))
  mouse_types <- sort(unique(d$mouse_cell_type))

  null_stats <- numeric(N_PERM_GLOBAL)

  for (b in seq_len(N_PERM_GLOBAL)) {
    perm_map <- data.table(
      mouse_cell_type = mouse_types,
      permuted_mouse_cell_type = sample(mouse_types, length(mouse_types), replace = FALSE)
    )

    perm_d <- merge(d, perm_map, by = "mouse_cell_type", all.x = TRUE, sort = FALSE)
    perm_d[, matched_cell_type := human_cell_type == permuted_mouse_cell_type]
    null_stats[b] <- global_stat_from_long(perm_d)
  }

  global_results[[m]] <- data.table(
    mouse_model = m,
    observed_matched_minus_mismatched = observed,
    n_human_celltypes = length(human_types),
    n_mouse_celltypes = length(mouse_types),
    n_permutations = N_PERM_GLOBAL,
    permutation_type = "mouse cell-type label permutation within model correlation matrix",
    p_value_greater = (1 + sum(null_stats >= observed, na.rm = TRUE)) / (N_PERM_GLOBAL + 1),
    p_value_two_sided = (1 + sum(abs(null_stats) >= abs(observed), na.rm = TRUE)) / (N_PERM_GLOBAL + 1),
    null_mean = mean(null_stats, na.rm = TRUE),
    null_sd = stats::sd(null_stats, na.rm = TRUE),
    null_q025 = as.numeric(stats::quantile(null_stats, 0.025, na.rm = TRUE)),
    null_q975 = as.numeric(stats::quantile(null_stats, 0.975, na.rm = TRUE))
  )

  global_null_list[[m]] <- data.table(
    mouse_model = m,
    iteration = seq_len(N_PERM_GLOBAL),
    null_matched_minus_mismatched = null_stats
  )
}

global_permutation_results <- rbindlist(global_results, use.names = TRUE)
global_permutation_results[, q_value_BH_greater := stats::p.adjust(p_value_greater, method = "BH")]
global_permutation_results[, q_value_BH_two_sided := stats::p.adjust(p_value_two_sided, method = "BH")]
setorder(global_permutation_results, p_value_greater, mouse_model)

global_null_distribution <- rbindlist(global_null_list, use.names = TRUE)

fwrite(pairwise_permutation_results, file.path(OUT_DIR, "pairwise_pathway_permutation_results.csv"))
fwrite(matched_celltype_significance, file.path(OUT_DIR, "matched_celltype_significance.csv"))
fwrite(global_permutation_results, file.path(OUT_DIR, "global_matched_vs_mismatched_permutation.csv"))
fwrite(global_null_distribution, file.path(OUT_DIR, "global_null_distribution.csv"))

saveRDS(
  list(
    pairwise_pathway_permutation_results = pairwise_permutation_results,
    matched_celltype_significance = matched_celltype_significance,
    global_matched_vs_mismatched_permutation = global_permutation_results,
    global_null_distribution = global_null_distribution,
    n_perm_pairwise = N_PERM_PAIRWISE,
    n_perm_global = N_PERM_GLOBAL,
    pairwise_rule = "shuffle mouse Hallmark pathway labels within each hsRA-mouse cell-type pair",
    global_rule = "permute mouse cell-type labels within each model correlation matrix"
  ),
  file.path(OUT_DIR, "step4_permutation_results.rds")
)

message("Step 4 complete.")
message("Step 2 input: ", normalizePath(profile_file, winslash = "/", mustWork = FALSE))
message("Step 3 input: ", normalizePath(correlation_file, winslash = "/", mustWork = FALSE))
message("Output directory: ", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE))
message("Pairwise permutations per pair: ", N_PERM_PAIRWISE)
message("Global permutations per model: ", N_PERM_GLOBAL)
message("Check these files:")
message("  - pairwise_pathway_permutation_results.csv")
message("  - matched_celltype_significance.csv")
message("  - global_matched_vs_mismatched_permutation.csv")
message("  - global_null_distribution.csv")
# Step 5 is chained here so that Steps 4 and 5 run as one unit, which is how
# all reported results were generated. To run the steps independently, set
#   RUN_STEP5_AFTER <- FALSE
# in the global environment before sourcing this script.
if (!exists("RUN_STEP5_AFTER") || isTRUE(RUN_STEP5_AFTER)) {
  source(file.path(FIG6_DIR, "RA-Scores", "05_step5_make_figures.R"))
}
