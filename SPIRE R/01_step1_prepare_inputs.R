options(stringsAsFactors = FALSE)

# Step 1 only
# 1. AUCell tables have already been cell-filtered upstream.
#    This script does not remove any cell subgroup after AUCell scoring.
# 2. AUCell means are calculated with equal sample weight:
#    cell-level scores -> sample-level mean -> equal-weight mean across samples.
# 3. GSE254560 uses sample == "PBS".

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
DATA_DIR <- file.path(FIG6_DIR, "AUCell", "\u7eb3\u5165\u7edf\u8ba1")
OUT_DIR <- file.path(FIG6_DIR, "RA-Scores", "01_step1_prepare_inputs")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

find_one_file <- function(pattern, data_dir = DATA_DIR) {
  if (!dir.exists(data_dir)) {
    stop(
      "DATA_DIR does not exist:\n", data_dir,
      "\nPlease check the DATA_DIR path in this script.",
      call. = FALSE
    )
  }

  all_csv <- list.files(
    data_dir,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  hits <- all_csv[grepl(pattern, basename(all_csv), ignore.case = TRUE, perl = TRUE)]

  if (length(hits) != 1L) {
    stop(
      "Expected exactly one file for pattern: ", pattern,
      "\nFound ", length(hits), " file(s):\n",
      if (length(hits) > 0) paste(basename(hits), collapse = "\n") else "(none)",
      "\n\nDATA_DIR used by this script:\n",
      data_dir,
      "\n\nCSV files visible under DATA_DIR:\n",
      if (length(all_csv) > 0) paste(basename(all_csv), collapse = "\n") else "(none)",
      call. = FALSE
    )
  }
  hits
}

extract_regex <- function(x, pattern, replacement = "\\1") {
  hit <- grepl(pattern, x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  out[hit] <- sub(pattern, replacement, x[hit], perl = TRUE)
  out
}

assert_cols <- function(d, cols, file_label) {
  missing <- setdiff(cols, names(d))
  if (length(missing) > 0) {
    stop(
      file_label, " missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

hallmark_cols <- function(d, file_label) {
  cols <- grep("^HALLMARK_", names(d), value = TRUE)
  if (length(cols) != 50L) {
    stop(
      file_label, " should contain exactly 50 HALLMARK_ columns; found ",
      length(cols),
      call. = FALSE
    )
  }
  cols
}

standardize_cell_type <- function(raw_type, celltype_map = NULL) {
  raw_type <- as.character(raw_type)
  mapped <- rep(NA_character_, length(raw_type))

  if (!is.null(celltype_map)) {
    mapped <- unname(celltype_map[raw_type])
  }

  # Keep unmapped cell types with original labels.
  mapped[is.na(mapped) | mapped == ""] <- raw_type[is.na(mapped) | mapped == ""]
  mapped
}

standardize_scores <- function(path, model_id, source_id, celltype_col, sample_col = NULL,
                               celltype_map = NULL, keep_fun = NULL, sample_fun = NULL) {
  file_label <- basename(path)
  d <- fread(path, showProgress = FALSE)

  assert_cols(d, c("cell_id", celltype_col), file_label)
  pathways <- hallmark_cols(d, file_label)

  if (!is.null(keep_fun)) {
    d <- keep_fun(d)
  }

  if (nrow(d) == 0L) {
    stop(file_label, " has zero cells after sample filtering. Check sample rules.", call. = FALSE)
  }

  if (!is.null(sample_fun)) {
    sample_id <- sample_fun(d)
  } else {
    assert_cols(d, sample_col, file_label)
    sample_id <- as.character(d[[sample_col]])
  }

  raw_cell_type <- as.character(d[[celltype_col]])
  cell_type <- standardize_cell_type(raw_cell_type, celltype_map)

  out <- d[, c("cell_id", pathways), with = FALSE]
  out[, `:=`(
    model_id = model_id,
    source_id = source_id,
    sample_id = sample_id,
    raw_cell_type = raw_cell_type,
    cell_type = cell_type
  )]

  bad_sample <- is.na(out$sample_id) | out$sample_id == ""
  if (any(bad_sample)) {
    stop(
      file_label, " has ", sum(bad_sample),
      " cells with unresolved sample_id. Check sample parsing before continuing.",
      call. = FALSE
    )
  }

  out[, cell_id := paste(source_id, cell_id, sep = "::")]
  setcolorder(out, c(
    "model_id", "source_id", "sample_id",
    "raw_cell_type", "cell_type", "cell_id",
    pathways
  ))
  out[]
}

common_map <- c(
  "Fibroblast" = "FLS",
  "Fibroblast0" = "FLS",
  "Fibroblast1" = "FLS",
  "Fibroblast2" = "FLS",
  "Macrophage" = "Macrophage",
  "Monocyte" = "Monocyte",
  "T cell" = "T cell",
  "DC" = "DC",
  "Neutrophil" = "Neutrophil",
  "B cell" = "B cell",
  "B" = "B cell",
  "NK" = "NK",
  "Endothelial" = "Endothelial",
  "Pericyte" = "Pericyte",
  "Osteoblast" = "Osteoblast",
  "Chondrocyte" = "Chondrocyte",
  "Muscle" = "Muscle",
  "RBC" = "RBC",
  "Schwann" = "Schwann"
)

inputs <- list(
  standardize_scores(
    path = find_one_file("^GSE109449.*hallmark50.*\\.csv$"),
    model_id = "hsRA",
    source_id = "GSE109449_FLS",
    celltype_col = "celltype",
    celltype_map = common_map,
    keep_fun = function(d) {
      d[grepl("^(RA8|RA9)(_|$)", as.character(cell_id), perl = TRUE)]
    },
    sample_fun = function(d) {
      extract_regex(as.character(d$cell_id), "^(RA8|RA9)(?:_|$).*")
    }
  ),

  standardize_scores(
    path = find_one_file("^GAS003842.*hallmark50.*\\.csv$"),
    model_id = "hsRA",
    source_id = "GAS003842_infiltrating",
    celltype_col = "celltype",
    sample_col = "sample",
    celltype_map = common_map,
    keep_fun = function(d) d[as.character(sample) %in% paste0("RA", 1:6)]
  ),

  standardize_scores(
    path = find_one_file("^GSE129087.*hallmark50.*\\.csv$"),
    model_id = "STIA",
    source_id = "GSE129087_FLS",
    celltype_col = "celltype_manual",
    sample_col = "sample",
    celltype_map = common_map,
    keep_fun = function(d) d[as.character(sample) %in% paste0("inflamed_R", 1:3)]
  ),

  standardize_scores(
    path = find_one_file("^GSE254560.*hallmark50.*\\.csv$"),
    model_id = "STIA",
    source_id = "GSE254560_infiltrating",
    celltype_col = "celltype",
    sample_col = "sample",
    celltype_map = common_map,
    keep_fun = function(d) d[as.character(sample) == "PBS"]
  ),

  standardize_scores(
    path = find_one_file("^GSE184609.*hallmark50.*\\.csv$"),
    model_id = "GIA",
    source_id = "GSE184609",
    celltype_col = "celltype",
    sample_col = "sample",
    celltype_map = common_map,
    keep_fun = function(d) d[as.character(sample) %in% c("D6", "D14", "D25")]
  ),

  standardize_scores(
    path = find_one_file("^GSE192504.*hallmark50.*\\.csv$"),
    model_id = "CIA",
    source_id = "GSE192504",
    celltype_col = "celltype",
    sample_col = "sample",
    celltype_map = common_map,
    keep_fun = function(d) d[as.character(sample) == "CIA"]
  ),

  standardize_scores(
    path = find_one_file("C.*AIA.*Disease.*hallmark50.*\\.csv$"),
    model_id = "C&AIA",
    source_id = "C&AIA",
    celltype_col = "DefineTypes1",
    celltype_map = common_map,
    keep_fun = function(d) {
      if ("sample_id" %in% names(d)) {
        d[as.character(sample_id) %in% paste0("B", 1:3)]
      } else if ("sample" %in% names(d)) {
        d[as.character(sample) %in% paste0("B", 1:3)]
      } else {
        d[grepl("B[123]", as.character(cell_id), perl = TRUE)]
      }
    },
    sample_fun = function(d) {
      if ("sample_id" %in% names(d)) {
        as.character(d$sample_id)
      } else if ("sample" %in% names(d)) {
        as.character(d$sample)
      } else {
        extract_regex(as.character(d$cell_id), ".*(B[123]).*")
      }
    }
  )
)

cells <- rbindlist(inputs, use.names = TRUE)
pathways <- grep("^HALLMARK_", names(cells), value = TRUE)

sample_celltype_counts <- cells[, .(
  n_cells = .N
), by = .(model_id, source_id, sample_id, raw_cell_type, cell_type)]
setorder(sample_celltype_counts, model_id, source_id, cell_type, raw_cell_type, sample_id)

source_celltype_qc <- sample_celltype_counts[, .(
  n_samples = uniqueN(sample_id),
  n_cells = sum(n_cells),
  min_cells_per_sample = min(n_cells),
  max_cells_per_sample = max(n_cells)
), by = .(model_id, source_id, raw_cell_type, cell_type)]
setorder(source_celltype_qc, model_id, source_id, cell_type, raw_cell_type)

model_celltype_availability <- cells[, .(
  n_sources = uniqueN(source_id),
  n_samples = uniqueN(paste(source_id, sample_id, sep = "::")),
  n_cells = .N,
  raw_cell_types = paste(sort(unique(raw_cell_type)), collapse = "; ")
), by = .(model_id, cell_type)]
setorder(model_celltype_availability, model_id, cell_type)

sample_level_mean_AUCell <- cells[, lapply(.SD, mean, na.rm = TRUE),
  by = .(model_id, source_id, sample_id, cell_type),
  .SDcols = pathways
]
setorder(sample_level_mean_AUCell, model_id, source_id, cell_type, sample_id)

dataset_equal_sample_mean_AUCell <- sample_level_mean_AUCell[, lapply(.SD, mean, na.rm = TRUE),
  by = .(model_id, source_id, cell_type),
  .SDcols = pathways
]
setorder(dataset_equal_sample_mean_AUCell, model_id, source_id, cell_type)

model_equal_source_mean_AUCell <- dataset_equal_sample_mean_AUCell[, lapply(.SD, mean, na.rm = TRUE),
  by = .(model_id, cell_type),
  .SDcols = pathways
]
setorder(model_equal_source_mean_AUCell, model_id, cell_type)

metadata_assumptions <- data.table(
  model_id = c("hsRA", "hsRA", "STIA", "STIA", "GIA", "CIA", "C&AIA"),
  source_id = c(
    "GSE109449_FLS", "GAS003842_infiltrating", "GSE129087_FLS",
    "GSE254560_infiltrating", "GSE184609", "GSE192504", "C&AIA"
  ),
  samples_kept = c(
    "cell_id prefix RA8/RA9",
    "sample RA1-RA6",
    "sample inflamed_R1/inflamed_R2/inflamed_R3",
    "sample PBS",
    "sample D6/D14/D25",
    "sample CIA",
    "sample_id or sample B1/B2/B3; fallback parses B1/B2/B3 from cell_id"
  ),
  nominal_n = c(2L, 6L, 3L, 1L, 3L, 1L, 3L),
  cell_filtering_after_AUCell = "none; all cells in selected samples are retained",
  averaging_rule = "cell means within sample, then equal-weight mean across samples"
)

fwrite(cells, file.path(OUT_DIR, "harmonized_cell_level_AUCell.csv"))
saveRDS(cells, file.path(OUT_DIR, "harmonized_cell_level_AUCell.rds"))
fwrite(sample_celltype_counts, file.path(OUT_DIR, "sample_celltype_counts.csv"))
fwrite(source_celltype_qc, file.path(OUT_DIR, "source_celltype_qc.csv"))
fwrite(model_celltype_availability, file.path(OUT_DIR, "model_celltype_availability.csv"))
fwrite(sample_level_mean_AUCell, file.path(OUT_DIR, "sample_level_mean_AUCell.csv"))
fwrite(dataset_equal_sample_mean_AUCell, file.path(OUT_DIR, "dataset_equal_sample_mean_AUCell.csv"))
fwrite(model_equal_source_mean_AUCell, file.path(OUT_DIR, "model_equal_source_mean_AUCell.csv"))
fwrite(metadata_assumptions, file.path(OUT_DIR, "metadata_assumptions.csv"))

message("Step 1 complete.")
message("Output directory: ", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE))
message("Cells retained: ", nrow(cells))
message("Hallmark pathways: ", length(pathways))
message("No cell subgroup was removed by cell count or target-cell-type filtering.")
message("Check these files before Step 2:")
message("  - sample_celltype_counts.csv")
message("  - source_celltype_qc.csv")
message("  - model_celltype_availability.csv")
message("  - sample_level_mean_AUCell.csv")
message("  - dataset_equal_sample_mean_AUCell.csv")
message("  - model_equal_source_mean_AUCell.csv")
message("  - metadata_assumptions.csv")
