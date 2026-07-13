#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

parse_args <- function(x) {
  out <- list(); i <- 1L
  while (i <= length(x)) {
    key <- sub("^--", "", x[[i]])
    if (i == length(x) || startsWith(x[[i + 1L]], "--")) {
      out[[key]] <- TRUE; i <- i + 1L
    } else {
      out[[key]] <- x[[i + 1L]]; i <- i + 2L
    }
  }
  out
}

read_any <- function(path) {
  if (grepl("\\.gz$", path)) fread(cmd = paste("gzip -dc", shQuote(path))) else fread(path)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("rds", "config", "out", "cell-id-col")
missing <- required[!required %in% names(args)]
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))

cfg <- fromJSON(args$config, simplifyVector = FALSE)
if (!length(cfg$programs)) stop("config must contain programs")
obj <- readRDS(args$rds)
assay <- if (is.null(args$assay)) DefaultAssay(obj) else args$assay
layer <- if (is.null(args$layer)) "counts" else args$layer
mat <- tryCatch(
  LayerData(obj[[assay]], layer = layer),
  error = function(e) GetAssayData(obj[[assay]], slot = layer)
)
cell_col <- args$`cell-id-col`
quantiles <- c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1)
default_total_thresholds <- c(0, 100, 200, 300)

metric_outputs <- list()
marker_outputs <- list()
group_outputs <- list()
grid_outputs <- list()

for (i in seq_along(cfg$programs)) {
  spec <- cfg$programs[[i]]
  membership <- read_any(spec$membership)
  if (!cell_col %in% names(membership)) stop("cell ID column missing from ", spec$membership)
  ids <- intersect(as.character(membership[[cell_col]]), colnames(obj))
  requested_pos <- unique(as.character(unlist(spec$positive_markers)))
  requested_anti <- unique(as.character(unlist(spec$anti_markers)))
  pos <- intersect(requested_pos, rownames(mat))
  anti <- intersect(requested_anti, rownames(mat))
  if (!length(ids) || !length(pos)) stop("empty cells or positive marker set for ", spec$label)
  x <- mat[, ids, drop = FALSE]

  metrics <- list(
    n_positive_detected = Matrix::colSums(x[pos, , drop = FALSE] > 0),
    n_anti_detected = if (length(anti)) Matrix::colSums(x[anti, , drop = FALSE] > 0) else rep(0, length(ids)),
    positive_score = Matrix::colSums(x[pos, , drop = FALSE]),
    total_counts = Matrix::colSums(x)
  )
  metric_outputs[[i]] <- rbindlist(lapply(names(metrics), function(metric) {
    q <- quantile(metrics[[metric]], probs = quantiles, na.rm = TRUE, names = FALSE)
    data.table(
      label = spec$label,
      membership = spec$membership,
      n_cells = length(ids),
      metric = metric,
      quantile = quantiles,
      value = as.numeric(q),
      positive_markers_present = paste(pos, collapse = ","),
      anti_markers_present = paste(anti, collapse = ",")
    )
  }))

  requested <- rbindlist(list(
    data.table(marker = requested_pos, marker_role = "positive"),
    data.table(marker = requested_anti, marker_role = "anti")
  ))
  requested[, `:=`(
    label = spec$label,
    membership = spec$membership,
    n_cells = length(ids),
    feature_present = marker %in% rownames(mat),
    detection_fraction = NA_real_,
    mean_counts = NA_real_
  )]
  present_markers <- requested[feature_present == TRUE, marker]
  if (length(present_markers)) {
    requested[feature_present == TRUE, detection_fraction := as.numeric(Matrix::rowMeans(x[present_markers, , drop = FALSE] > 0))]
    requested[feature_present == TRUE, mean_counts := as.numeric(Matrix::rowMeans(x[present_markers, , drop = FALSE]))]
  }
  marker_outputs[[i]] <- requested

  groups <- spec$required_positive_groups
  group_pass <- rep(TRUE, length(ids))
  if (!is.null(groups) && length(groups)) {
    group_rows <- list()
    for (j in seq_along(groups)) {
      g <- groups[[j]]
      group_name <- if (is.null(g$name)) paste0("group_", j) else as.character(g$name)
      group_requested <- unique(as.character(unlist(g$markers)))
      group_present <- intersect(group_requested, rownames(mat))
      min_detected <- if (is.null(g$min_detected)) 1L else as.integer(g$min_detected)
      values <- if (length(group_present)) Matrix::colSums(x[group_present, , drop = FALSE] > 0) else rep(0, length(ids))
      group_pass <- group_pass & values >= min_detected
      q <- quantile(values, probs = quantiles, na.rm = TRUE, names = FALSE)
      group_rows[[j]] <- data.table(
        label = spec$label,
        membership = spec$membership,
        group_name = group_name,
        n_cells = length(ids),
        min_detected = min_detected,
        markers_present = paste(group_present, collapse = ","),
        markers_missing = paste(setdiff(group_requested, group_present), collapse = ","),
        quantile = quantiles,
        value = as.numeric(q),
        n_group_pass = sum(values >= min_detected),
        group_pass_fraction = mean(values >= min_detected)
      )
    }
    group_outputs[[i]] <- rbindlist(group_rows)
  }

  positive_thresholds <- if (is.null(spec$profile_positive_thresholds)) {
    seq_len(min(length(pos), 6L))
  } else as.integer(unlist(spec$profile_positive_thresholds))
  anti_thresholds <- if (is.null(spec$profile_anti_thresholds)) {
    seq.int(0L, min(length(anti), 3L))
  } else as.integer(unlist(spec$profile_anti_thresholds))
  total_thresholds <- if (is.null(spec$profile_total_count_thresholds)) {
    default_total_thresholds
  } else as.numeric(unlist(spec$profile_total_count_thresholds))
  grid <- CJ(
    min_positive_detected = positive_thresholds,
    max_anti_detected = anti_thresholds,
    min_total_counts = total_thresholds,
    unique = TRUE
  )
  grid[, c("n_pass", "pass_fraction") := {
    keep <- group_pass &
      metrics$n_positive_detected >= min_positive_detected &
      metrics$n_anti_detected <= max_anti_detected &
      metrics$total_counts >= min_total_counts
    list(sum(keep), mean(keep))
  }, by = .(min_positive_detected, max_anti_detected, min_total_counts)]
  grid[, `:=`(
    label = spec$label,
    membership = spec$membership,
    n_cells = length(ids),
    required_groups_applied = !is.null(groups) && length(groups) > 0
  )]
  grid_outputs[[i]] <- grid
}

out <- rbindlist(metric_outputs, fill = TRUE)
marker_out <- rbindlist(marker_outputs, fill = TRUE)
group_out <- if (length(group_outputs)) rbindlist(group_outputs, fill = TRUE) else data.table()
grid_out <- rbindlist(grid_outputs, fill = TRUE)
dir.create(dirname(args$out), recursive = TRUE, showWarnings = FALSE)
fwrite(out, args$out, sep = "\t")
fwrite(marker_out, paste0(args$out, ".marker_detection.tsv"), sep = "\t")
fwrite(group_out, paste0(args$out, ".required_group_detection.tsv"), sep = "\t")
fwrite(grid_out, paste0(args$out, ".threshold_grid.tsv"), sep = "\t")
cat(toJSON(list(
  status = "PASS",
  output = normalizePath(args$out),
  marker_detection = normalizePath(paste0(args$out, ".marker_detection.tsv")),
  required_group_detection = normalizePath(paste0(args$out, ".required_group_detection.tsv")),
  threshold_grid = normalizePath(paste0(args$out, ".threshold_grid.tsv")),
  n_programs = length(metric_outputs)
), auto_unbox = TRUE, pretty = TRUE), "\n")
