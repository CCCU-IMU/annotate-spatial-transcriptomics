#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratObject)
  library(data.table)
  library(jsonlite)
})

parse_args <- function(values) {
  result <- list()
  i <- 1L
  while (i <= length(values)) {
    if (!startsWith(values[[i]], "--") || i == length(values)) {
      stop("arguments must use --name value pairs")
    }
    result[[substring(values[[i]], 3L)]] <- values[[i + 1L]]
    i <- i + 2L
  }
  result
}
sha256 <- function(path) {
  output <- system2("sha256sum", shQuote(path), stdout = TRUE)
  strsplit(output[[1L]], "\\s+")[[1L]][[1L]]
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("rds", "mapping", "out", "input-sha256")
missing <- required[!required %in% names(args)]
if (length(missing)) stop("missing: ", paste(missing, collapse = ", "))
input <- normalizePath(args$rds, mustWork = TRUE)
mapping_path <- normalizePath(args$mapping, mustWork = TRUE)
out <- args$out
dir.create(out, recursive = TRUE, showWarnings = FALSE)

mapping <- fread(mapping_path, colClasses = "character")
if (!all(c("resolution", "cluster_column") %in% names(mapping))) {
  stop("mapping must contain resolution and cluster_column")
}
mapping[, resolution_numeric := as.numeric(resolution)]
if (nrow(mapping) < 3L || any(!is.finite(mapping$resolution_numeric)) ||
    anyDuplicated(mapping$resolution_numeric)) {
  stop("partition mapping requires at least three unique numeric resolutions")
}
setorder(mapping, resolution_numeric)

obj <- readRDS(input)
if (!inherits(obj, "Seurat")) stop("input is not a Seurat object")
metadata <- as.data.table(obj[[]], keep.rownames = "cell_id")
if (!all(mapping$cluster_column %in% names(metadata))) {
  stop("Seurat metadata lacks mapped cluster columns")
}
rows <- vector("list", nrow(mapping))
for (i in seq_len(nrow(mapping))) {
  cluster <- as.character(metadata[[mapping$cluster_column[[i]]]])
  if (anyNA(cluster) || any(!nzchar(cluster))) {
    stop("cluster column has missing values: ", mapping$cluster_column[[i]])
  }
  rows[[i]] <- data.table(
    cell_id = metadata$cell_id,
    boundary_id = "whole_tissue",
    resolution = mapping$resolution_numeric[[i]],
    cluster = cluster,
    resolution_role = ifelse(i == 1L, "selected", "grid")
  )
}
partition_path <- file.path(out, "partition_grid.tsv.gz")
fwrite(rbindlist(rows), partition_path, sep = "\t")

layer_audit <- list()
for (assay in Assays(obj)) {
  layer_audit[[assay]] <- Layers(obj[[assay]])
}
spatial_data_counts_identical_sample <- NA
if ("Spatial" %in% Assays(obj) &&
    all(c("counts", "data") %in% Layers(obj[["Spatial"]]))) {
  counts <- LayerData(obj[["Spatial"]], layer = "counts")
  data <- LayerData(obj[["Spatial"]], layer = "data")
  feature_index <- seq_len(min(1000L, nrow(counts), nrow(data)))
  cell_index <- seq_len(min(1000L, ncol(counts), ncol(data)))
  spatial_data_counts_identical_sample <- isTRUE(all.equal(
    as.matrix(counts[feature_index, cell_index, drop = FALSE]),
    as.matrix(data[feature_index, cell_index, drop = FALSE]),
    check.attributes = FALSE
  ))
}
manifest <- list(
  status = "PASS",
  schema_version = "2.2",
  input_rds = input,
  input_sha256 = args[["input-sha256"]],
  mapping = list(path = mapping_path, sha256 = sha256(mapping_path)),
  partition_grid = list(
    path = normalizePath(partition_path, mustWork = TRUE),
    sha256 = sha256(partition_path),
    resolutions = mapping$resolution_numeric,
    n_observations = nrow(metadata)
  ),
  expression_layer_audit = list(
    assays = layer_audit,
    spatial_data_counts_identical_sample =
      spatial_data_counts_identical_sample
  ),
  output_columns = c(
    "cell_id", "boundary_id", "resolution", "cluster", "resolution_role"
  ),
  historical_annotation_exported = FALSE
)
write_json(
  manifest, file.path(out, "partition_grid_manifest.json"),
  auto_unbox = TRUE, pretty = TRUE, null = "null"
)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
writeLines("PASS", file.path(out, "RUN_COMPLETE.tsv"))
