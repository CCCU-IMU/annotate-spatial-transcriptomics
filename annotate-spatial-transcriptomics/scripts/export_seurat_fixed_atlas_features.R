#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

parse_args <- function(values) {
  result <- list(); index <- 1L
  while (index <= length(values)) {
    if (!startsWith(values[[index]], "--") || index == length(values)) {
      stop("arguments must use --name value pairs")
    }
    result[[substring(values[[index]], 3L)]] <- values[[index + 1L]]
    index <- index + 2L
  }
  result
}
read_tsv <- function(path) {
  if (grepl("\\.gz$", path)) fread(cmd = paste("gzip -dc", shQuote(path))) else fread(path)
}
sha256 <- function(path) {
  value <- system2("sha256sum", shQuote(path), stdout = TRUE)
  strsplit(value[[1L]], "\\s+")[[1L]][[1L]]
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("rds", "membership", "features", "out")
missing <- required[!required %in% names(args)]
if (length(missing)) stop("missing: ", paste(missing, collapse = ", "))

input <- normalizePath(args$rds, mustWork = TRUE)
membership_path <- normalizePath(args$membership, mustWork = TRUE)
feature_path <- normalizePath(args$features, mustWork = TRUE)
out <- args$out
dir.create(out, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input)
if (!inherits(obj, "Seurat")) stop("--rds is not a Seurat object")
membership <- read_tsv(membership_path)
if (!"cell_id" %in% names(membership) || anyDuplicated(membership$cell_id) ||
    any(!nzchar(as.character(membership$cell_id)))) {
  stop("membership must contain unique nonempty cell_id")
}
ids <- as.character(membership$cell_id)
if (!all(ids %in% colnames(obj))) stop("membership contains cells absent from the Seurat object")

features <- as.character(fread(feature_path)$gene)
if (!length(features) || anyDuplicated(features) || any(!nzchar(features))) {
  stop("fixed Atlas features must be unique and nonempty")
}
assays <- Assays(obj)
count_backed <- assays[vapply(assays, function(assay) {
  assay != "SCT" && "counts" %in% Layers(obj[[assay]])
}, logical(1))]
preferred <- c("RNA", "Spatial")
assay <- preferred[preferred %in% count_backed][1L]
if (!length(assay) || is.na(assay)) assay <- count_backed[1L]
if (!length(assay) || is.na(assay) || identical(assay, "SCT")) {
  stop("fixed Atlas export requires a project-local non-SCT raw-count assay")
}
counts <- LayerData(obj[[assay]], layer = "counts")
if (any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-10)) {
  stop("raw-count layer is not nonnegative integer-like")
}
present <- features %in% rownames(counts)
if (sum(present) < 500L) stop("fewer than 500 fixed Atlas features are present")
matrix_value <- counts[features[present], match(ids, colnames(counts)), drop = FALSE]
if (any(!present)) {
  zero <- Matrix(0, nrow = sum(!present), ncol = length(ids), sparse = TRUE)
  rownames(zero) <- features[!present]
  colnames(zero) <- ids
  matrix_value <- rbind(matrix_value, zero)[features, , drop = FALSE]
}
if (!identical(rownames(matrix_value), features) || !identical(colnames(matrix_value), ids)) {
  stop("fixed-feature export order is not exact")
}

matrix_path <- file.path(out, "query_cells_by_features.mtx")
writeMM(t(matrix_value), matrix_path)
fwrite(data.table(cell_id = ids), file.path(out, "cells.tsv"), sep = "\t")
fwrite(data.table(gene = features), file.path(out, "features.tsv"), sep = "\t")
write_json(list(
  status = "PASS",
  schema_version = "2.2",
  input_rds = input,
  membership = list(path = membership_path, sha256 = sha256(membership_path)),
  raw_count_assay = assay,
  n_query = length(ids),
  n_features = length(features),
  n_present = sum(present),
  n_missing_zero_filled = sum(!present),
  query_reference_joint_retraining = FALSE,
  historical_labels_read = FALSE
), file.path(out, "export_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
