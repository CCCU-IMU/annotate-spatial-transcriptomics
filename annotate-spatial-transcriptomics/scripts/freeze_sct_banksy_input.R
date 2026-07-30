#!/usr/bin/env Rscript

## Freeze a project-local SCT+BANKSY Seurat input without reading annotations.
suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("Usage: freeze_sct_banksy_input.R <rds> <sha256sum> <cluster_mapping.tsv|AUTO> <preprocess_manifest.json|AUTO> <sample_id> <out>")
}
input_rds <- normalizePath(args[[1L]], mustWork = TRUE)
sha256_record <- normalizePath(args[[2L]], mustWork = TRUE)
mapping_arg <- args[[3L]]
preprocess_manifest_arg <- args[[4L]]
sample_id <- args[[5L]]
out <- args[[6L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) {
  value <- system2("sha256sum", shQuote(path), stdout = TRUE)
  strsplit(value[[1L]], "\\s+")[[1L]][[1L]]
}
write_tsv <- function(value, path) {
  if (grepl("\\.gz$", path) && nrow(value) == 0L) {
    con <- gzfile(path, open = "wt")
    on.exit(close(con), add = TRUE)
    writeLines(paste(names(value), collapse = "\t"), con)
  } else {
    fwrite(value, path, sep = "\t")
  }
}

sha_parts <- strsplit(trimws(readLines(sha256_record, warn = FALSE)[[1L]]), "\\s+")[[1L]]
if (length(sha_parts) < 2L || nchar(sha_parts[[1L]]) != 64L ||
    normalizePath(sub("^\\*", "", sha_parts[[length(sha_parts)]]), mustWork = TRUE) != input_rds) {
  stop("invalid or mismatched input sha256sum record")
}

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("input is not a Seurat object")
metadata <- as.data.table(obj[[]], keep.rownames = "cell_id")

infer_banksy_mapping <- function(columns) {
  candidates <- columns[
    grepl("banksy", columns, ignore.case = TRUE) &
      grepl("res|resolution|cluster|snn", columns, ignore.case = TRUE)
  ]
  if (length(candidates) < 3L) {
    stop(paste0(
      "AUTO cluster mapping found fewer than three BANKSY resolution columns; ",
      "provide an explicit mapping TSV with resolution and cluster_column"
    ))
  }
  values <- vapply(candidates, function(column) {
    hits <- regmatches(column, gregexpr("[0-9]+(?:\\.[0-9]+)?", column, perl = TRUE))[[1L]]
    if (!length(hits)) return(NA_real_)
    as.numeric(tail(hits, 1L))
  }, numeric(1L))
  keep <- is.finite(values) & values > 0 & values <= 5
  candidates <- candidates[keep]; values <- values[keep]
  if (length(candidates) < 3L || anyDuplicated(values)) {
    stop(paste0(
      "AUTO BANKSY resolution inference is absent or ambiguous; ",
      "provide an explicit mapping TSV"
    ))
  }
  data.table(resolution = as.character(values), cluster_column = candidates)
}

if (identical(toupper(mapping_arg), "AUTO")) {
  mapping <- infer_banksy_mapping(names(metadata))
  mapping_path <- file.path(out, "inferred_banksy_cluster_mapping.tsv")
  fwrite(mapping, mapping_path, sep = "\t")
} else {
  mapping_path <- normalizePath(mapping_arg, mustWork = TRUE)
  mapping <- fread(mapping_path, colClasses = "character")
}
if (!all(c("resolution", "cluster_column") %in% names(mapping))) {
  stop("cluster mapping must contain resolution and cluster_column")
}
mapping[, resolution_numeric := as.numeric(resolution)]
if (nrow(mapping) < 3L || any(!is.finite(mapping$resolution_numeric)) ||
    anyDuplicated(mapping$resolution_numeric)) stop("invalid whole-tissue resolution grid")
setorder(mapping, resolution_numeric)

if (!all(mapping$cluster_column %in% names(metadata))) {
  stop("input metadata lacks one or more bound BANKSY cluster columns")
}
coordinate_pairs <- list(c("x", "y"), c("sdimx", "sdimy"),
                         c("spatial_x", "spatial_y"))
coordinate_pair <- NULL
for (pair in coordinate_pairs) {
  if (all(pair %in% names(metadata))) {
    coordinate_pair <- pair
    break
  }
}
if (is.null(coordinate_pair)) stop("input metadata lacks a supported x/y coordinate pair")
x_values <- as.numeric(metadata[[coordinate_pair[[1L]]]])
y_values <- as.numeric(metadata[[coordinate_pair[[2L]]]])
assays <- Assays(obj)
count_backed <- assays[vapply(assays, function(assay) {
  assay != "SCT" && "counts" %in% Layers(obj[[assay]])
}, logical(1))]
preferred <- c("RNA", "Spatial")
raw_assay <- preferred[preferred %in% count_backed][1L]
if (!length(raw_assay) || is.na(raw_assay)) raw_assay <- count_backed[1L]
if (!length(raw_assay) || is.na(raw_assay) || identical(raw_assay, "SCT")) {
  stop("no project-local non-SCT raw-count assay is available")
}
counts <- LayerData(obj[[raw_assay]], layer = "counts")
if (!inherits(counts, "sparseMatrix")) stop("raw counts are not sparse")
if (any(!is.finite(counts@x)) || any(counts@x < 0) ||
    any(abs(counts@x - round(counts@x)) > 1e-10)) {
  stop("raw counts are not finite nonnegative integers")
}
if (!identical(colnames(counts), metadata$cell_id)) {
  stop("raw-count columns and Seurat metadata are not identically ordered")
}

ncount <- Matrix::colSums(counts)
nfeature <- Matrix::colSums(counts > 0)
finite_xy <- is.finite(x_values) & is.finite(y_values)
analysis_flag <- ncount >= 100 & nfeature >= 75 & finite_xy
analysis <- data.table(
  cell_id = metadata$cell_id[analysis_flag], x = x_values[analysis_flag],
  y = y_values[analysis_flag], analysis_scope = "analysis_set"
)
excluded <- data.table(
  cell_id = metadata$cell_id[!analysis_flag], x = x_values[!analysis_flag],
  y = y_values[!analysis_flag], analysis_scope = "excluded_initial_qc",
  exclusion_reason = ifelse(
    !finite_xy[!analysis_flag], "invalid_spatial_coordinate",
    ifelse(ncount[!analysis_flag] < 100, "nCount_below_100", "nFeature_below_75")
  )
)
if (!nrow(analysis)) stop("analysis set is empty")
write_tsv(analysis, file.path(out, "analysis_membership.tsv.gz"))
write_tsv(excluded, file.path(out, "excluded_initial_qc.tsv.gz"))
write_tsv(rbindlist(list(analysis, excluded[, .(cell_id, x, y, analysis_scope)])),
          file.path(out, "analysis_scope.tsv.gz"))

partitions <- rbindlist(lapply(seq_len(nrow(mapping)), function(index) {
  cluster <- as.character(metadata[[mapping$cluster_column[[index]]]])
  if (anyNA(cluster) || any(!nzchar(cluster))) stop("BANKSY cluster column contains missing labels")
  data.table(
    cell_id = metadata$cell_id[analysis_flag], boundary_id = "whole_tissue",
    resolution = mapping$resolution_numeric[[index]], cluster = cluster[analysis_flag],
    resolution_role = "grid"
  )
}))
write_tsv(partitions, file.path(out, "partition_grid.tsv.gz"))

if (identical(toupper(preprocess_manifest_arg), "AUTO")) {
  preprocess_manifest_path <- file.path(out, "inferred_upstream_preprocess_manifest.json")
  write_json(list(
    schema_version = "1.0", sample_id = sample_id,
    source_rds = input_rds, source_rds_sha256 = tolower(sha_parts[[1L]]),
    declared_preprocessing = "SCT+BANKSY already present in supplied object",
    inference_scope = "assay/layer and BANKSY metadata audit only",
    biological_labels_read = FALSE
  ), preprocess_manifest_path, pretty = TRUE, auto_unbox = TRUE)
} else {
  preprocess_manifest_path <- normalizePath(preprocess_manifest_arg, mustWork = TRUE)
}
preprocess_manifest <- fromJSON(preprocess_manifest_path, simplifyVector = FALSE)
if (!identical(as.character(preprocess_manifest$sample_id), sample_id)) {
  stop("preprocessing manifest sample differs from runtime sample")
}
per_resolution <- lapply(seq_len(nrow(mapping)), function(index) list(
  resolution = mapping$resolution_numeric[[index]],
  cluster_column = mapping$cluster_column[[index]],
  n_clusters = uniqueN(metadata[[mapping$cluster_column[[index]]]][analysis_flag])
))
write_json(list(
  schema_version = "2.2", candidate_resolutions = mapping$resolution_numeric,
  method = "BANKSY", source = "bound_upstream_input",
  cluster_columns = mapping$cluster_column,
  preprocessing_manifest = list(path = preprocess_manifest_path, sha256 = sha256(preprocess_manifest_path))
), file.path(out, "whole_tissue_grid.json"), pretty = TRUE, auto_unbox = TRUE)

label_pattern <- paste(c(
  "celltype", "cell_type", "annotation", "broad", "fine", "subtype",
  "historical", "repair", "atlas", "predicted", "label"
), collapse = "|")
write_json(list(
  status = "PASS", schema_version = "2.2", sample_id = sample_id,
  blind_runtime = TRUE, input_rds = input_rds,
  input_sha256 = tolower(sha_parts[[1L]]), n_observations = nrow(metadata),
  analysis_set_n = nrow(analysis), excluded_initial_qc_n = nrow(excluded),
  n_features_raw = nrow(counts), raw_count_assay = raw_assay,
  default_assay = DefaultAssay(obj), assays = assays,
  layers = setNames(lapply(assays, function(assay) Layers(obj[[assay]])), assays),
  coordinate_columns = coordinate_pair, whole_tissue_partitions = per_resolution,
  metadata_columns = names(metadata),
  historical_like_metadata_columns = grep(label_pattern, names(metadata), ignore.case = TRUE, value = TRUE),
  scoring_exports_historical_columns = FALSE,
  expression_boundary = "project-local non-SCT raw counts"
), file.path(out, "input_audit_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines("PASS", file.path(out, "RUN_COMPLETE.tsv"))
