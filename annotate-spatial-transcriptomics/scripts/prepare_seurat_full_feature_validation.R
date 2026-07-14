#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(digest)
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

a <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("rds", "membership", "out", "sample")
missing <- required[!required %in% names(a)]
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))

assay <- ifelse(is.null(a$assay), "Spatial", a$assay)
count_layer <- ifelse(is.null(a$`count-layer`), "counts", a$`count-layer`)
data_layer <- ifelse(is.null(a$`data-layer`), "data", a$`data-layer`)
cell_id_col <- ifelse(is.null(a$`cell-id-col`), "cell_id", a$`cell-id-col`)
include_col <- ifelse(is.null(a$`include-col`), "analysis_set", a$`include-col`)
include_values <- strsplit(ifelse(is.null(a$`include-values`), "TRUE", a$`include-values`), ",", fixed = TRUE)[[1]]
scale_factor <- as.numeric(ifelse(is.null(a$`scale-factor`), 10000, a$`scale-factor`))
if (!is.finite(scale_factor) || scale_factor != 10000) {
  stop("The validation contract requires --scale-factor 10000")
}

out <- normalizePath(a$out, mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "tables"), showWarnings = FALSE)
dir.create(file.path(out, "provenance"), showWarnings = FALSE)

read_any <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) fread(cmd = paste("gzip -dc", shQuote(path))) else fread(path)
}
hash_file <- function(path) digest(path, algo = "sha256", file = TRUE)
hash_ids <- function(ids) digest(paste(sort(as.character(ids)), collapse = "\n"), algo = "sha256")
layer_matrix <- function(object, assay_name, layer_name) {
  tryCatch(LayerData(object[[assay_name]], layer = layer_name),
           error = function(e) GetAssayData(object[[assay_name]], slot = layer_name))
}
sparse_exact_equal <- function(left, right) {
  left <- methods::as(left, "dgCMatrix"); right <- methods::as(right, "dgCMatrix")
  identical(dim(left), dim(right)) && identical(left@p, right@p) &&
    identical(left@i, right@i) && identical(as.numeric(left@x), as.numeric(right@x))
}

message("Reading immutable raw-count Seurat container")
source_object <- readRDS(a$rds)
if (!inherits(source_object, "Seurat")) stop("--rds is not a Seurat object")
if (!assay %in% Assays(source_object)) stop("Missing assay: ", assay)
membership <- read_any(a$membership)
if (!cell_id_col %in% names(membership) || !include_col %in% names(membership)) {
  stop("Membership lacks required columns: ", cell_id_col, ", ", include_col)
}
membership[[cell_id_col]] <- as.character(membership[[cell_id_col]])
if (anyDuplicated(membership[[cell_id_col]])) stop("Membership contains duplicate observation IDs")
if (nrow(membership) != ncol(source_object) ||
    !setequal(membership[[cell_id_col]], colnames(source_object))) {
  stop("Membership must account for every observation in the immutable source object")
}
selected <- tolower(as.character(membership[[include_col]])) %in% tolower(include_values)
analysis_ids <- membership[selected, get(cell_id_col)]
analysis_ids <- colnames(source_object)[colnames(source_object) %in% analysis_ids]
if (length(analysis_ids) < 3L) stop("Fewer than three analysis-set observations")

raw_counts <- layer_matrix(source_object, assay, count_layer)[, analysis_ids, drop = FALSE]
metadata <- source_object[[]][analysis_ids, , drop = FALSE]
validation <- CreateSeuratObject(counts = raw_counts, assay = assay, meta.data = metadata)
if (!identical(rownames(validation), rownames(raw_counts)) ||
    !identical(colnames(validation), colnames(raw_counts))) {
  stop("Creating the validation object changed feature or observation identities")
}
validation@reductions <- list()
DefaultAssay(validation) <- assay
message("Building full-feature LogNormalize validation layer")
validation <- NormalizeData(validation, assay = assay, normalization.method = "LogNormalize",
                            scale.factor = scale_factor, verbose = TRUE)
normalized <- layer_matrix(validation, assay, data_layer)
if (sparse_exact_equal(layer_matrix(validation, assay, count_layer), normalized)) {
  stop("LogNormalize output is exactly identical to counts")
}
normalized_sparse <- methods::as(normalized, "dgCMatrix")
if (!length(normalized_sparse@x) || !any(abs(normalized_sparse@x - round(normalized_sparse@x)) > 1e-8)) {
  stop("LogNormalize output lacks non-integer normalized values")
}

selected_membership <- data.table(cell_id = colnames(validation), analysis_set = TRUE)
membership_path <- file.path(out, "tables", "validation_analysis_set.tsv.gz")
fwrite(selected_membership, membership_path, sep = "\t")
analysis_sha <- hash_ids(selected_membership$cell_id)
source_membership_sha <- hash_file(a$membership)
input_sha <- hash_file(a$rds)

validation@misc$full_feature_validation_contract <- list(
  role = "full_feature_deg_marker_validation_only",
  clustering_eligible = FALSE,
  source_input_sha256 = input_sha,
  source_membership_sha256 = source_membership_sha,
  analysis_set_sha256 = analysis_sha,
  normalization_method = "LogNormalize",
  scale_factor = scale_factor,
  assay = assay,
  count_layer = count_layer,
  data_layer = data_layer
)
object_path <- file.path(out, paste0(a$sample, "_seurat_full_feature_lognormalized_validation.rds"))
saveRDS(validation, object_path, compress = FALSE)
object_sha <- hash_file(object_path)

manifest <- list(
  status = "PASS",
  sample = a$sample,
  role = "full_feature_deg_marker_validation_only",
  clustering_eligible = FALSE,
  source_rds = normalizePath(a$rds),
  source_rds_sha256 = input_sha,
  source_membership = normalizePath(a$membership),
  source_membership_sha256 = source_membership_sha,
  analysis_set_sha256 = analysis_sha,
  assay = assay,
  count_layer = count_layer,
  data_layer = data_layer,
  normalization_method = "LogNormalize",
  scale_factor = scale_factor,
  n_features = nrow(validation),
  n_observations = ncol(validation),
  counts_equal_data = FALSE,
  reductions_removed = TRUE,
  normalized_object = normalizePath(object_path),
  normalized_object_sha256 = object_sha
)
write_json(manifest, file.path(out, "provenance", "full_feature_validation_manifest.json"),
           pretty = TRUE, auto_unbox = TRUE)
capture.output(sessionInfo(), file = file.path(out, "provenance", "sessionInfo.txt"))
writeLines(c("status\tPASS", paste0("analysis_set_sha256\t", analysis_sha),
             paste0("normalized_object_sha256\t", object_sha)), file.path(out, "RUN_COMPLETE.tsv"))
message("Full-feature LogNormalize validation object complete")
