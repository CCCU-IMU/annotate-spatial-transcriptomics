#!/usr/bin/env Rscript

## Materialize the v2 broad-family evidence matrix from the current project
## expression object.  This script reports raw absolute measurements only; it
## does not choose labels or infer absence from a centered score.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(digest)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot resolve script directory")
source(file.path(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))),
                 "seurat_validation_layer.R"))

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
  if (grepl("\\.gz$", path, ignore.case = TRUE)) fread(cmd = paste("gzip -dc", shQuote(path))) else fread(path)
}

sha256_file <- function(path) digest(file = path, algo = "sha256", serialize = FALSE)

resolve_dotted <- function(value, dotted) {
  for (key in strsplit(dotted, "\\.", fixed = FALSE)[[1]]) {
    if (is.null(value[[key]])) stop("Cannot resolve profile program: ", dotted)
    value <- value[[key]]
  }
  value
}

a <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("rds", "clusters", "profile", "catalog", "validation-manifest", "out", "cell-id-col", "cluster-col")
missing <- required[!required %in% names(a)]
if (length(missing)) stop("Missing: ", paste(missing, collapse = ", "))

obj <- readRDS(a$rds)
if (!inherits(obj, "Seurat")) stop("--rds must contain a Seurat object")
validation_record <- fromJSON(a$`validation-manifest`, simplifyVector = TRUE)
if (!identical(as.character(validation_record$status), "PASS")) {
  stop("--validation-manifest must be a PASS full-feature validation record")
}
assay <- if (is.null(a$assay)) DefaultAssay(obj) else a$assay
count_layer <- if (is.null(a$`count-layer`)) "counts" else a$`count-layer`
data_layer <- if (is.null(a$`data-layer`)) "data" else a$`data-layer`
assert_seurat_validation_layer(
  obj, assay = assay, data_layer = data_layer, count_layer = count_layer,
  manifest_path = a$`validation-manifest`, object_path = a$rds
)

counts <- tryCatch(LayerData(obj[[assay]], layer = count_layer),
                   error = function(e) GetAssayData(obj[[assay]], slot = count_layer))
normalized <- tryCatch(LayerData(obj[[assay]], layer = data_layer),
                       error = function(e) GetAssayData(obj[[assay]], slot = data_layer))
clusters <- read_any(a$clusters)
cell_col <- a$`cell-id-col`; cluster_col <- a$`cluster-col`
if (!all(c(cell_col, cluster_col) %in% names(clusters))) stop("cluster table lacks requested columns")
clusters[[cell_col]] <- as.character(clusters[[cell_col]])
if (uniqueN(clusters[[cell_col]]) != nrow(clusters)) stop("cluster cell IDs are not unique")
if (setequal(clusters[[cell_col]], colnames(obj)) == FALSE) stop("cluster membership must match the expression object exactly")
if (!is.null(validation_record$n_observations) && as.integer(validation_record$n_observations) != ncol(obj)) {
  stop("validation manifest observation count differs from the expression object")
}
clusters <- clusters[match(colnames(obj), get(cell_col))]

profile <- fromJSON(a$profile, simplifyVector = FALSE)
catalog <- fromJSON(a$catalog, simplifyVector = FALSE)
specs <- catalog$candidate_boundaries
if (is.null(specs) || length(specs) < 2L) stop("candidate catalog is empty")

rows <- list(); cursor <- 1L
cluster_ids <- sort(unique(as.character(clusters[[cluster_col]])))
candidate_ids <- character()
for (candidate in specs) {
  if (!isTRUE(candidate$review_required)) next
  candidate_id <- as.character(candidate$candidate_id)
  candidate_ids <- c(candidate_ids, candidate_id)
  program <- resolve_dotted(profile, as.character(candidate$profile_program))
  families <- program$positive_families
  if (is.null(families) || length(families) < 2L) {
    stop(candidate_id, " lacks two explicit positive_families in the bound profile")
  }
  for (cluster_id in cluster_ids) {
    ids <- clusters[get(cluster_col) == cluster_id, get(cell_col)]
    for (family_name in names(families)) {
      requested <- unique(as.character(unlist(families[[family_name]])))
      available <- intersect(requested, rownames(counts))
      missing_markers <- setdiff(requested, available)
      if (length(available)) {
        raw <- counts[available, ids, drop = FALSE]
        dat <- normalized[available, ids, drop = FALSE]
        hits_per_obs <- Matrix::colSums(raw > 0)
        marker_totals <- Matrix::rowSums(raw)
        pseudobulk <- sum(marker_totals)
        mean_norm <- as.numeric(Matrix::sum(dat)) / (length(available) * length(ids))
        any_fraction <- mean(hits_per_obs > 0)
        median_hits <- median(hits_per_obs)
        detected_marker_count <- sum(marker_totals > 0)
      } else {
        pseudobulk <- 0; mean_norm <- 0; any_fraction <- 0; median_hits <- 0; detected_marker_count <- 0L
      }
      rows[[cursor]] <- data.table(
        cluster = cluster_id,
        n_observations = length(ids),
        candidate_id = candidate_id,
        candidate_family = as.character(candidate$family),
        release_level = as.character(candidate$release_level),
        family_name = family_name,
        requested_markers = paste(requested, collapse = ","),
        available_markers = paste(available, collapse = ","),
        missing_markers = paste(missing_markers, collapse = ","),
        available_marker_count = length(available),
        detected_marker_count = detected_marker_count,
        any_detection_fraction = any_fraction,
        median_detected_markers_per_observation = as.numeric(median_hits),
        pseudobulk_sum = as.numeric(pseudobulk),
        mean_normalized_expression = as.numeric(mean_norm),
        comparative_score_role = "comparative_not_absence_gate"
      )
      cursor <- cursor + 1L
    }
  }
}

outdir <- normalizePath(a$out, mustWork = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
table_path <- file.path(outdir, "broad_family_evidence.tsv.gz")
manifest_path <- file.path(outdir, "broad_family_evidence_manifest.json")
fwrite(rbindlist(rows), table_path, sep = "\t")
manifest <- list(
  schema_version = "2.0",
  feature_scope = "full_feature",
  source_object = list(path = normalizePath(a$rds), sha256 = sha256_file(a$rds)),
  cluster_membership = list(path = normalizePath(a$clusters), sha256 = sha256_file(a$clusters)),
  biological_profile = list(path = normalizePath(a$profile), sha256 = sha256_file(a$profile)),
  candidate_catalog = list(path = normalizePath(a$catalog), sha256 = sha256_file(a$catalog)),
  evidence_table = list(path = normalizePath(table_path), sha256 = sha256_file(table_path)),
  validation_manifest = list(path = normalizePath(a$`validation-manifest`), sha256 = sha256_file(a$`validation-manifest`)),
  n_observations = ncol(obj),
  n_clusters = length(cluster_ids),
  candidate_ids = unique(candidate_ids),
  complete_cartesian_product = TRUE,
  decision_performed = FALSE,
  centered_scores_can_establish_absence = FALSE
)
write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(manifest, pretty = TRUE, auto_unbox = TRUE), "\n")
