#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit[length(hit)] == length(args)) stop("missing value for ", flag)
  args[hit[length(hit)] + 1]
}

rds_path <- value("--rds")
membership_path <- value("--analysis-membership")
marker_path <- value("--marker-manifest")
out_dir <- value("--out")
requested_assay <- value("--assay", "")
if (any(vapply(list(rds_path, membership_path, marker_path, out_dir), is.null, logical(1)))) {
  stop("required: --rds --analysis-membership --marker-manifest --out")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(rds_path)
available_assays <- Assays(obj)
if (nzchar(requested_assay)) {
  if (!requested_assay %in% available_assays) stop("requested assay is absent: ", requested_assay)
  assay <- requested_assay
} else if ("RNA" %in% available_assays) {
  assay <- "RNA"
} else if ("Spatial" %in% available_assays) {
  assay <- "Spatial"
} else {
  stop("no project-local RNA or Spatial raw-count assay is available")
}
if (assay == "SCT") stop("SCT corrected counts cannot be used for cell-type review")

counts <- tryCatch(
  GetAssayData(obj, assay = assay, layer = "counts"),
  error = function(e) GetAssayData(obj, assay = assay, slot = "counts")
)
if (nrow(counts) == 0 || ncol(counts) == 0) stop("selected raw-count assay is empty")

metadata <- obj[[]]
coordinate_candidates <- list(c("x", "y"), c("X", "Y"), c("imagecol", "imagerow"))
coordinate_columns <- NULL
for (candidate in coordinate_candidates) {
  if (all(candidate %in% colnames(metadata))) {
    coordinate_columns <- candidate
    break
  }
}
if (is.null(coordinate_columns)) {
  stop("cell-type review requires project-local x/y coordinates in Seurat metadata")
}

read_character_table <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    connection <- gzfile(path, open = "rt")
    on.exit(close(connection), add = TRUE)
    return(as.data.table(read.delim(
      connection, header = TRUE, sep = "\t", quote = "", comment.char = "",
      check.names = FALSE, stringsAsFactors = FALSE,
      colClasses = "character"
    )))
  }
  fread(path, colClasses = "character")
}

membership <- read_character_table(membership_path)
if (!"cell_id" %in% names(membership)) stop("analysis membership lacks cell_id")
if (anyDuplicated(membership$cell_id)) stop("analysis membership duplicates cell_id")
missing_cells <- setdiff(membership$cell_id, colnames(counts))
if (length(missing_cells) > 0) stop("analysis membership contains cells absent from raw counts")
counts <- counts[, membership$cell_id, drop = FALSE]
coordinates <- data.table(
  cell_id = membership$cell_id,
  x = as.numeric(metadata[membership$cell_id, coordinate_columns[[1]]]),
  y = as.numeric(metadata[membership$cell_id, coordinate_columns[[2]]])
)
if (any(!is.finite(coordinates$x)) || any(!is.finite(coordinates$y))) {
  stop("cell-type review coordinates contain non-finite values")
}

marker_manifest <- unique(fread(marker_path, colClasses = "character"))
if (!all(c("gene", "candidate_id", "broad_label", "evidence_role", "family_id") %in% names(marker_manifest))) {
  stop("marker manifest lacks required columns")
}
features <- rownames(counts)
feature_by_upper <- split(features, toupper(features))
requested_genes <- sort(unique(marker_manifest$gene))
matched <- vapply(requested_genes, function(gene) {
  values <- feature_by_upper[[toupper(gene)]]
  if (is.null(values)) "" else sort(values)[1]
}, character(1))
gene_map <- data.table(
  requested_gene = requested_genes,
  matched_feature = unname(matched),
  status = ifelse(nzchar(matched), "matched", "missing")
)
present <- gene_map[nzchar(matched_feature)]
if (nrow(present) == 0) stop("none of the catalog review markers occur in the raw-count assay")
marker_counts <- counts[present$matched_feature, , drop = FALSE]
rownames(marker_counts) <- present$requested_gene

# The full sparse raw-count matrix is immutable across the serial per-broad
# review.  Persist it once so later DEG/pseudobulk steps do not deserialize the
# complete Seurat object (assays, reductions and graphs) for every cell type.
full_count_path <- file.path(out_dir, "cell_type_review_full_raw_counts.rds")
saveRDS(counts, full_count_path, compress = FALSE)

matrix_path <- file.path(out_dir, "cell_type_review_marker_counts.mtx")
writeMM(marker_counts, matrix_path)
system2("gzip", c("-f", matrix_path))
fwrite(gene_map, file.path(out_dir, "cell_type_review_gene_map.tsv"), sep = "\t")
fwrite(data.table(cell_index = seq_along(colnames(marker_counts)), cell_id = colnames(marker_counts)),
       file.path(out_dir, "cell_type_review_cells.tsv"), sep = "\t")
fwrite(data.table(
  cell_id = colnames(counts),
  total_raw_counts = as.numeric(Matrix::colSums(counts)),
  detected_raw_genes = as.numeric(Matrix::colSums(counts > 0))
), file.path(out_dir, "cell_type_review_library_size.tsv.gz"), sep = "\t")
fwrite(coordinates, file.path(out_dir, "cell_type_review_coordinates.tsv"), sep = "\t")

manifest <- list(
  schema_version = "2.2",
  artifact_role = "query_raw_count_cell_type_review_export",
  source_rds = normalizePath(rds_path),
  raw_count_assay = assay,
  assay_ancestry = "project_local_non_SCT_raw_counts",
  analysis_observation_n = ncol(counts),
  full_feature_n = nrow(counts),
  requested_marker_n = length(requested_genes),
  matched_marker_n = nrow(present),
  missing_marker_n = sum(!nzchar(matched)),
  marker_matrix = normalizePath(paste0(matrix_path, ".gz")),
  full_count_matrix = normalizePath(full_count_path),
  gene_map = normalizePath(file.path(out_dir, "cell_type_review_gene_map.tsv")),
  cells = normalizePath(file.path(out_dir, "cell_type_review_cells.tsv")),
  library_size = normalizePath(file.path(out_dir, "cell_type_review_library_size.tsv.gz")),
  coordinates = normalizePath(file.path(out_dir, "cell_type_review_coordinates.tsv")),
  coordinate_columns = coordinate_columns
)
write_json(manifest, file.path(out_dir, "cell_type_review_count_export_manifest.json"),
           pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(manifest, pretty = TRUE, auto_unbox = TRUE), "\n")
