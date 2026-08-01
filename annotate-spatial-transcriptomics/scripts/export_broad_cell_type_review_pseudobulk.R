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
values <- function(flag) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(character())
  if (any(hit == length(args))) stop("missing value for ", flag)
  unique(args[hit + 1])
}
read_character_table <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    connection <- gzfile(path, open = "rt")
    on.exit(close(connection), add = TRUE)
    return(as.data.table(read.delim(
      connection, header = TRUE, sep = "\t", quote = "", comment.char = "",
      check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character"
    )))
  }
  fread(path, colClasses = "character")
}

rds_path <- value("--rds")
count_cache_path <- value("--count-cache", "")
raw_count_assay <- value("--raw-count-assay", "")
membership_path <- value("--membership")
recall_path <- value("--recall-membership")
out_dir <- value("--out")
requested_assay <- value("--assay", "")
active_broad <- value("--active-broad-label", "")
comparison_broads <- values("--comparison-broad-label")
if (any(vapply(list(membership_path, recall_path, out_dir), is.null, logical(1))) ||
    (!nzchar(count_cache_path) && is.null(rds_path))) {
  stop("required: (--count-cache or --rds) --membership --recall-membership --out")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (nzchar(count_cache_path)) {
  if (!file.exists(count_cache_path)) stop("raw-count cache is absent")
  counts <- readRDS(count_cache_path)
  if (!inherits(counts, "sparseMatrix")) stop("raw-count cache is not sparse")
  assay <- raw_count_assay
  if (!nzchar(assay) || toupper(assay) == "SCT") {
    stop("raw-count cache lacks a non-SCT assay identity")
  }
} else {
  obj <- readRDS(rds_path)
  available <- Assays(obj)
  if (nzchar(requested_assay)) {
    if (!requested_assay %in% available) stop("requested assay is absent")
    assay <- requested_assay
  } else if ("RNA" %in% available) {
    assay <- "RNA"
  } else if ("Spatial" %in% available) {
    assay <- "Spatial"
  } else stop("no RNA/Spatial raw-count assay")
  if (assay == "SCT") stop("SCT corrected counts cannot be pseudobulk input")
  counts <- tryCatch(
    GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "counts")
  )
}
membership <- read_character_table(membership_path)
recall <- read_character_table(recall_path)
if (!all(c("cell_id", "final_broad_label") %in% names(membership))) stop("membership lacks broad label")
if (!all(c("cell_id", "broad_label") %in% names(recall))) stop("recall membership lacks target broad")
if (anyDuplicated(membership$cell_id)) stop("membership duplicates cell_id")
if (length(setdiff(membership$cell_id, colnames(counts))) > 0) stop("membership has foreign cells")

groups <- list()
present_broads <- sort(unique(
  membership$final_broad_label[nzchar(membership$final_broad_label)]
))
if (nzchar(active_broad)) {
  present_broads <- intersect(
    present_broads, unique(c(active_broad, comparison_broads))
  )
}
for (broad in present_broads) {
  groups[[paste0("current::", broad)]] <- unique(membership[final_broad_label == broad, cell_id])
  groups[[paste0("outside_current::", broad)]] <- unique(membership[final_broad_label != broad, cell_id])
}
for (broad in sort(unique(recall$broad_label[nzchar(recall$broad_label)]))) {
  if (nzchar(active_broad) && broad != active_broad) next
  target_ids <- unique(recall[broad_label == broad, cell_id])
  groups[[paste0("recall_question::", broad)]] <- target_ids
  origins <- membership[cell_id %in% target_ids, .(cell_id, final_broad_label)]
  for (origin in sort(unique(origins$final_broad_label[nzchar(origins$final_broad_label)]))) {
    origin_ids <- origins[final_broad_label == origin, cell_id]
    if (length(origin_ids) >= 5) {
      groups[[paste0("recall_question::", broad, "::from::", origin)]] <- origin_ids
    }
  }
}
groups[["whole_analysis_set"]] <- membership$cell_id

groups <- lapply(groups, function(ids) unique(ids[ids %in% colnames(counts)]))
groups <- groups[lengths(groups) > 0]
cell_index <- setNames(seq_len(ncol(counts)), colnames(counts))
design_i <- unlist(lapply(groups, function(ids) unname(cell_index[ids])), use.names = FALSE)
design_j <- rep(seq_along(groups), lengths(groups))
design <- sparseMatrix(
  i = design_i, j = design_j, x = 1,
  dims = c(ncol(counts), length(groups)),
  dimnames = list(colnames(counts), names(groups))
)
# One sparse multiplication computes every current-label and recall-question
# pseudobulk. This avoids repeatedly slicing the same 12 GB object once per
# cell type or origin label.
sums_matrix <- counts %*% design
detection_matrix <- (counts > 0) %*% design
group_n <- as.numeric(Matrix::colSums(design))
group_total <- as.numeric(Matrix::colSums(sums_matrix))
features <- rownames(counts)
rows <- lapply(seq_along(groups), function(index) {
  sums <- as.numeric(sums_matrix[, index])
  data.table(
    gene = features, group_id = names(groups)[index],
    n_observations = group_n[index], sum_counts = sums,
    cpm = if (group_total[index] > 0) 1e6 * sums / group_total[index] else 0,
    detection_fraction = as.numeric(detection_matrix[, index]) / group_n[index]
  )
})
census <- lapply(seq_along(groups), function(index) data.table(
  group_id = names(groups)[index], n_observations = group_n[index],
  total_raw_counts = group_total[index]
))
pseudobulk_path <- file.path(out_dir, "broad_cell_type_review_pseudobulk.tsv.gz")
census_path <- file.path(out_dir, "broad_cell_type_review_pseudobulk_census.tsv")
fwrite(rbindlist(rows), pseudobulk_path, sep = "\t", compress = "gzip")
fwrite(rbindlist(census), census_path, sep = "\t")
manifest <- list(
  schema_version = "2.2",
  status = "PASS_EVIDENCE_ONLY",
  artifact_role = "broad_cell_type_full_transcriptome_pseudobulk",
  formal_membership_written = FALSE,
  source_rds = if (!is.null(rds_path)) normalizePath(rds_path) else "",
  source_count_cache = if (nzchar(count_cache_path)) normalizePath(count_cache_path) else "",
  source_membership = normalizePath(membership_path),
  recall_membership = normalizePath(recall_path),
  raw_count_assay = assay,
  assay_ancestry = "project_local_non_SCT_raw_counts",
  feature_n = nrow(counts),
  analysis_observation_n = nrow(membership),
  group_n = length(rows),
  active_broad_label = active_broad,
  comparison_broad_labels = comparison_broads,
  active_broad_only = nzchar(active_broad),
  pseudobulk = normalizePath(pseudobulk_path),
  census = normalizePath(census_path)
)
write_json(manifest, file.path(out_dir, "broad_cell_type_review_pseudobulk_manifest.json"),
           pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(manifest, pretty = TRUE, auto_unbox = TRUE), "\n")
