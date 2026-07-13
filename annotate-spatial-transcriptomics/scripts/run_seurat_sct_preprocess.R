≠rá^—f•ñÿ¶{ø,y 'v√Æ∂õ≠#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
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

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("rds", "out", "sample")
missing <- required[!required %in% names(args)]
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))

value_or <- function(name, default) if (is.null(args[[name]])) default else args[[name]]
as_num <- function(name, default) as.numeric(value_or(name, default))
as_int <- function(name, default) as.integer(value_or(name, default))

input_rds <- normalizePath(args$rds, mustWork = TRUE)
out_dir <- args$out
sample_id <- args$sample
assay_name <- value_or("assay", "Spatial")
count_layer <- value_or("count-layer", "counts")
min_counts <- as_num("min-counts", 100)
min_features <- as_num("min-features", 75)
nfeatures <- as_int("nfeatures", 3000)
sct_ncells_cap <- as_int("sct-ncells", 50000)
pca_npcs <- as_int("pca-npcs", 50)
dims_n <- as_int("dims", 30)
k_param <- as_int("k", 30)
annoy_trees <- as_int("annoy-trees", 50)
umap_min_dist <- as_num("umap-min-dist", 0.3)
seed <- as_int("seed", 20260706)
resolutions <- as.numeric(strsplit(value_or("resolutions", "0.1,0.2,0.3,0.4,0.6"), ",", fixed = TRUE)[[1]])
future_globals_max_gb <- as_num("future-globals-max-gb", 100)

if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
  stop("glmGamPoi is required by the frozen SCT batch profile; refusing a silent method fallback")
}
if (pca_npcs < dims_n) stop("--pca-npcs must be >= --dims")
if (any(!is.finite(resolutions)) || !length(resolutions)) stop("Invalid --resolutions")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(out_dir, "provenance"), showWarnings = FALSE)
options(future.globals.maxSize = future_globals_max_gb * 1024^3)
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")
set.seed(seed)

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("--rds is not a Seurat object")
if (!assay_name %in% Assays(obj)) stop("Assay not found: ", assay_name)
DefaultAssay(obj) <- assay_name

counts <- tryCatch(
  LayerData(obj[[assay_name]], layer = count_layer),
  error = function(e) GetAssayData(obj[[assay_name]], slot = count_layer)
)
if (!nrow(counts) || !ncol(counts)) stop("Selected count layer is empty")
n_counts <- Matrix::colSums(counts)
n_features <- Matrix::colSums(counts > 0)
obj[[paste0("nCount_", assay_name)]] <- n_counts[colnames(obj)]
obj[[paste0("nFeature_", assay_name)]] <- n_features[colnames(obj)]

count_col <- paste0("nCount_", assay_name)
feature_col <- paste0("nFeature_", assay_name)
metadata <- obj[[]]
qc_pass <- is.finite(metadata[[count_col]]) & is.finite(metadata[[feature_col]]) &
  metadata[[count_col]] >= min_counts & metadata[[feature_col]] >= min_features
high_count_cutoff <- as.numeric(quantile(metadata[[count_col]], 0.999, na.rm = TRUE))
high_count_flag <- metadata[[count_col]] > high_count_cutoff

scope <- data.table(
  cell_id = rownames(metadata),
  n_counts = as.numeric(metadata[[count_col]]),
  n_features = as.numeric(metadata[[feature_col]]),
  analysis_set = qc_pass,
  initial_scope = ifelse(qc_pass, "analysis_set", "excluded_initial_qc"),
  high_count_review_flag = high_count_flag
)
fwrite(scope, file.path(out_dir, "tables", "analysis_scope.tsv.gz"), sep = "\t")
analysis_ids <- scope[analysis_set == TRUE, cell_id]
if (length(analysis_ids) < 3L) stop("Fewer than three observations passed entry QC")
obj <- subset(obj, cells = analysis_ids)
rm(counts); invisible(gc())

sct_ncells <- min(sct_ncells_cap, ncol(obj))
obj <- SCTransform(
  obj,
  assay = assay_name,
  new.assay.name = "SCT",
  vst.flavor = "v2",
  method = "glmGamPoi",
  variable.features.n = nfeatures,
  ncells = sct_ncells,
  conserve.memory = TRUE,
  return.only.var.genes = TRUE,
  seed.use = seed,
  verbose = TRUE
)
DefaultAssay(obj) <- "SCT"
obj <- RunPCA(
  obj,
  assay = "SCT",
  npcs = pca_npcs,
  features = VariableFeatures(obj),
  seed.use = seed,
  verbose = TRUE
)
obj <- FindNeighbors(
  obj,
  reduction = "pca",
  dims = seq_len(dims_n),
  k.param = k_param,
  nn.method = "annoy",
  n.trees = annoy_trees,
  annoy.metric = "cosine",
  graph.name = c("SCT_nn", "SCT_snn"),
  verbose = TRUE
)
obj <- FindClusters(
  obj,
  graph.name = "SCT_snn",
  algorithm = 4,
  resolution = resolutions,
  random.seed = seed,
  verbose = TRUE
)
obj <- RunUMAP(
  obj,
  reduction = "pca",
  dims = seq_len(dims_n),
  n.neighbors = k_param,
  min.dist = umap_min_dist,
  metric = "cosine",
  seed.use = seed,
  verbose = TRUE
)

resolution_cols <- paste0("SCT_snn_res.", resolutions)
if (!all(resolution_cols %in% colnames(obj[[]]))) {
  stop("Missing resolution columns: ", paste(setdiff(resolution_cols, colnames(obj[[]])), collapse = ", "))
}
for (i in seq_along(resolutions)) {
  resolution <- resolutions[[i]]
  cluster_col <- resolution_cols[[i]]
  labels <- as.character(obj[[cluster_col, drop = TRUE]])
  tag <- gsub("\\.", "p", as.character(resolution))
  membership <- data.table(cell_id = colnames(obj), resolution = resolution, cluster = labels)
  counts_table <- membership[, .(n_observations = .N), by = .(resolution, cluster)]
  counts_table[, `:=`(
    fraction = n_observations / sum(n_observations),
    small_cluster_review = n_observations < 100L
  )]
  fwrite(membership, file.path(out_dir, "tables", paste0("clusters_res", tag, ".tsv.gz")), sep = "\t")
  fwrite(counts_table, file.path(out_dir, "tables", paste0("cluster_counts_res", tag, ".tsv")), sep = "\t")
}

umap <- Embeddings(obj, "umap")
fwrite(
  data.table(cell_id = rownames(umap), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2]),
  file.path(out_dir, "tables", "umap_coordinates.tsv.gz"), sep = "\t"
)

metadata_names <- colnames(obj[[]])
coordinate_candidates <- list(c("x", "y"), c("sdimx", "sdimy"), c("imagecol", "imagerow"), c("col", "row"))
coordinate_hits <- vapply(coordinate_candidates, function(pair) all(pair %in% metadata_names), logical(1))
spatial_pair <- if (any(coordinate_hits)) coordinate_candidates[[which(coordinate_hits)[[1]]]] else character()
if (length(spatial_pair)) {
  md <- obj[[]]
  fwrite(
    data.table(
      cell_id = rownames(md),
      x = as.numeric(md[[spatial_pair[[1]]]]),
      y = as.numeric(md[[spatial_pair[[2]]]])
    ),
    file.path(out_dir, "tables", "spatial_coordinates.tsv.gz"), sep = "\t"
  )
}

hash_file <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) digest::digest(path, algo = "sha256", file = TRUE) else NA_character_
}
hash_ids <- function(ids) {
  ids <- sort(as.character(ids))
  if (requireNamespace("digest", quietly = TRUE)) digest::digest(paste(ids, collapse = "\n"), algo = "sha256") else NA_character_
}
manifest <- data.table(
  parameter = c(
    "sample", "input_rds", "input_sha256", "assay", "count_layer",
    "entry_qc", "high_count_policy", "mitochondrial_hard_filter", "doublet_hard_filter",
    "n_input", "n_analysis", "analysis_set_sha256", "sct_vst_flavor", "sct_method",
    "sct_variable_features", "sct_ncells", "sct_conserve_memory", "sct_return_only_var_genes",
    "pca_npcs", "neighbor_dims", "neighbor_k", "neighbor_method", "neighbor_trees",
    "neighbor_metric", "cluster_algorithm", "resolution_grid", "resolution_selection",
    "umap_neighbors", "umap_min_dist", "umap_metric", "seed", "future_plan",
    "tiny_cluster_policy", "spatial_x_col", "spatial_y_col"
  ),
  value = as.character(c(
    sample_id, input_rds, hash_file(input_rds), assay_name, count_layer,
    paste0(count_col, ">=", min_counts, " AND ", feature_col, ">=", min_features),
    paste0("flag_only_above_q99.9=", signif(high_count_cutoff, 8)), "none", "none",
    nrow(scope), ncol(obj), hash_ids(colnames(obj)), "v2", "glmGamPoi",
    nfeatures, sct_ncells, TRUE, TRUE, pca_npcs, dims_n, k_param, "annoy", annoy_trees,
    "cosine", 4, paste(resolutions, collapse = ","), "not_selected_adaptive_review_required",
    k_param, umap_min_dist, "cosine", seed, "sequential",
    "flag_below_100_never_auto_reassign", if (length(spatial_pair)) spatial_pair[[1]] else "",
    if (length(spatial_pair)) spatial_pair[[2]] else ""
  ))
)
fwrite(manifest, file.path(out_dir, "provenance", "preprocessing_manifest.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(out_dir, "provenance", "sessionInfo.txt"))
saveRDS(obj, file.path(out_dir, paste0(sample_id, "_seurat_sct_preprocessed.rds")), compress = FALSE)
writeLines(c("status\tPASS", paste0("completed_at\t", format(Sys.time(), tz = "UTC", usetz = TRUE))), file.path(out_dir, "RUN_COMPLETE.tsv"))
