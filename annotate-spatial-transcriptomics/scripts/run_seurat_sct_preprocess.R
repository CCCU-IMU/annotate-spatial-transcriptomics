#!/usr/bin/env Rscript

# Compatibility entry point retained for existing workflows. The verified
# same-batch whole-tissue implementation is now SCT v2 + BANKSY, not Seurat SNN.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(Banksy)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
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

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("rds", "out", "sample")
missing <- required[!required %in% names(args)]
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))

retired <- intersect(names(args), c(
  "pca-npcs", "dims", "k", "annoy-trees",
  "resolution-workers", "resolution-future-plan"
))
if (length(retired)) {
  stop(
    "Retired Seurat-SNN whole-tissue arguments supplied: ", paste(retired, collapse = ", "),
    ". Use the fixed SCT+BANKSY parameters or an explicit BANKSY batch exception."
  )
}

value_or <- function(name, default) if (is.null(args[[name]])) default else args[[name]]
as_num <- function(name, default) as.numeric(value_or(name, default))
as_int <- function(name, default) as.integer(value_or(name, default))
fmt <- function(x) format(x, scientific = FALSE, trim = TRUE)

input_rds <- normalizePath(args$rds, mustWork = TRUE)
out_dir <- args$out
sample_id <- args$sample
assay_name <- value_or("assay", "Spatial")
count_layer <- value_or("count-layer", "counts")
min_counts <- as_num("min-counts", 100)
min_features <- as_num("min-features", 75)
nfeatures <- as_int("nfeatures", 4000)
sct_ncells_cap <- as_int("sct-ncells", 50000)
banksy_npcs <- as_int("banksy-npcs", 30)
k_geom <- as_int("banksy-k-geom", 30)
lambda <- as_num("banksy-lambda", 0.2)
k_neighbors <- as_int("banksy-k-neighbors", 50)
resolutions <- as.numeric(strsplit(value_or("resolutions", "0.2,0.4,0.6,0.8"), ",", fixed = TRUE)[[1]])
umap_neighbors <- as_int("umap-neighbors", 30)
umap_min_dist <- as_num("umap-min-dist", 0.3)
umap_spread <- as_num("umap-spread", 1)
umap_epochs <- as_int("umap-epochs", 300)
seed <- as_int("seed", 20260717)
future_globals_max_gb <- as_num("future-globals-max-gb", 700)

detect_scheduler_cpus <- function() {
  for (name in c("LSB_DJOB_NUMPROC", "SLURM_CPUS_PER_TASK", "NSLOTS", "AIP_CPUS")) {
    value <- suppressWarnings(as.integer(Sys.getenv(name, unset = "")))
    if (length(value) && is.finite(value) && value > 0L) return(list(cpus = value, source = name))
  }
  list(cpus = NA_integer_, source = "not_detected")
}
scheduler_cpu <- detect_scheduler_cpus()
analysis_threads_default <- if (is.finite(scheduler_cpu$cpus)) scheduler_cpu$cpus else 1L
analysis_threads <- as_int("analysis-threads", analysis_threads_default)

fixed_profile <- list(
  `min-counts` = 100, `min-features` = 75, nfeatures = 4000,
  `sct-ncells` = 50000, `banksy-npcs` = 30, `banksy-k-geom` = 30,
  `banksy-lambda` = 0.2, `banksy-k-neighbors` = 50,
  resolutions = "0.2,0.4,0.6,0.8", `umap-neighbors` = 30,
  `umap-min-dist` = 0.3, `umap-spread` = 1, `umap-epochs` = 300,
  seed = 20260717
)
observed_profile <- list(
  `min-counts` = min_counts, `min-features` = min_features, nfeatures = nfeatures,
  `sct-ncells` = sct_ncells_cap, `banksy-npcs` = banksy_npcs,
  `banksy-k-geom` = k_geom, `banksy-lambda` = lambda,
  `banksy-k-neighbors` = k_neighbors,
  resolutions = paste(resolutions, collapse = ","),
  `umap-neighbors` = umap_neighbors, `umap-min-dist` = umap_min_dist,
  `umap-spread` = umap_spread, `umap-epochs` = umap_epochs, seed = seed
)
different <- names(fixed_profile)[vapply(names(fixed_profile), function(name) {
  !identical(as.character(fixed_profile[[name]]), as.character(observed_profile[[name]]))
}, logical(1))]
exception_allowed <- isTRUE(args$`allow-batch-exception`)
exception_reason <- as.character(value_or("batch-exception-reason", ""))
if (length(different) && (!exception_allowed || nchar(trimws(exception_reason)) < 20L)) {
  stop(
    "Frozen StereoPy cellbin_PPed SCT+BANKSY parameters changed (",
    paste(different, collapse = ", "),
    "). Use --allow-batch-exception with --batch-exception-reason of at least 20 characters."
  )
}

for (package in c("glmGamPoi", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required by the frozen SCT+BANKSY batch profile")
  }
}
if (any(!is.finite(resolutions)) || !length(resolutions) || anyDuplicated(resolutions)) {
  stop("Invalid --resolutions")
}
if (banksy_npcs < 2L || k_geom < 1L || k_neighbors < 2L) stop("Invalid BANKSY graph parameters")
if (!is.finite(analysis_threads) || analysis_threads < 1L) stop("--analysis-threads must be >= 1")
if (is.finite(scheduler_cpu$cpus) && analysis_threads > scheduler_cpu$cpus) {
  stop("--analysis-threads exceeds scheduler allocation from ", scheduler_cpu$source)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(out_dir, "provenance"), showWarnings = FALSE)
options(future.globals.maxSize = future_globals_max_gb * 1024^3)
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")
set.seed(seed)

source_obj <- readRDS(input_rds)
if (!inherits(source_obj, "Seurat")) stop("--rds is not a Seurat object")
if (!assay_name %in% names(source_obj@assays)) stop("Assay not found: ", assay_name)
counts <- tryCatch(
  LayerData(source_obj[[assay_name]], layer = count_layer),
  error = function(e) GetAssayData(source_obj[[assay_name]], slot = count_layer)
)
if (!nrow(counts) || !ncol(counts)) stop("Selected count layer is empty")
if (!inherits(counts, "sparseMatrix") || any(counts@x < 0) ||
    any(abs(counts@x - round(counts@x)) > 1e-10)) {
  stop("The selected input must be a nonnegative sparse integer count matrix")
}

cell_ids <- colnames(counts)
metadata <- source_obj@meta.data[cell_ids, , drop = FALSE]
coordinate_candidates <- list(c("x", "y"), c("sdimx", "sdimy"), c("imagecol", "imagerow"), c("col", "row"))
coordinate_hits <- vapply(coordinate_candidates, function(pair) all(pair %in% colnames(metadata)), logical(1))
if (any(coordinate_hits)) {
  spatial_pair <- coordinate_candidates[[which(coordinate_hits)[[1]]]]
  coords <- as.matrix(metadata[, spatial_pair, drop = FALSE])
} else if ("spatial" %in% Reductions(source_obj)) {
  spatial_pair <- c("spatial_1", "spatial_2")
  coords <- Embeddings(source_obj, "spatial")[cell_ids, 1:2, drop = FALSE]
} else {
  stop("No x/y-like metadata or spatial reduction found")
}
storage.mode(coords) <- "double"
rownames(coords) <- cell_ids
colnames(coords) <- c("x", "y")

n_counts <- Matrix::colSums(counts)
n_features <- Matrix::colSums(counts > 0)
qc_pass <- is.finite(n_counts) & is.finite(n_features) &
  n_counts >= min_counts & n_features >= min_features &
  is.finite(coords[, 1]) & is.finite(coords[, 2])
high_count_cutoff <- as.numeric(quantile(n_counts, 0.999, na.rm = TRUE))
scope <- data.table(
  cell_id = cell_ids,
  n_counts = as.numeric(n_counts),
  n_features = as.numeric(n_features),
  analysis_set = qc_pass,
  initial_scope = ifelse(qc_pass, "analysis_set", "excluded_initial_qc"),
  high_count_review_flag = n_counts > high_count_cutoff
)
fwrite(scope, file.path(out_dir, "tables", "analysis_scope.tsv.gz"), sep = "\t")
analysis_ids <- scope[analysis_set == TRUE, cell_id]
if (length(analysis_ids) < 3L) stop("Fewer than three observations passed entry QC")

clean_meta <- data.frame(
  sample_id = sample_id,
  x = coords[analysis_ids, 1],
  y = coords[analysis_ids, 2],
  nCount_source = as.numeric(n_counts[analysis_ids]),
  nFeature_source = as.numeric(n_features[analysis_ids]),
  qc_pass = TRUE,
  high_count_review_flag = n_counts[analysis_ids] > high_count_cutoff,
  row.names = analysis_ids,
  check.names = FALSE
)
obj <- CreateSeuratObject(
  counts = counts[, analysis_ids, drop = FALSE], assay = assay_name,
  project = paste0(sample_id, "_SCT_BANKSY"), min.cells = 1,
  min.features = 0, meta.data = clean_meta
)
rm(source_obj, counts, metadata, n_counts, n_features); invisible(gc())

spatial_embedding <- as.matrix(obj@meta.data[, c("x", "y"), drop = FALSE])
colnames(spatial_embedding) <- c("coord_1", "coord_2")
obj[["spatial"]] <- CreateDimReducObject(
  embeddings = spatial_embedding, key = "coord_", assay = assay_name
)

# Keep all-gene LogNormalize expression for marker access, but never feed it to BANKSY.
obj <- NormalizeData(
  obj, assay = assay_name, normalization.method = "LogNormalize",
  scale.factor = 10000, verbose = TRUE
)
sct_ncells <- min(sct_ncells_cap, ncol(obj))
obj <- SCTransform(
  obj, assay = assay_name, new.assay.name = "SCT",
  vst.flavor = "v2", method = "glmGamPoi",
  variable.features.n = nfeatures, ncells = sct_ncells,
  conserve.memory = TRUE, return.only.var.genes = TRUE,
  seed.use = seed, verbose = TRUE
)
banksy_features <- VariableFeatures(obj, assay = "SCT")
banksy_matrix <- LayerData(obj[["SCT"]], layer = "scale.data")[banksy_features, , drop = FALSE]
if (nrow(banksy_matrix) != nfeatures || ncol(banksy_matrix) != ncol(obj) ||
    !identical(colnames(banksy_matrix), colnames(obj))) {
  stop("SCT Pearson-residual BANKSY matrix failed dimensions or identity checks")
}

spe <- SpatialExperiment(
  assays = SimpleList(banksy_input = banksy_matrix),
  spatialCoords = as.matrix(obj@meta.data[, c("x", "y"), drop = FALSE])
)
colnames(spe) <- colnames(obj)
rm(banksy_matrix); invisible(gc())

spe <- computeBanksy(
  spe, assay_name = "banksy_input", compute_agf = FALSE,
  k_geom = k_geom, spatial_mode = "kNN_median", seed = seed,
  dimensions = "all", center = TRUE, verbose = TRUE
)
spe <- runBanksyPCA(
  spe, use_agf = FALSE, lambda = lambda, npcs = banksy_npcs,
  assay_name = "banksy_input", scale = TRUE, seed = seed
)
spe <- clusterBanksy(
  spe, use_agf = FALSE, lambda = lambda, use_pcs = TRUE,
  npcs = banksy_npcs, assay_name = "banksy_input", algo = "leiden",
  k_neighbors = k_neighbors, resolution = resolutions, seed = seed
)
spe <- runBanksyUMAP(
  spe, use_agf = FALSE, lambda = lambda, use_pcs = TRUE,
  npcs = banksy_npcs, assay_name = "banksy_input",
  n_neighbors = umap_neighbors, spread = umap_spread,
  min_dist = umap_min_dist, n_epochs = umap_epochs,
  seed = seed, n_threads = analysis_threads
)

pca_name <- grep(paste0("^PCA_M0_lam", gsub("\\.", "\\\\.", fmt(lambda)), "$"), reducedDimNames(spe), value = TRUE)
if (length(pca_name) != 1L) pca_name <- grep("^PCA", reducedDimNames(spe), value = TRUE)[1]
umap_name <- grep(paste0("^UMAP_M0_lam", gsub("\\.", "\\\\.", fmt(lambda)), "$"), reducedDimNames(spe), value = TRUE)
if (length(umap_name) != 1L) umap_name <- grep("^UMAP", reducedDimNames(spe), value = TRUE)[1]
if (length(pca_name) != 1L || length(umap_name) != 1L) stop("BANKSY PCA/UMAP reductions not found")

pca_embedding <- as.matrix(reducedDim(spe, pca_name))
umap_embedding <- as.matrix(reducedDim(spe, umap_name))
rownames(pca_embedding) <- rownames(umap_embedding) <- colnames(obj)
colnames(pca_embedding) <- paste0("BANKSYPC_", seq_len(ncol(pca_embedding)))
colnames(umap_embedding) <- paste0("BANKSYUMAP_", seq_len(ncol(umap_embedding)))
obj[["banksy_pca"]] <- CreateDimReducObject(
  embeddings = pca_embedding, key = "BANKSYPC_", assay = "SCT"
)
obj[["banksy_umap"]] <- CreateDimReducObject(
  embeddings = umap_embedding, key = "BANKSYUMAP_", assay = "SCT"
)

cluster_cols <- grep("^clust_", colnames(colData(spe)), value = TRUE)
if (length(cluster_cols) != length(resolutions)) {
  stop("BANKSY cluster grid does not match the requested resolution grid")
}
for (i in seq_along(resolutions)) {
  resolution <- resolutions[[i]]
  hits <- cluster_cols[endsWith(cluster_cols, paste0("res", fmt(resolution)))]
  if (length(hits) != 1L) stop("Cannot bind BANKSY cluster column for resolution ", resolution)
  cluster_col <- hits[[1]]
  output_col <- paste0("banksy_", cluster_col)
  labels <- as.character(colData(spe)[[cluster_col]])
  obj[[output_col]] <- as.factor(labels)
  tag <- gsub("\\.", "p", fmt(resolution))
  membership <- data.table(cell_id = colnames(obj), resolution = resolution, cluster = labels)
  counts_table <- membership[, .(n_observations = .N), by = .(resolution, cluster)]
  counts_table[, `:=`(
    fraction = n_observations / sum(n_observations),
    small_cluster_review = n_observations < 100L
  )]
  fwrite(membership, file.path(out_dir, "tables", paste0("clusters_res", tag, ".tsv.gz")), sep = "\t")
  fwrite(counts_table, file.path(out_dir, "tables", paste0("cluster_counts_res", tag, ".tsv")), sep = "\t")
}

obj$banksy_preprocess_method <- "sct"
obj$banksy_input_source <- "SCT v2 glmGamPoi Pearson residuals, 4000 HVG"
obj$banksy_compute_agf <- FALSE
obj$banksy_k_geom <- k_geom
obj$banksy_lambda <- lambda
obj$banksy_k_neighbors <- k_neighbors
DefaultAssay(obj) <- "SCT"

fwrite(
  data.table(cell_id = rownames(umap_embedding), UMAP_1 = umap_embedding[, 1], UMAP_2 = umap_embedding[, 2]),
  file.path(out_dir, "tables", "umap_coordinates.tsv.gz"), sep = "\t"
)
md <- obj[[]]
fwrite(
  data.table(cell_id = rownames(md), x = as.numeric(md$x), y = as.numeric(md$y)),
  file.path(out_dir, "tables", "spatial_coordinates.tsv.gz"), sep = "\t"
)

hash_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
hash_ids <- function(ids) digest::digest(paste(sort(as.character(ids)), collapse = "\n"), algo = "sha256")
grid_artifact <- list(
  method = "SCT+BANKSY",
  candidate_resolutions = as.list(resolutions),
  cluster_columns = paste0("banksy_", cluster_cols),
  parameters = list(M = 0, compute_agf = FALSE, k_geom = k_geom,
                    spatial_mode = "kNN_median", lambda = lambda,
                    npcs = banksy_npcs, k_neighbors = k_neighbors)
)
write_json(grid_artifact, file.path(out_dir, "provenance", "banksy_grid.json"), pretty = TRUE, auto_unbox = TRUE)

parameter_manifest <- list(
  sample_id = sample_id,
  input_rds = input_rds,
  input_sha256 = hash_file(input_rds),
  source_policy = "raw Spatial/counts and x/y only; imported normalization/PCA/UMAP/clusters/labels ignored",
  analysis_set_sha256 = hash_ids(colnames(obj)),
  qc = list(min_counts = min_counts, min_features = min_features,
            high_count_policy = paste0("flag_only_above_q99.9=", signif(high_count_cutoff, 8)),
            mitochondrial_hard_filter = FALSE, doublet_hard_filter = FALSE,
            retained_cells = ncol(obj), retained_features = nrow(obj)),
  normalization = list(primary = "SCTransform", vst_flavor = "v2", method = "glmGamPoi",
                       ncells = sct_ncells, n_hvg = nfeatures, conserve_memory = TRUE,
                       return_only_var_genes = TRUE,
                       auxiliary_marker_layer = paste0(assay_name, " LogNormalize scale.factor=10000; not used by BANKSY")),
  banksy = list(version = as.character(packageVersion("Banksy")),
                input = "SCT scale.data: v2 glmGamPoi Pearson residuals, 4000 HVG",
                compute_agf = FALSE, M = 0, k_geom = k_geom,
                spatial_mode = "kNN_median", lambda = lambda, npcs = banksy_npcs,
                k_neighbors = k_neighbors, resolution_grid = as.list(resolutions),
                umap = list(n_neighbors = umap_neighbors, min_dist = umap_min_dist,
                            spread = umap_spread, n_epochs = umap_epochs)),
  seed = seed,
  batch_profile_status = if (length(different)) "approved_recorded_exception" else "frozen_profile_exact",
  batch_exception_reason = if (length(different)) exception_reason else "",
  output_reductions = c("spatial", "banksy_pca", "banksy_umap"),
  output_cluster_columns = paste0("banksy_", cluster_cols)
)
write_json(parameter_manifest, file.path(out_dir, "provenance", "preprocessing_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
fwrite(
  data.table(
    parameter = c(
      "sample", "input_rds", "input_sha256", "source_policy", "analysis_set_sha256",
      "entry_qc", "high_count_policy", "mitochondrial_hard_filter", "doublet_hard_filter",
      "sct_vst_flavor", "sct_method", "sct_variable_features", "sct_ncells",
      "sct_conserve_memory", "sct_return_only_var_genes", "banksy_input",
      "banksy_compute_agf", "banksy_M", "banksy_k_geom", "banksy_spatial_mode",
      "banksy_lambda", "banksy_npcs", "banksy_k_neighbors", "resolution_grid",
      "resolution_selection", "umap_neighbors", "umap_min_dist", "umap_spread",
      "umap_epochs", "seed", "analysis_threads", "scheduler_cpus_detected",
      "scheduler_cpu_source", "batch_profile_status", "batch_exception_reason"
    ),
    value = as.character(c(
      sample_id, input_rds, hash_file(input_rds), parameter_manifest$source_policy,
      hash_ids(colnames(obj)), paste0("nCount_", assay_name, ">=", min_counts,
      " AND nFeature_", assay_name, ">=", min_features),
      parameter_manifest$qc$high_count_policy, "none", "none", "v2", "glmGamPoi",
      nfeatures, sct_ncells, TRUE, TRUE,
      "SCT_scale.data_Pearson_residuals_4000_HVG", FALSE, 0, k_geom,
      "kNN_median", lambda, banksy_npcs, k_neighbors, paste(resolutions, collapse = ","),
      "not_selected_adaptive_review_required", umap_neighbors, umap_min_dist,
      umap_spread, umap_epochs, seed, analysis_threads,
      if (is.finite(scheduler_cpu$cpus)) scheduler_cpu$cpus else "", scheduler_cpu$source,
      parameter_manifest$batch_profile_status, parameter_manifest$batch_exception_reason
    ))
  ),
  file.path(out_dir, "provenance", "preprocessing_manifest.tsv"), sep = "\t"
)
writeLines(banksy_features, file.path(out_dir, "provenance", "banksy_features.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "provenance", "sessionInfo.txt"))

output_rds <- file.path(out_dir, paste0(sample_id, "_sct_BANKSY_preprocessed_seurat.rds"))
saveRDS(obj, output_rds, compress = "gzip")
writeLines(
  c("status\tPASS", paste0("completed_at\t", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("output_rds\t", output_rds)),
  file.path(out_dir, "RUN_COMPLETE.tsv")
)
