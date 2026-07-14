# Shared fail-closed checks for Seurat marker/DEG evidence layers.

read_validation_manifest <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("A readable --validation-manifest is required for Spatial DEG/marker evidence. ",
         "Build it with prepare_seurat_full_feature_validation.R")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required")
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

hash_observation_ids <- function(ids) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required")
  digest::digest(paste(sort(as.character(ids)), collapse = "\n"), algo = "sha256")
}

layer_matrix <- function(object, assay, layer) {
  tryCatch(
    SeuratObject::LayerData(object[[assay]], layer = layer),
    error = function(e) SeuratObject::GetAssayData(object[[assay]], slot = layer)
  )
}

sparse_exact_equal <- function(left, right) {
  if (!all(dim(left) == dim(right))) return(FALSE)
  left <- methods::as(left, "dgCMatrix")
  right <- methods::as(right, "dgCMatrix")
  identical(left@p, right@p) && identical(left@i, right@i) &&
    identical(as.numeric(left@x), as.numeric(right@x))
}

assert_seurat_validation_layer <- function(object, assay, data_layer = "data",
                                           count_layer = "counts",
                                           manifest_path = NULL,
                                           object_path = NULL) {
  if (!assay %in% SeuratObject::Assays(object)) stop("Missing assay: ", assay)
  counts <- layer_matrix(object, assay, count_layer)
  data <- layer_matrix(object, assay, data_layer)
  if (!nrow(counts) || !ncol(counts) || !nrow(data) || !ncol(data)) {
    stop("Counts or validation data layer is empty")
  }
  if (!identical(rownames(counts), rownames(data)) || !identical(colnames(counts), colnames(data))) {
    stop("Counts and validation data layers have different feature/observation identities")
  }
  if (sparse_exact_equal(counts, data)) {
    stop("The selected Seurat data layer is exactly identical to raw counts and is not eligible ",
         "for Wilcoxon DEG. Build a separate full-feature LogNormalize validation object with ",
         "prepare_seurat_full_feature_validation.R; do not modify the SCT clustering object.")
  }
  data_sparse <- methods::as(data, "dgCMatrix")
  if (!length(data_sparse@x) || !any(abs(data_sparse@x - round(data_sparse@x)) > 1e-8)) {
    stop("The selected validation data layer lacks non-integer normalized values")
  }
  manifest <- NULL
  if (identical(tolower(assay), "spatial")) {
    manifest <- read_validation_manifest(manifest_path)
    required <- c("status", "role", "clustering_eligible", "normalization_method", "scale_factor", "assay",
                  "count_layer", "data_layer", "n_features", "n_observations",
                  "analysis_set_sha256", "counts_equal_data", "reductions_removed",
                  "normalized_object", "normalized_object_sha256")
    missing <- setdiff(required, names(manifest))
    if (length(missing)) stop("Validation manifest lacks fields: ", paste(missing, collapse = ", "))
    if (!identical(as.character(manifest$status), "PASS") ||
        !identical(as.character(manifest$role), "full_feature_deg_marker_validation_only") ||
        !identical(as.logical(manifest$clustering_eligible), FALSE) ||
        !identical(as.character(manifest$normalization_method), "LogNormalize") ||
        as.numeric(manifest$scale_factor) != 10000 ||
        !identical(as.character(manifest$assay), assay) ||
        !identical(as.character(manifest$count_layer), count_layer) ||
        !identical(as.character(manifest$data_layer), data_layer) ||
        as.integer(manifest$n_features) != nrow(counts) ||
        as.integer(manifest$n_observations) != ncol(counts) ||
        !identical(as.logical(manifest$counts_equal_data), FALSE) ||
        !identical(as.logical(manifest$reductions_removed), TRUE) ||
        !identical(as.character(manifest$analysis_set_sha256),
                   hash_observation_ids(colnames(counts)))) {
      stop("Validation manifest does not match the requested LogNormalize object/layers")
    }
    contract <- object@misc$full_feature_validation_contract
    if (is.null(contract) ||
        !identical(as.character(contract$role), as.character(manifest$role)) ||
        !identical(as.logical(contract$clustering_eligible), FALSE) ||
        !identical(as.character(contract$analysis_set_sha256), as.character(manifest$analysis_set_sha256)) ||
        !identical(as.character(contract$normalization_method), "LogNormalize") ||
        as.numeric(contract$scale_factor) != 10000 ||
        !identical(as.character(contract$assay), assay) ||
        !identical(as.character(contract$count_layer), count_layer) ||
        !identical(as.character(contract$data_layer), data_layer) ||
        length(SeuratObject::Reductions(object)) != 0L) {
      stop("Validation object-internal contract does not match its manifest")
    }
    if (!is.null(object_path) &&
        normalizePath(as.character(manifest$normalized_object), mustWork = TRUE) !=
          normalizePath(object_path, mustWork = TRUE)) {
      stop("Validation manifest points to a different normalized object")
    }
  }
  invisible(list(counts = counts, data = data, manifest = manifest))
}
