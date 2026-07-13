#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: inspect_r_object.R INPUT.rds OUTPUT.json")
obj <- readRDS(args[[1]])
info <- list(path = normalizePath(args[[1]]), class = class(obj))
if (inherits(obj, "Seurat")) {
  suppressPackageStartupMessages({library(Seurat); library(SeuratObject)})
  info$type <- "Seurat"; info$n_features <- nrow(obj); info$n_observations <- ncol(obj)
  info$assays <- Assays(obj); info$reductions <- Reductions(obj); info$metadata_columns <- colnames(obj[[]])
  info$default_assay <- DefaultAssay(obj)
} else if (inherits(obj, "SingleCellExperiment") || inherits(obj, "SummarizedExperiment")) {
  suppressPackageStartupMessages({library(SingleCellExperiment); library(SummarizedExperiment)})
  info$type <- "SingleCellExperiment"; info$n_features <- nrow(obj); info$n_observations <- ncol(obj)
  info$assays <- assayNames(obj); info$reduced_dims <- reducedDimNames(obj); info$coldata_columns <- colnames(colData(obj))
  if (inherits(obj, "SpatialExperiment")) info$has_spatial_coords <- !is.null(spatialCoords(obj))
} else {
  info$type <- "other"; info$length <- length(obj); info$names <- names(obj)
}
write_json(info, args[[2]], pretty = TRUE, auto_unbox = TRUE, null = "null")
cat(args[[2]], "\n")

