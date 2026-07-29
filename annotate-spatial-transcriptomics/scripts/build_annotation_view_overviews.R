#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(data.table)
  library(ggplot2); library(scattermore)
})
parse_args <- function(x) {
  out <- list(); i <- 1L
  while (i <= length(x)) {
    key <- sub("^--", "", x[[i]])
    if (i == length(x) || startsWith(x[[i + 1L]], "--")) {
      out[[key]] <- TRUE; i <- i + 1L
    } else {out[[key]] <- x[[i + 1L]]; i <- i + 2L}
  }
  out
}
a <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("rds", "metadata", "out", "cell-id-col", "final-cell-type-col")
missing <- required[!required %in% names(a)]
if (length(missing)) stop("Missing: ", paste(missing, collapse = ", "))
dir.create(file.path(a$out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(a$out, "tables"), recursive = TRUE, showWarnings = FALSE)
read_any <- function(path) if (grepl("\\.gz$", path)) {
  fread(cmd = paste("gzip -dc", shQuote(path)))
} else fread(path)
save_both <- function(plot, stem) {
  ggsave(paste0(stem, ".png"), plot, width = 9, height = 7, dpi = 400,
         bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = 9, height = 7,
         device = cairo_pdf, bg = "white", limitsize = FALSE)
}
object <- readRDS(a$rds); metadata <- read_any(a$metadata)
cell_col <- a$`cell-id-col`; type_col <- a$`final-cell-type-col`
if (!all(c(cell_col, type_col) %in% names(metadata))) stop("missing final_cell_type")
metadata[[cell_col]] <- as.character(metadata[[cell_col]])
if (uniqueN(metadata[[cell_col]]) != nrow(metadata)) stop("duplicate metadata IDs")
cells <- intersect(colnames(object), metadata[[cell_col]])
metadata <- metadata[match(cells, get(cell_col))]
if ("analysis_scope" %in% names(metadata)) {
  keep <- metadata$analysis_scope == "analysis_set"
  cells <- cells[keep]; metadata <- metadata[keep]
}
labels <- as.character(metadata[[type_col]])
labels[is.na(labels) | !nzchar(labels)] <- "QC/Unknown"
if (!is.null(a$umap)) {
  umap <- read_any(a$umap); umap[[cell_col]] <- as.character(umap[[cell_col]])
  umap <- umap[match(cells, get(cell_col))]
  ux <- intersect(c("UMAP_1", "umap_1", "UMAP1"), names(umap))[[1L]]
  uy <- intersect(c("UMAP_2", "umap_2", "UMAP2"), names(umap))[[1L]]
  umap_data <- data.table(x = umap[[ux]], y = umap[[uy]], label = labels)
} else {
  reduction <- ifelse(is.null(a$reduction), "umap", a$reduction)
  embedding <- Embeddings(object, reduction)[cells, 1:2, drop = FALSE]
  umap_data <- data.table(x = embedding[, 1], y = embedding[, 2], label = labels)
}
plot_one <- function(data, title, stem, spatial = FALSE) {
  palette <- setNames(hcl.colors(length(sort(unique(data$label))), "Dark 3"),
                      sort(unique(data$label)))
  plot <- ggplot(data, aes(x, y, colour = label)) +
    scattermore::geom_scattermore(pointsize = ifelse(spatial, .5, .65),
                                  pixels = c(2200, 2200)) +
    scale_colour_manual(values = palette) + labs(title = title, colour = NULL)
  if (spatial) plot <- plot + scale_y_reverse() + coord_equal() + theme_void()
  else plot <- plot + theme_classic(base_size = 8)
  save_both(plot, stem)
}
umap_stem <- file.path(a$out, "figures", "final_cell_type_UMAP")
plot_one(umap_data, "Final cell type UMAP", umap_stem, FALSE)
spatial_stem <- ""
if (!is.null(a$coordinates)) {
  coordinates <- read_any(a$coordinates)
  coordinates[[cell_col]] <- as.character(coordinates[[cell_col]])
  coordinates <- coordinates[match(cells, get(cell_col))]
  x_col <- intersect(c("sdimx", "x", "spatial_x"), names(coordinates))[[1L]]
  y_col <- intersect(c("sdimy", "y", "spatial_y"), names(coordinates))[[1L]]
  spatial <- data.table(x = coordinates[[x_col]], y = coordinates[[y_col]], label = labels)
  spatial_stem <- file.path(a$out, "figures", "final_cell_type_spatial")
  plot_one(spatial, "Final cell type spatial", spatial_stem, TRUE)
}
fwrite(data.table(
  view = "final", level = "cell_type", n_labeled = length(labels),
  umap_png = paste0(umap_stem, ".png"), umap_pdf = paste0(umap_stem, ".pdf"),
  spatial_png = if (nzchar(spatial_stem)) paste0(spatial_stem, ".png") else "",
  spatial_pdf = if (nzchar(spatial_stem)) paste0(spatial_stem, ".pdf") else ""
), file.path(a$out, "tables", "final_annotation_overview_asset_index.tsv"), sep = "\t")
