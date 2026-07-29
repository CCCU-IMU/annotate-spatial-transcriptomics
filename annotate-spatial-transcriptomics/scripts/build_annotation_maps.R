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
    } else {
      out[[key]] <- x[[i + 1L]]; i <- i + 2L
    }
  }
  out
}

a <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("rds", "metadata", "out", "cell-id-col", "final-cell-type-col")
missing <- required[!required %in% names(a)]
if (length(missing)) stop("Missing: ", paste(missing, collapse = ", "))
dir.create(a$out, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(a$out, "figures")
node_dir <- file.path(a$out, "spatial_nodes")
table_dir <- file.path(a$out, "tables")
dir.create(figure_dir, showWarnings = FALSE)
dir.create(node_dir, showWarnings = FALSE)
dir.create(table_dir, showWarnings = FALSE)

read_any <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    fread(cmd = paste("gzip -dc", shQuote(path)))
  } else fread(path)
}
save_both <- function(plot, stem, width = 9, height = 7, background = "white") {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
         dpi = 400, bg = background, limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
         device = cairo_pdf, bg = background, limitsize = FALSE)
}
safe <- function(x) substr(gsub("_+", "_", gsub("[^A-Za-z0-9_.-]+", "_", x)), 1, 150)

preferred <- c(
  "Oocyte" = "#FFD60A", "Granulosa" = "#FF375F", "Theca" = "#FF9F0A",
  "Luteal" = "#FF453A", "Stromal/mesenchymal" = "#64D2FF",
  "Smooth muscle" = "#30D158", "Endothelial" = "#00E5FF",
  "Lymphatic endothelial" = "#00A6A6", "Pericyte/mural" = "#BF5AF2",
  "Immune" = "#5E5CE6", "Epithelial/mesothelial" = "#FF2D55",
  "Glial/Schwann-like" = "#AC8E68", "QC/Unknown" = "#8E8E93"
)
palette_for <- function(levels) {
  missing_levels <- setdiff(levels, names(preferred))
  extra <- if (length(missing_levels)) {
    setNames(hcl.colors(length(missing_levels), "Dark 3"), missing_levels)
  } else character()
  c(preferred[intersect(names(preferred), levels)], extra)[levels]
}

object <- readRDS(a$rds)
metadata <- read_any(a$metadata)
cell_col <- a$`cell-id-col`; type_col <- a$`final-cell-type-col`
if (!all(c(cell_col, type_col) %in% names(metadata))) {
  stop("metadata lacks cell ID or final_cell_type")
}
metadata[[cell_col]] <- as.character(metadata[[cell_col]])
if (uniqueN(metadata[[cell_col]]) != nrow(metadata)) stop("duplicate metadata IDs")
cells <- intersect(colnames(object), metadata[[cell_col]])
if (!length(cells)) stop("No matching IDs")
metadata <- metadata[match(cells, get(cell_col))]
labels <- as.character(metadata[[type_col]])
labels[is.na(labels) | !nzchar(labels)] <- "QC/Unknown"

plot_map <- function(data, title, stem, spatial = FALSE) {
  data[, label := labels]
  levels <- sort(unique(data$label))
  plot <- ggplot(data, aes(x, y, colour = label)) +
    scattermore::geom_scattermore(pointsize = if (spatial) .82 else .75,
                                  pixels = c(2200, 2200)) +
    scale_colour_manual(values = palette_for(levels)) + labs(title = title, colour = NULL)
  if (spatial) {
    plot <- plot + scale_y_reverse() + coord_equal() + theme_void() +
      theme(plot.background = element_rect(fill = "black", colour = NA),
            panel.background = element_rect(fill = "black", colour = NA),
            legend.background = element_rect(fill = "black", colour = NA),
            legend.key = element_rect(fill = "black", colour = NA),
            legend.text = element_text(colour = "white"),
            plot.title = element_text(colour = "white"))
    save_both(plot, stem, background = "black")
  } else {
    plot <- plot + theme_classic(base_size = 8)
    save_both(plot, stem)
  }
}

if (is.null(a$`skip-umap`)) {
  if (!is.null(a$umap)) {
    umap <- read_any(a$umap); umap[[cell_col]] <- as.character(umap[[cell_col]])
    umap <- umap[match(cells, get(cell_col))]
    ux <- intersect(c("UMAP_1", "umap_1", "UMAP1"), names(umap))[[1L]]
    uy <- intersect(c("UMAP_2", "umap_2", "UMAP2"), names(umap))[[1L]]
    umap_data <- data.table(cell_id = cells, x = umap[[ux]], y = umap[[uy]])
  } else {
    reduction <- ifelse(is.null(a$reduction), "umap", a$reduction)
    embedding <- Embeddings(object, reduction)[cells, 1:2, drop = FALSE]
    umap_data <- data.table(cell_id = cells, x = embedding[, 1], y = embedding[, 2])
  }
  plot_map(umap_data, "Final cell type UMAP",
           file.path(figure_dir, "final_cell_type_UMAP"), FALSE)
}

spatial <- NULL
if (!is.null(a$coordinates)) {
  coordinates <- read_any(a$coordinates)
  coordinates[[cell_col]] <- as.character(coordinates[[cell_col]])
  coordinates <- coordinates[match(cells, get(cell_col))]
  x_col <- intersect(c("sdimx", "x", "spatial_x"), names(coordinates))[[1L]]
  y_col <- intersect(c("sdimy", "y", "spatial_y"), names(coordinates))[[1L]]
  spatial <- data.table(cell_id = cells, x = coordinates[[x_col]], y = coordinates[[y_col]])
} else if (inherits(object, "Seurat") && all(c("x", "y") %in% colnames(object[[]]))) {
  md <- object[[]]
  spatial <- data.table(cell_id = cells, x = as.numeric(md[cells, "x"]),
                        y = as.numeric(md[cells, "y"]))
}

node_rows <- list()
if (!is.null(spatial)) {
  plot_map(spatial, "Final cell type spatial",
           file.path(figure_dir, "final_cell_type_spatial"), TRUE)
  for (label in sort(unique(labels))) {
    selected <- labels == label
    plot <- ggplot(spatial, aes(x, y)) +
      scattermore::geom_scattermore(data = spatial[!selected], colour = "#595959",
                                    pointsize = .55, pixels = c(1600, 1600)) +
      scattermore::geom_scattermore(data = spatial[selected], colour = "#FF2D20",
                                    pointsize = .9, pixels = c(1600, 1600)) +
      scale_y_reverse() + coord_equal() + theme_void() +
      labs(title = paste0(label, " (n=", sum(selected), ")")) +
      theme(plot.background = element_rect(fill = "black", colour = NA),
            panel.background = element_rect(fill = "black", colour = NA),
            plot.title = element_text(colour = "white"))
    stem <- file.path(node_dir, paste0("cell_type__", safe(label)))
    save_both(plot, stem, background = "black")
    node_rows[[length(node_rows) + 1L]] <- data.table(
      level = "cell_type", parent_label = "", label = label,
      n_observations = sum(selected), png = paste0(stem, ".png"),
      pdf = paste0(stem, ".pdf")
    )
  }
}
if (length(node_rows)) {
  fwrite(rbindlist(node_rows), file.path(table_dir, "spatial_node_asset_index.tsv"),
         sep = "\t")
}
