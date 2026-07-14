#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scattermore)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot resolve script directory")
source(file.path(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))),
                 "seurat_validation_layer.R"))

parse_args <- function(x) {
  out <- list(); i <- 1L
  while (i <= length(x)) {
    key <- sub("^--", "", x[[i]])
    if (i == length(x) || startsWith(x[[i + 1L]], "--")) { out[[key]] <- TRUE; i <- i + 1L }
    else { out[[key]] <- x[[i + 1L]]; i <- i + 2L }
  }; out
}
a <- parse_args(commandArgs(trailingOnly = TRUE))
need <- c("rds", "clusters", "out", "cell-id-col", "cluster-col")
miss <- need[!need %in% names(a)]; if (length(miss)) stop("Missing: ", paste(miss, collapse = ", "))
dir.create(a$out, recursive = TRUE, showWarnings = FALSE)
tab <- file.path(a$out, "tables"); fig <- file.path(a$out, "figures"); dir.create(tab, showWarnings = FALSE); dir.create(fig, showWarnings = FALSE)
min_pct <- as.numeric(ifelse(is.null(a$`min-pct`), 0.05, a$`min-pct`))
lfc <- as.numeric(ifelse(is.null(a$`logfc-threshold`), 0.1, a$`logfc-threshold`))
max_cells <- as.integer(ifelse(is.null(a$`max-cells-per-ident`), 3000, a$`max-cells-per-ident`))
seed <- as.integer(ifelse(is.null(a$seed), 20260713, a$seed)); set.seed(seed)

read_any <- function(path) if (grepl("\\.gz$", path, ignore.case = TRUE)) fread(cmd = paste("gzip -dc", shQuote(path))) else fread(path)
save_both <- function(p, stem, w, h) {
  ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 420, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h, device = cairo_pdf, bg = "white", limitsize = FALSE)
}
message("Reading Seurat object")
obj <- readRDS(a$rds); stopifnot(inherits(obj, "Seurat"))
cl <- read_any(a$clusters); cl[[a$`cell-id-col`]] <- as.character(cl[[a$`cell-id-col`]])
stopifnot(uniqueN(cl[[a$`cell-id-col`]]) == nrow(cl))
common <- intersect(colnames(obj), cl[[a$`cell-id-col`]])
if (length(common) != ncol(obj) || length(common) != nrow(cl)) stop("Object/cluster IDs do not match exactly")
cl <- cl[match(colnames(obj), get(a$`cell-id-col`))]; stopifnot(all(cl[[a$`cell-id-col`]] == colnames(obj)))
obj$framework_cluster <- factor(as.character(cl[[a$`cluster-col`]]), levels = sort(unique(as.character(cl[[a$`cluster-col`]]))))
assay <- if (!is.null(a$assay)) a$assay else DefaultAssay(obj)
data_layer <- if (!is.null(a$layer)) a$layer else "data"
count_layer <- if (!is.null(a$`count-layer`)) a$`count-layer` else "counts"
if (!assay %in% Assays(obj)) stop("Missing assay: ", assay)
DefaultAssay(obj) <- assay
assert_seurat_validation_layer(
  obj, assay = assay, data_layer = data_layer, count_layer = count_layer,
  manifest_path = a$`validation-manifest`, object_path = a$rds
)
Idents(obj) <- "framework_cluster"

counts <- data.table(cluster = levels(obj$framework_cluster), n_observations = as.integer(table(obj$framework_cluster)))
fwrite(counts, file.path(tab, "cluster_counts.tsv"), sep = "\t")

message("Running one-vs-rest DEG")
deg_path <- file.path(tab, "cluster_DEG_one_vs_rest_all.tsv")
sig_path <- file.path(tab, "cluster_DEG_significant_positive.tsv")
top_path <- file.path(tab, "cluster_DEG_top100_positive.tsv")
if (all(file.exists(c(deg_path, sig_path, top_path)))) {
  message("Reusing completed DEG tables")
  deg <- fread(deg_path); sig <- fread(sig_path)
} else {
  deg <- as.data.table(FindAllMarkers(obj, assay = assay, slot = data_layer, only.pos = TRUE, test.use = "wilcox",
                                      min.pct = min_pct, logfc.threshold = lfc, return.thresh = 1,
                                      max.cells.per.ident = max_cells, random.seed = seed, densify = FALSE, verbose = TRUE))
  lfc_col <- intersect(c("avg_log2FC", "avg_logFC"), names(deg))[1]
  fwrite(deg, deg_path, sep = "\t", quote = FALSE)
  sig <- deg[p_val_adj < 0.05 & get(lfc_col) > 0]
  fwrite(sig, sig_path, sep = "\t", quote = FALSE)
  fwrite(sig[order(cluster, p_val_adj, -get(lfc_col)), head(.SD, 100), by = cluster], top_path, sep = "\t", quote = FALSE)
}

pal <- setNames(hcl.colors(length(levels(obj$framework_cluster)), "Dynamic"), levels(obj$framework_cluster))
if (all(c("sdimx", "sdimy") %in% names(cl))) {
  sp <- data.table(x = as.numeric(cl$sdimx), y = as.numeric(cl$sdimy), cluster = obj$framework_cluster)
} else if (all(c("x", "y") %in% colnames(obj[[]]))) {
  sp <- data.table(x = as.numeric(obj$x), y = as.numeric(obj$y), cluster = obj$framework_cluster)
} else stop("No spatial coordinates in cluster table or object metadata")
p_sp <- ggplot(sp, aes(x, y, colour = cluster)) + scattermore::geom_scattermore(pointsize = 0.45, pixels = c(2400, 2400)) +
  scale_colour_manual(values = pal) + scale_y_reverse() + coord_equal() + theme_void() + labs(title = "Selected clustering: spatial", colour = "Cluster")
save_both(p_sp, file.path(fig, "selected_clustering_spatial"), 10.5, 8)

if (!is.null(a$umap) && file.exists(a$umap)) {
  um <- read_any(a$umap); um[[a$`cell-id-col`]] <- as.character(um[[a$`cell-id-col`]])
  um <- um[match(colnames(obj), get(a$`cell-id-col`))]
  ux <- intersect(c("UMAP_1", "umap_1", "UMAP1"), names(um))[1]; uy <- intersect(c("UMAP_2", "umap_2", "UMAP2"), names(um))[1]
  if (!is.na(ux) && !is.na(uy)) {
    ud <- data.table(x = um[[ux]], y = um[[uy]], cluster = obj$framework_cluster)
    p_um <- ggplot(ud, aes(x, y, colour = cluster)) + scattermore::geom_scattermore(pointsize = 0.6, pixels = c(2200, 2200)) +
      scale_colour_manual(values = pal) + theme_classic(base_size = 8) + labs(title = "Selected clustering: UMAP", colour = "Cluster")
    save_both(p_um, file.path(fig, "selected_clustering_UMAP"), 10, 8)
  }
}

message("Generating per-cluster spatial highlight grid")
plots <- lapply(levels(obj$framework_cluster), function(k) {
  d <- copy(sp); d[, selected := cluster == k]
  ggplot(d, aes(x, y)) +
    scattermore::geom_scattermore(data = d[selected == FALSE], colour = "#E2E2E2", pointsize = 0.35, pixels = c(1100, 1100)) +
    scattermore::geom_scattermore(data = d[selected == TRUE], colour = "#D62728", pointsize = 0.75, pixels = c(1100, 1100)) +
    scale_y_reverse() + coord_equal() + theme_void() + labs(title = paste0("c", k, " (n=", sum(d$selected), ")")) +
    theme(plot.title = element_text(size = 7, hjust = 0.5))
})
grid <- wrap_plots(plots, ncol = 4) + plot_annotation(title = "Per-cluster spatial highlights")
save_both(grid, file.path(fig, "selected_clustering_highlight_grid"), 13, max(5, ceiling(length(plots) / 4) * 3.1))

manifest <- data.table(
  parameter = c("n_observations", "n_clusters", "assay", "count_layer", "data_layer",
                "validation_manifest", "min_pct", "logfc_threshold", "max_cells_per_ident", "seed"),
  value = c(ncol(obj), length(levels(obj$framework_cluster)), assay, count_layer, data_layer,
            ifelse(is.null(a$`validation-manifest`), "", normalizePath(a$`validation-manifest`)),
            min_pct, lfc, max_cells, seed)
)
fwrite(manifest, file.path(a$out, "evidence_run_manifest.tsv"), sep = "\t")
message("Initial cluster evidence complete")
