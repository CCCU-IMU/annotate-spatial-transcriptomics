#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

parse_args <- function(x) {
  out <- list(); i <- 1L
  while (i <= length(x)) {
    key <- sub("^--", "", x[[i]])
    if (i == length(x) || startsWith(x[[i + 1L]], "--")) { out[[key]] <- TRUE; i <- i + 1L }
    else { out[[key]] <- x[[i + 1L]]; i <- i + 2L }
  }
  out
}
arg <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("rds", "metadata", "markers", "out", "cell-id-col", "broad-col", "subtype-col")
missing <- required[!required %in% names(arg)]
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))

dir.create(arg$out, recursive = TRUE, showWarnings = FALSE)
arg$out <- normalizePath(arg$out, mustWork = TRUE)
min_n <- as.integer(ifelse(is.null(arg$`min-n`), 20, arg$`min-n`))
clip <- as.numeric(ifelse(is.null(arg$clip), 2.5, arg$clip))
analysis_view <- ifelse(is.null(arg$`analysis-view`), "provided_labels", arg$`analysis-view`)
evidence_cohort <- ifelse(is.null(arg$`evidence-cohort`), "provided_metadata", arg$`evidence-cohort`)

read_any <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) fread(cmd = paste("gzip -dc", shQuote(path))) else fread(path)
}
safe_layer <- function(assay, layer) {
  tryCatch(LayerData(assay, layer = layer), error = function(e) GetAssayData(assay, slot = layer))
}
safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
wrap_marker_group <- function(values) vapply(values, function(value) {
  words <- strwrap(gsub("([/-])", "\\1 ", value), width = 18)
  paste(words, collapse = "\n")
}, character(1))
save_both <- function(p, stem, width, height) {
  ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = 420, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height, device = cairo_pdf, bg = "white", limitsize = FALSE)
}
hclust_segments <- function(hc) {
  n <- length(hc$order)
  if (n < 2L) return(data.table(x = numeric(), xend = numeric(), y = numeric(), yend = numeric()))
  leaf_x <- numeric(n); leaf_x[hc$order] <- seq_len(n)
  node_x <- node_y <- numeric(n - 1L); seg <- list(); z <- 1L
  child_xy <- function(child) if (child < 0) c(x = leaf_x[-child], y = 0) else c(x = node_x[child], y = node_y[child])
  for (i in seq_len(n - 1L)) {
    a <- child_xy(hc$merge[i, 1]); b <- child_xy(hc$merge[i, 2]); h <- hc$height[i]
    seg[[z]] <- data.table(x = a["x"], xend = a["x"], y = a["y"], yend = h); z <- z + 1L
    seg[[z]] <- data.table(x = b["x"], xend = b["x"], y = b["y"], yend = h); z <- z + 1L
    seg[[z]] <- data.table(x = a["x"], xend = b["x"], y = h, yend = h); z <- z + 1L
    node_x[i] <- mean(c(a["x"], b["x"])); node_y[i] <- h
  }
  rbindlist(seg)
}

message("Reading object and metadata")
obj <- readRDS(arg$rds)
meta <- read_any(arg$metadata)
cell_col <- arg$`cell-id-col`
stopifnot(cell_col %in% names(meta), arg$`broad-col` %in% names(meta), arg$`subtype-col` %in% names(meta))
meta[[cell_col]] <- as.character(meta[[cell_col]])
stopifnot(uniqueN(meta[[cell_col]]) == nrow(meta))
if(!is.null(arg$`include-col`)){
  if(!arg$`include-col`%in%names(meta))stop("include-col absent from metadata")
  allowed<-strsplit(arg$`include-values`,",",fixed=TRUE)[[1]]
  meta<-meta[tolower(as.character(get(arg$`include-col`)))%in%tolower(allowed)]
}
common <- intersect(colnames(obj), meta[[cell_col]])
if (!length(common)) stop("No overlapping cell IDs")
meta <- meta[match(common, get(cell_col))]
stopifnot(all(meta[[cell_col]] == common))

if (inherits(obj, "Seurat")) {
  assay_name <- if (!is.null(arg$assay)) arg$assay else DefaultAssay(obj)
  if (!assay_name %in% Assays(obj)) stop("Assay not present: ", assay_name)
  assay <- obj[[assay_name]]
  data_layer <- ifelse(is.null(arg$`data-layer`), "data", arg$`data-layer`)
  count_layer <- ifelse(is.null(arg$`count-layer`), "counts", arg$`count-layer`)
  data_mat <- safe_layer(assay, data_layer)[, common, drop = FALSE]
  count_mat <- safe_layer(assay, count_layer)[, common, drop = FALSE]
} else if (inherits(obj, "SingleCellExperiment") || inherits(obj, "SummarizedExperiment")) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("SummarizedExperiment package is required for SCE input")
  available <- SummarizedExperiment::assayNames(obj)
  data_name <- if (!is.null(arg$`data-assay`)) arg$`data-assay` else if ("normcounts" %in% available) "normcounts" else if ("logcounts" %in% available) "logcounts" else stop("Specify --data-assay")
  count_name <- if (!is.null(arg$`count-assay`)) arg$`count-assay` else if ("counts" %in% available) "counts" else stop("Specify --count-assay")
  data_mat <- SummarizedExperiment::assay(obj, data_name)[, common, drop = FALSE]
  count_mat <- SummarizedExperiment::assay(obj, count_name)[, common, drop = FALSE]
} else stop("Unsupported R object class: ", paste(class(obj), collapse = ","))
markers <- read_any(arg$markers)
stopifnot(all(c("gene", "marker_group") %in% names(markers)))
if (!"panel" %in% names(markers)) markers[, panel := "canonical"]
if (!"level" %in% names(markers)) markers[, level := "both"]
if (any(!markers$level %in% c("broad", "subtype", "both"))) stop("marker level must be broad, subtype or both")
markers <- unique(markers[gene %in% rownames(data_mat), .(gene, marker_group, panel, level)])
setorder(markers, panel, level, marker_group, gene)
# A marker may legitimately support more than one current label group. Keep the
# gene/group pairing so the report can facet by the biological label it supports.
markers <- markers[!duplicated(markers, by = c("panel", "level", "gene", "marker_group"))]
if (!nrow(markers)) stop("No marker genes overlap expression features")

make_source <- function(label_col, level_name, panel_name) {
  level_view <- if (!is.null(arg[[paste0(level_name, "-analysis-view")]])) arg[[paste0(level_name, "-analysis-view")]] else analysis_view
  level_cohort <- if (!is.null(arg[[paste0(level_name, "-evidence-cohort")]])) arg[[paste0(level_name, "-evidence-cohort")]] else evidence_cohort
  panel_dt <- markers[panel == panel_name & level %in% c(level_name, "both")]
  if (!nrow(panel_dt)) return(NULL)
  labels <- as.character(meta[[label_col]])
  valid <- !is.na(labels) & nzchar(labels)
  counts <- table(labels[valid]); keep_labels <- names(counts[counts >= min_n])
  valid <- valid & labels %in% keep_labels
  labels <- labels[valid]; cells <- common[valid]
  if (!length(cells) || length(unique(labels)) < 2L) stop(level_name, " has fewer than two eligible labels")
  genes <- unique(panel_dt[["gene"]])
  dm <- data_mat[genes, cells, drop = FALSE]; cm <- count_mat[genes, cells, drop = FALSE]
  rows <- vector("list", length(unique(labels))); z <- 1L
  for (lab in sort(unique(labels))) {
    use <- labels == lab
    rows[[z]] <- data.table(
      gene = genes, label = lab,
      avg_expression = as.numeric(Matrix::rowMeans(dm[, use, drop = FALSE])),
      pct_expressed_absolute = 100 * as.numeric(Matrix::rowMeans(cm[, use, drop = FALSE] > 0)),
      n_observations = sum(use)
    ); z <- z + 1L
  }
  ds <- merge(rbindlist(rows), panel_dt[, .(gene, marker_group)], by = "gene", all.x = TRUE, sort = FALSE)
  ds[!is.finite(avg_expression), avg_expression := 0]
  ds[!is.finite(pct_expressed_absolute), pct_expressed_absolute := 0]
  ds[, avg_expression_scaled_within_gene := {
    s <- sd(avg_expression)
    if (is.finite(s) && s > 0) pmax(pmin((avg_expression - mean(avg_expression)) / s, clip), -clip) else 0
  }, by = gene]
  ds[, pct_expressed_scaled_within_gene := {
    m <- max(pct_expressed_absolute)
    if (is.finite(m) && m > 0) pmax(pmin(100 * pct_expressed_absolute / m, 100), 0) else 0
  }, by = gene]
  ds[, `:=`(analysis_view=level_view,evidence_cohort=level_cohort)]

  tree_ds <- unique(ds[, .(label, gene, avg_expression_scaled_within_gene)])
  mat <- dcast(tree_ds, label ~ gene, value.var = "avg_expression_scaled_within_gene", fill = 0)
  labs <- mat$label; x <- as.matrix(mat[, -1]); rownames(x) <- labs
  hc <- hclust(dist(x), method = "ward.D2"); ord <- hc$labels[hc$order]
  ds[, label := factor(label, levels = ord)]
  ds[, gene := factor(gene, levels = unique(panel_dt[["gene"]]))]
  seg <- hclust_segments(hc)
  # Put the label tree on the left and marker genes on the x axis. Faceting
  # genes by marker_group makes the cell type/program supported by each marker
  # explicit above the x axis, while the tree and label rows remain aligned.
  p_tree <- ggplot(seg) +
    geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend), linewidth = 0.3) +
    scale_y_continuous(limits = c(0.5, length(ord) + 0.5), expand = c(0, 0)) +
    theme_void()
  p_dot <- ggplot(ds[pct_expressed_scaled_within_gene > 0], aes(gene, label)) +
    geom_point(aes(size = pct_expressed_scaled_within_gene, colour = avg_expression_scaled_within_gene), alpha = 0.9) +
    scale_size_area(max_size = 5.5, limits = c(0, 100), oob = scales::squish, breaks = c(25, 50, 75, 100), name = "Within-gene detection") +
    scale_colour_gradient2(low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426", midpoint = 0,
                           limits = c(-clip, clip), oob = scales::squish, name = "Within-gene expression") +
    facet_grid(. ~ marker_group, scales = "free_x", space = "free_x",
               labeller = labeller(marker_group = as_labeller(wrap_marker_group))) +
    scale_x_discrete(drop = TRUE) + scale_y_discrete(drop = FALSE) +
    labs(x = "Marker genes grouped by supported cell type/program", y = paste(level_name, "label"),
         subtitle = "Point size and color are normalized within each gene; absolute values are retained in the source TSV") +
    theme_bw(base_size = 7) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1), panel.grid = element_blank(),
          strip.text.x = element_text(face = "bold", size = 7))
  combined <- (p_tree | p_dot) + plot_layout(widths = c(1.2, 8))
  stem <- file.path(arg$out, paste0(level_name, "_", panel_name, "_marker_dotplot_with_tree"))
  fwrite(ds, paste0(stem, "_source.tsv"), sep = "\t", quote = FALSE)
  fwrite(data.table(tree_order = ord), paste0(stem, "_tree_order.tsv"), sep = "\t", quote = FALSE)
  save_both(combined, stem, max(12, 0.17 * length(unique(panel_dt[["gene"]])) + 7), max(7, 0.32 * length(ord) + 4))
  data.table(level = level_name, panel = panel_name, n_labels = length(ord), n_genes = length(unique(panel_dt[["gene"]])),
             analysis_view=level_view,evidence_cohort=level_cohort,
             png = paste0(stem, ".png"), pdf = paste0(stem, ".pdf"), source = paste0(stem, "_source.tsv"))
}

assets <- list(); z <- 1L
for (lev in c("broad", "subtype")) {
  col <- if (lev == "broad") arg$`broad-col` else arg$`subtype-col`
  for (panel_name in unique(markers$panel)) {
    message("Rendering ", lev, " / ", panel_name)
    ans <- make_source(col, lev, panel_name)
    if (!is.null(ans)) { assets[[z]] <- ans; z <- z + 1L }
  }
}
asset_dt <- rbindlist(assets, fill = TRUE)
fwrite(asset_dt, file.path(arg$out, "marker_dotplot_asset_index.tsv"), sep = "\t", quote = FALSE)
required_levels <- c("broad", "subtype")
if (!all(required_levels %in% asset_dt$level)) stop("Both broad and subtype dotplots are required")
message("Completed marker dotplots: ", nrow(asset_dt), " assets")
