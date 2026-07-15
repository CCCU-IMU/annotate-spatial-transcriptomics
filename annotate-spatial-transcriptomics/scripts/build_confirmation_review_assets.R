#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(data.table)
  library(ggplot2); library(scattermore); library(jsonlite)
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
read_any <- function(path) if (grepl("\\.gz$", path, ignore.case=TRUE)) fread(cmd=paste("gzip -dc", shQuote(path))) else fread(path)
safe_layer <- function(assay, layer) tryCatch(LayerData(assay, layer=layer), error=function(e) GetAssayData(assay, slot=layer))
file_sha256 <- function(path) {
  ans <- system2("sha256sum", shQuote(path), stdout=TRUE, stderr=TRUE)
  status <- attr(ans, "status")
  if (!length(ans) || (!is.null(status) && status != 0)) stop("sha256sum failed for ", path)
  strsplit(ans[[1]], "[[:space:]]+")[[1]][[1]]
}

a <- parse_args(commandArgs(trailingOnly=TRUE))
required <- c("project-root", "rds", "metadata", "markers", "out", "cell-id-col", "broad-col")
missing <- required[!required %in% names(a)]
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse=", "))
if (a$`broad-col` != "primary_broad_label") stop("confirmation assets must use --broad-col primary_broad_label")
root <- normalizePath(a$`project-root`, mustWork=TRUE)
out <- normalizePath(a$out, mustWork=FALSE); dir.create(out, recursive=TRUE, showWarnings=FALSE); out <- normalizePath(out)
obj <- readRDS(a$rds); meta <- read_any(a$metadata); markers <- read_any(a$markers)
if (!inherits(obj, "Seurat")) stop("lightweight confirmation asset builder currently requires a Seurat RDS")
cell_col <- a$`cell-id-col`; broad_col <- a$`broad-col`
if (!all(c(cell_col, broad_col) %in% names(meta))) stop("metadata lacks cell/broad columns")
meta[[cell_col]] <- as.character(meta[[cell_col]])
if (anyDuplicated(meta[[cell_col]])) stop("metadata contains duplicate observation IDs")
cells <- intersect(colnames(obj), meta[[cell_col]])
if (!length(cells)) stop("no metadata IDs overlap the R object")
meta <- meta[match(cells, get(cell_col))]
labels <- as.character(meta[[broad_col]]); broad_labels <- sort(unique(labels[!is.na(labels) & nzchar(labels)]))
if (!length(broad_labels)) stop("no accepted broad labels in primary_broad_label")

palette_bank <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#332288",
  "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677",
  "#882255", "#AA4499", "#661100", "#6699CC", "#EE6677", "#228833",
  "#CCBB44", "#66CCEE", "#AA3377", "#4477AA", "#EE8866", "#BBBBBB"
)
if (length(broad_labels) > length(palette_bank)) stop("too many broad labels for the audited review palette; simplify taxonomy")
palette <- setNames(palette_bank[seq_along(broad_labels)], broad_labels)
palette_path <- file.path(out, "broad_palette.tsv")
fwrite(data.table(label=broad_labels, color=unname(palette)), palette_path, sep="\t")

coordinates <- NULL
if (!is.null(a$coordinates)) {
  coordinates <- read_any(a$coordinates); coordinates[[cell_col]] <- as.character(coordinates[[cell_col]])
  coordinates <- coordinates[match(cells, get(cell_col))]
} else if (all(c("x", "y") %in% names(meta))) {
  coordinates <- copy(meta)
} else if (inherits(obj, "Seurat") && all(c("x", "y") %in% colnames(obj[[]]))) {
  md <- obj[[]]; coordinates <- data.table(cell_id=cells, x=as.numeric(md[cells,"x"]), y=as.numeric(md[cells,"y"]))
  setnames(coordinates, "cell_id", cell_col)
} else stop("spatial coordinates are required for the pre-confirmation spatial review")
x_candidates <- intersect(c("sdimx","x","spatial_x"), names(coordinates))
y_candidates <- intersect(c("sdimy","y","spatial_y"), names(coordinates))
x_col <- if (is.null(a$`x-col`) && length(x_candidates)) x_candidates[[1]] else a$`x-col`
y_col <- if (is.null(a$`y-col`) && length(y_candidates)) y_candidates[[1]] else a$`y-col`
if (is.null(x_col) || is.null(y_col) || is.na(x_col) || is.na(y_col) || !all(c(x_col,y_col) %in% names(coordinates))) stop("cannot resolve spatial x/y columns")
spatial <- data.table(x=as.numeric(coordinates[[x_col]]), y=as.numeric(coordinates[[y_col]]), label=labels)
spatial[, accepted := !is.na(label) & nzchar(label)]
p_spatial <- ggplot(spatial, aes(x,y)) +
  scattermore::geom_scattermore(data=spatial[accepted==FALSE], colour="#D9D9D9", pointsize=.35, pixels=c(2200,2200)) +
  scattermore::geom_scattermore(data=spatial[accepted==TRUE], aes(colour=label), pointsize=.62, pixels=c(2200,2200)) +
  scale_colour_manual(values=palette, drop=FALSE) + scale_y_reverse() + coord_equal() + theme_void(base_size=9) +
  theme(legend.position="right", legend.key.height=grid::unit(.45,"cm"), legend.text=element_text(size=8)) +
  guides(colour=guide_legend(override.aes=list(size=3))) + labs(title="Frozen broad annotation for confirmation", colour="Broad class")
spatial_png <- file.path(out, "broad_spatial_review.png")
ggsave(spatial_png, p_spatial, width=12, height=9, dpi=320, bg="white", limitsize=FALSE)

if (!all(c("gene","marker_group") %in% names(markers))) stop("marker table requires gene and marker_group")
if ("panel" %in% names(markers)) markers <- markers[panel=="canonical"]
if ("level" %in% names(markers)) markers <- markers[level %in% c("broad","both")]
markers <- unique(markers[marker_group %in% broad_labels, .(gene, marker_group)])
missing_groups <- setdiff(broad_labels, unique(markers$marker_group))
if (length(missing_groups)) stop("canonical marker groups missing for broad labels: ", paste(missing_groups, collapse=", "))
assay_name <- ifelse(is.null(a$assay), DefaultAssay(obj), a$assay)
if (!assay_name %in% Assays(obj)) stop("assay missing: ", assay_name)
data_layer <- ifelse(is.null(a$`data-layer`), "data", a$`data-layer`)
count_layer <- ifelse(is.null(a$`count-layer`), "counts", a$`count-layer`)
data_mat <- safe_layer(obj[[assay_name]], data_layer)
count_mat <- safe_layer(obj[[assay_name]], count_layer)
markers <- markers[gene %in% rownames(data_mat)]
missing_groups <- setdiff(broad_labels, unique(markers$marker_group))
if (length(missing_groups)) stop("no expressed canonical marker remains for: ", paste(missing_groups, collapse=", "))
genes <- unique(markers$gene); valid <- !is.na(labels) & nzchar(labels); dot_cells <- cells[valid]; dot_labels <- labels[valid]
rows <- lapply(broad_labels, function(label) {
  use <- dot_labels == label
  data.table(gene=genes, label=label,
             avg_expression=as.numeric(Matrix::rowMeans(data_mat[genes,dot_cells[use],drop=FALSE])),
             pct_expressed_absolute=100*as.numeric(Matrix::rowMeans(count_mat[genes,dot_cells[use],drop=FALSE] > 0)),
             n_observations=sum(use))
})
dot <- merge(rbindlist(rows), markers, by="gene", all.x=TRUE, allow.cartesian=TRUE)
dot[, avg_expression_scaled_within_gene := { s<-sd(avg_expression); if(is.finite(s)&&s>0) pmax(pmin((avg_expression-mean(avg_expression))/s,2.5),-2.5) else 0 }, by=gene]
dot[, pct_expressed_scaled_within_gene := { m<-max(pct_expressed_absolute); if(is.finite(m)&&m>0) 100*pct_expressed_absolute/m else 0 }, by=gene]
setorder(markers, marker_group, gene); dot[, gene:=factor(gene, levels=unique(markers$gene))]; dot[, label:=factor(label, levels=rev(broad_labels))]
p_dot <- ggplot(dot[pct_expressed_scaled_within_gene>0], aes(gene,label)) +
  geom_point(aes(size=pct_expressed_scaled_within_gene, colour=avg_expression_scaled_within_gene), alpha=.9) +
  scale_size_area(max_size=5.5, limits=c(0,100), name="Within-gene detection") +
  scale_colour_gradient2(low="#3B4CC0", mid="#F7F7F7", high="#B40426", midpoint=0, limits=c(-2.5,2.5), oob=scales::squish, name="Within-gene expression") +
  facet_grid(.~marker_group, scales="free_x", space="free_x") + theme_bw(base_size=7) +
  theme(axis.text.x=element_text(angle=55,hjust=1), panel.grid=element_blank(), strip.text=element_text(face="bold")) +
  labs(x="Canonical markers grouped by supported broad class", y="Frozen broad label", subtitle="Point size and color are normalized within each gene; absolute values are retained in the source TSV")
dotplot_png <- file.path(out, "broad_canonical_marker_dotplot_review.png")
ggsave(dotplot_png, p_dot, width=max(12, .18*length(unique(markers$gene))+6), height=max(7,.34*length(broad_labels)+3), dpi=320, bg="white", limitsize=FALSE)
dot_source <- file.path(out, "broad_canonical_marker_dotplot_review_source.tsv")
dot[, `:=`(analysis_view="preconfirmation_final_candidate", evidence_cohort="all_accepted_final_broad")]
fwrite(dot, dot_source, sep="\t")

rel <- function(path) {
  value <- normalizePath(path)
  prefix <- paste0(root, .Platform$file.sep)
  if (startsWith(value, prefix)) substring(value, nchar(prefix) + 1L) else value
}
manifest <- list(
  status="PASS", label_column=broad_col, marker_panel="canonical", broad_labels=broad_labels,
  spatial_png=rel(spatial_png), spatial_png_sha256=file_sha256(spatial_png),
  dotplot_png=rel(dotplot_png), dotplot_png_sha256=file_sha256(dotplot_png),
  dotplot_source=rel(dot_source), dotplot_source_sha256=file_sha256(dot_source),
  palette_tsv=rel(palette_path), palette_tsv_sha256=file_sha256(palette_path),
  scope="preconfirmation_lightweight_only_no_final_deg"
)
write_json(manifest, file.path(out,"review_asset_manifest.json"), auto_unbox=TRUE, pretty=TRUE)
capture.output(sessionInfo(), file=file.path(out,"confirmation_review_sessionInfo.txt"))
message("PASS: lightweight confirmation review assets")
