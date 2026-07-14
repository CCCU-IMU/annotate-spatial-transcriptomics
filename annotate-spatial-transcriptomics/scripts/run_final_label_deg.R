#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot resolve script directory")
source(file.path(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))),
                 "seurat_validation_layer.R"))

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

read_any <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) fread(cmd = paste("gzip -dc", shQuote(path))) else fread(path)
}

a <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("rds", "metadata", "out", "cell-id-col", "broad-col", "subtype-col")
missing <- required[!required %in% names(a)]
if (length(missing)) stop("Missing: ", paste(missing, collapse = ", "))

dir.create(file.path(a$out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(a$out, "provenance"), recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(a$rds)
meta <- read_any(a$metadata)
cell_id_col <- a$`cell-id-col`
meta[[cell_id_col]] <- as.character(meta[[cell_id_col]])
if (anyDuplicated(meta[[cell_id_col]])) stop("Metadata contains duplicate observation IDs")
cells <- intersect(colnames(obj), meta[[cell_id_col]])
meta <- meta[match(cells, get(cell_id_col))]
seed <- as.integer(ifelse(is.null(a$seed), 20260713, a$seed))
max_cells <- as.integer(ifelse(is.null(a$`max-cells-per-ident`), 3000, a$`max-cells-per-ident`))
set.seed(seed)

prefix <- ifelse(is.null(a$prefix), "", paste0(a$prefix, "_"))
analysis_view <- ifelse(is.null(a$`analysis-view`),
                        ifelse(nzchar(prefix), sub("_$", "", prefix), "provided_labels"),
                        a$`analysis-view`)
evidence_cohort <- ifelse(is.null(a$`evidence-cohort`), "provided_metadata", a$`evidence-cohort`)

if (inherits(obj, "Seurat")) {
  assay <- ifelse(is.null(a$assay), DefaultAssay(obj), a$assay)
  data_layer <- ifelse(is.null(a$layer), "data", a$layer)
  count_layer <- ifelse(is.null(a$`count-layer`), "counts", a$`count-layer`)
  if (!assay %in% Assays(obj)) stop("Missing assay: ", assay)
  assert_seurat_validation_layer(
    obj, assay = assay, data_layer = data_layer, count_layer = count_layer,
    manifest_path = a$`validation-manifest`, object_path = a$rds
  )
}

if (!is.null(a$`include-col`)) {
  if (is.null(a$`include-values`)) stop("--include-values is required with --include-col")
  allowed <- strsplit(a$`include-values`, ",", fixed = TRUE)[[1]]
  meta <- meta[tolower(as.character(get(a$`include-col`))) %in% tolower(allowed)]
  cells <- meta[[cell_id_col]]
}

write_level <- function(level, column) {
  level_view <- if (!is.null(a[[paste0(level, "-analysis-view")]])) a[[paste0(level, "-analysis-view")]] else analysis_view
  level_cohort <- if (!is.null(a[[paste0(level, "-evidence-cohort")]])) a[[paste0(level, "-evidence-cohort")]] else evidence_cohort
  labels <- as.character(meta[[column]])
  valid <- !is.na(labels) & nzchar(labels)
  labs <- sort(unique(labels[valid]))
  if (length(labs) < 2) stop(level, " has fewer than two labels")

  if (inherits(obj, "Seurat")) {
    query <- subset(obj, cells = cells[valid])
    DefaultAssay(query) <- assay
    label_map <- setNames(labels[valid], cells[valid])
    Idents(query) <- factor(label_map[colnames(query)], levels = labs)
    deg <- as.data.table(FindAllMarkers(
      query, assay = assay, slot = data_layer, only.pos = TRUE, test.use = "wilcox",
      min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1,
      max.cells.per.ident = max_cells, random.seed = seed, densify = FALSE, verbose = TRUE
    ))
    setnames(deg, "cluster", "label")
  } else if (inherits(obj, "SingleCellExperiment") || inherits(obj, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("SummarizedExperiment required")
    available <- SummarizedExperiment::assayNames(obj)
    data_name <- if (!is.null(a$`data-assay`)) a$`data-assay` else if ("normcounts" %in% available) "normcounts" else if ("logcounts" %in% available) "logcounts" else stop("Specify --data-assay")
    count_name <- if (!is.null(a$`count-assay`)) a$`count-assay` else "counts"
    normalized <- SummarizedExperiment::assay(obj, data_name)[, cells[valid], drop = FALSE]
    counts <- SummarizedExperiment::assay(obj, count_name)[, cells[valid], drop = FALSE]
    groups <- factor(labels[valid], levels = labs)
    design <- sparse.model.matrix(~0 + groups); colnames(design) <- labs
    n <- colSums(design)
    sums <- as.matrix(normalized %*% design)
    detected <- as.matrix((counts > 0) %*% design)
    average <- sweep(sums, 2, n, "/")
    pct <- 100 * sweep(detected, 2, n, "/")
    total <- Matrix::rowSums(normalized)
    deg <- rbindlist(lapply(seq_along(labs), function(j) {
      rest <- (total - sums[, j]) / (sum(valid) - n[j])
      data.table(gene = rownames(obj), label = labs[j], avg_expression = average[, j],
                 pct_expressed_absolute = pct[, j], avg_expression_rest = rest,
                 avg_log2FC = log2((average[, j] + 1e-4) / (rest + 1e-4)),
                 n_observations = n[j])
    }))
    setorder(deg, label, -avg_log2FC)
  } else {
    stop("Unsupported object class")
  }

  deg[, `:=`(analysis_view = level_view, evidence_cohort = level_cohort)]
  fwrite(deg, file.path(a$out, "tables", paste0(prefix, level, "_DEG_one_vs_rest_all.tsv")), sep = "\t")
  lfc_column <- intersect(c("avg_log2FC", "avg_logFC"), names(deg))[1]
  if ("p_val_adj" %in% names(deg)) {
    top <- deg[p_val_adj < 0.05 & get(lfc_column) > 0][order(label, p_val_adj, -get(lfc_column)), head(.SD, 100), by = label]
  } else {
    top <- deg[order(label, -get(lfc_column)), head(.SD, 100), by = label]
  }
  fwrite(top, file.path(a$out, "tables", paste0(prefix, level, "_DEG_top100.tsv")), sep = "\t")
}

write_level("broad", a$`broad-col`)
write_level("subtype", a$`subtype-col`)
capture.output(sessionInfo(), file = file.path(a$out, "provenance", "DEG_sessionInfo.txt"))
