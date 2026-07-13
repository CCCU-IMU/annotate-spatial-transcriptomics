#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
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

read_any <- function(path) {
  if (grepl("\\.gz$", path)) fread(cmd = paste("gzip -dc", shQuote(path))) else fread(path)
}

resolve_config <- function(path) {
  cfg <- fromJSON(path, simplifyVector = FALSE)
  if (is.null(cfg$base_config)) return(cfg)
  base <- resolve_config(as.character(cfg$base_config))
  overrides <- cfg$anchor_overrides
  if (!is.null(overrides) && length(overrides)) {
    for (override in overrides) {
      label <- as.character(override$label)
      hit <- which(vapply(base$anchors, function(spec) identical(as.character(spec$label), label), logical(1)))
      if (length(hit) != 1L) stop("anchor override label must match exactly once: ", label)
      for (field in setdiff(names(override), "label")) base$anchors[[hit]][[field]] <- override[[field]]
    }
  }
  additions <- cfg$additional_programs
  if (!is.null(additions) && length(additions)) {
    if (is.null(base$programs)) base$programs <- base$anchors
    base$programs <- c(base$programs, additions)
  }
  for (field in setdiff(names(cfg), c("base_config", "anchor_overrides", "additional_programs"))) base[[field]] <- cfg[[field]]
  base
}

group_definition <- function(groups, features) {
  if (is.null(groups) || !length(groups)) return(list())
  lapply(seq_along(groups), function(i) {
    group <- groups[[i]]
    requested <- unique(as.character(unlist(group$markers)))
    list(
      name = if (is.null(group$name)) paste0("group_", i) else as.character(group$name),
      present = intersect(requested, features),
      min_detected = if (is.null(group$min_detected)) 1L else as.integer(group$min_detected)
    )
  })
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("rds", "cluster-glob", "program-config", "out", "cell-id-col")
missing <- required[!required %in% names(args)]
if (length(missing)) stop("Missing: ", paste(missing, collapse = ", "))

cluster_paths <- sort(Sys.glob(args$`cluster-glob`))
if (!length(cluster_paths)) stop("cluster-glob matched no files")
cell_col <- args$`cell-id-col`
cluster_col <- if (is.null(args$`cluster-col`)) "cluster" else args$`cluster-col`
cfg <- resolve_config(args$`program-config`)
specs <- if (!is.null(cfg$programs)) cfg$programs else cfg$anchors
if (is.null(specs) || length(specs) < 2L) stop("program config must contain at least two programs/anchors")

obj <- readRDS(args$rds)
assay <- if (is.null(args$assay)) DefaultAssay(obj) else args$assay
layer <- if (is.null(args$layer)) "counts" else args$layer
mat <- tryCatch(
  LayerData(obj[[assay]], layer = layer),
  error = function(e) GetAssayData(obj[[assay]], slot = layer)
)

clusters <- lapply(cluster_paths, function(path) {
  frame <- read_any(path)
  if (!all(c(cell_col, cluster_col) %in% names(frame))) stop("cluster file lacks cell/cluster columns: ", path)
  frame[[cell_col]] <- as.character(frame[[cell_col]])
  frame[[cluster_col]] <- as.character(frame[[cluster_col]])
  if (uniqueN(frame[[cell_col]]) != nrow(frame)) stop("duplicate cell IDs in ", path)
  if (any(!frame[[cell_col]] %in% colnames(obj))) stop("cluster cell IDs missing from object: ", path)
  resolution <- if ("resolution" %in% names(frame)) unique(as.character(frame$resolution)) else {
    stem <- basename(path)
    hit <- regmatches(stem, regexpr("res[0-9]+p?[0-9]*", stem, ignore.case = TRUE))
    sub("^res", "", gsub("p", ".", hit, ignore.case = TRUE), ignore.case = TRUE)
  }
  if (length(resolution) != 1L || !nzchar(resolution)) stop("cannot derive one resolution from ", path)
  setnames(frame, c(cell_col, cluster_col), c("cell_id", "cluster"))
  frame[, `:=`(resolution = resolution, cluster_file = normalizePath(path))]
  frame[, .(cell_id, cluster, resolution, cluster_file)]
})

reference_ids <- clusters[[1]]$cell_id
if (any(vapply(clusters, function(frame) !setequal(frame$cell_id, reference_ids), logical(1)))) {
  stop("all cluster files must contain the same observation membership")
}
ids <- intersect(reference_ids, colnames(obj))
x <- mat[, ids, drop = FALSE]

program_metrics <- list()
program_manifest <- list()
for (i in seq_along(specs)) {
  spec <- specs[[i]]
  label <- as.character(spec$label)
  groups <- group_definition(spec$required_positive_groups, rownames(mat))
  requested_positive <- unique(c(
    as.character(unlist(spec$positive_markers)),
    unlist(lapply(groups, `[[`, "present"), use.names = FALSE)
  ))
  positive <- intersect(requested_positive, rownames(mat))
  anti <- intersect(as.character(unlist(spec$anti_markers)), rownames(mat))
  if (!length(positive)) stop("no assayed positive markers for ", label)
  hits <- as.integer(Matrix::colSums(x[positive, , drop = FALSE] > 0))
  anti_hits <- if (length(anti)) as.integer(Matrix::colSums(x[anti, , drop = FALSE] > 0)) else rep(0L, length(ids))
  groups_pass <- rep(TRUE, length(ids))
  group_text <- character()
  for (group in groups) {
    if (length(group$present) < group$min_detected) stop("required group lacks enough assayed genes for ", label, ": ", group$name)
    values <- as.integer(Matrix::colSums(x[group$present, , drop = FALSE] > 0))
    groups_pass <- groups_pass & values >= group$min_detected
    group_text <- c(group_text, paste0(group$name, ">=", group$min_detected, "[", paste(group$present, collapse = ","), "]"))
  }
  min_positive <- if (is.null(spec$min_positive_detected)) 1L else as.integer(spec$min_positive_detected)
  max_anti <- if (is.null(spec$max_anti_detected)) length(anti) else as.integer(spec$max_anti_detected)
  strong <- groups_pass & hits >= min_positive & anti_hits <= max_anti
  program_metrics[[i]] <- data.table(
    cell_id = ids,
    program = label,
    positive_hits = hits,
    positive_fraction = hits / length(positive),
    anti_hits = anti_hits,
    required_groups_pass = groups_pass,
    strong_program = strong
  )
  program_manifest[[i]] <- data.table(
    program = label,
    n_positive_markers = length(positive),
    positive_markers = paste(positive, collapse = ","),
    n_anti_markers = length(anti),
    anti_markers = paste(anti, collapse = ","),
    required_groups = paste(group_text, collapse = ";"),
    min_positive_detected = min_positive,
    max_anti_detected = max_anti
  )
}
metrics <- rbindlist(program_metrics)

cell_wide <- dcast(metrics, cell_id ~ program, value.var = "positive_fraction")
program_cols <- setdiff(names(cell_wide), "cell_id")
score_matrix <- as.matrix(cell_wide[, ..program_cols])
top_index <- max.col(score_matrix, ties.method = "first")
top_score <- score_matrix[cbind(seq_len(nrow(score_matrix)), top_index)]
second_score <- apply(score_matrix, 1, function(values) sort(values, decreasing = TRUE)[min(2L, length(values))])
cell_top <- data.table(
  cell_id = cell_wide$cell_id,
  top_program = program_cols[top_index],
  top_program_fraction = top_score,
  program_fraction_margin = top_score - second_score
)

summary_outputs <- list()
top_outputs <- list()
for (j in seq_along(clusters)) {
  cluster_frame <- clusters[[j]]
  joined <- merge(cluster_frame, metrics, by = "cell_id", allow.cartesian = TRUE, sort = FALSE)
  summary_outputs[[j]] <- joined[, .(
    n_observations = .N,
    mean_positive_hits = mean(positive_hits),
    median_positive_hits = as.numeric(median(positive_hits)),
    any_positive_fraction = mean(positive_hits >= 1L),
    required_groups_pass_fraction = mean(required_groups_pass),
    strong_program_fraction = mean(strong_program),
    mean_anti_hits = mean(anti_hits),
    anti_any_fraction = mean(anti_hits >= 1L),
    anti_dominant_fraction = mean(anti_hits > positive_hits)
  ), by = .(resolution, cluster, cluster_file, program)]
  top_joined <- merge(cluster_frame, cell_top, by = "cell_id", sort = FALSE)
  top_outputs[[j]] <- top_joined[, .(
    n = .N,
    mean_top_program_fraction = mean(top_program_fraction),
    mean_program_fraction_margin = mean(program_fraction_margin)
  ), by = .(resolution, cluster, cluster_file, top_program)][
    , top_program_cluster_fraction := n / sum(n), by = .(resolution, cluster)
  ]
}

out_dir <- args$out
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(rbindlist(summary_outputs), file.path(out_dir, "cluster_program_summary.tsv"), sep = "\t")
fwrite(rbindlist(top_outputs), file.path(out_dir, "cluster_top_program_composition.tsv"), sep = "\t")
fwrite(cell_top, file.path(out_dir, "cell_top_program.tsv.gz"), sep = "\t")
cell_program_metrics_path <- NULL
if (!is.null(args$`write-cell-program-metrics`)) {
  cell_program_metrics_path <- file.path(out_dir, "cell_program_metrics.tsv.gz")
  fwrite(metrics, cell_program_metrics_path, sep = "\t")
}
fwrite(rbindlist(program_manifest), file.path(out_dir, "program_manifest.tsv"), sep = "\t")
write_json(list(
  status = "PASS",
  n_observations = length(ids),
  n_resolutions = length(clusters),
  n_programs = length(specs),
  cluster_program_summary = file.path(out_dir, "cluster_program_summary.tsv"),
  cluster_top_program_composition = file.path(out_dir, "cluster_top_program_composition.tsv"),
  cell_program_metrics = cell_program_metrics_path,
  warning = "Program fractions are evidence summaries, not automatic labels; review DEG, anchor distance, source/QC composition and spatial morphology before writeback."
), file.path(out_dir, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)
cat("PASS\n")
