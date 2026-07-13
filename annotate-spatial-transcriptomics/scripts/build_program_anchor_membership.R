#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(digest)
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

group_definition <- function(groups, features) {
  if (is.null(groups) || !length(groups)) return(list())
  lapply(seq_along(groups), function(i) {
    g <- groups[[i]]
    requested <- unique(as.character(unlist(g$markers)))
    list(
      name = if (is.null(g$name)) paste0("group_", i) else as.character(g$name),
      requested = requested,
      present = intersect(requested, features),
      min_detected = if (is.null(g$min_detected)) 1L else as.integer(g$min_detected)
    )
  })
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required_args <- c("rds", "query-membership", "anchor-config", "out", "cell-id-col")
missing_args <- required_args[!required_args %in% names(args)]
if (length(missing_args)) stop("Missing: ", paste(missing_args, collapse = ", "))

cell_col <- args$`cell-id-col`
cfg <- fromJSON(args$`anchor-config`, simplifyVector = FALSE)
if (!is.null(cfg$base_config)) {
  base_path <- as.character(cfg$base_config)
  base_cfg <- fromJSON(base_path, simplifyVector = FALSE)
  overrides <- cfg$anchor_overrides
  if (!is.null(overrides) && length(overrides)) {
    for (override in overrides) {
      label <- as.character(override$label)
      hit <- which(vapply(base_cfg$anchors, function(spec) identical(as.character(spec$label), label), logical(1)))
      if (length(hit) != 1L) stop("anchor override label must match exactly once: ", label)
      for (field in setdiff(names(override), "label")) base_cfg$anchors[[hit]][[field]] <- override[[field]]
    }
  }
  for (field in setdiff(names(cfg), c("base_config", "anchor_overrides"))) base_cfg[[field]] <- cfg[[field]]
  cfg <- base_cfg
}
obj <- readRDS(args$rds)
assay <- if (is.null(args$assay)) DefaultAssay(obj) else args$assay
layer <- if (is.null(args$layer)) "counts" else args$layer
mat <- tryCatch(
  LayerData(obj[[assay]], layer = layer),
  error = function(e) GetAssayData(obj[[assay]], slot = layer)
)

query <- read_any(args$`query-membership`)
query[[cell_col]] <- as.character(query[[cell_col]])
if (uniqueN(query[[cell_col]]) != nrow(query) || any(!query[[cell_col]] %in% colnames(obj))) {
  stop("invalid query membership")
}
metadata_cols <- c("source_key", "state_tags", "spatial_tags", "qc_tags", "candidate_lineages")
for (nm in metadata_cols) if (!nm %in% names(query)) query[, (nm) := ""]
query[, query_or_anchor := "query"]
set(query, j = "anchor_label", value = rep("", nrow(query)))

specs <- cfg$anchors
if (is.null(specs) || length(specs) < 2L) stop("at least two anchor specifications are required")
selected <- list()
summaries <- list()
query_ids <- query[[cell_col]]
taken <- character()
allow_overlap <- isTRUE(cfg$allow_candidate_overlap_query)

for (i in seq_along(specs)) {
  spec <- specs[[i]]
  candidate <- read_any(spec$candidate_membership)
  candidate[[cell_col]] <- as.character(candidate[[cell_col]])
  ids_all <- intersect(candidate[[cell_col]], colnames(obj))
  ids <- setdiff(ids_all, taken)
  if (!allow_overlap) ids <- setdiff(ids, query_ids)

  groups <- group_definition(spec$required_positive_groups, rownames(mat))
  requested_positive <- unique(c(
    as.character(unlist(spec$positive_markers)),
    unlist(lapply(groups, `[[`, "requested"), use.names = FALSE)
  ))
  positive <- intersect(requested_positive, rownames(mat))
  anti <- intersect(as.character(unlist(spec$anti_markers)), rownames(mat))
  if (!length(positive)) stop("no positive markers for ", spec$label)
  if (any(vapply(groups, function(g) length(g$present) < g$min_detected, logical(1)))) {
    bad <- vapply(groups, function(g) length(g$present) < g$min_detected, logical(1))
    stop(
      "required positive group lacks enough assayed markers for ", spec$label, ": ",
      paste(vapply(groups[bad], `[[`, character(1), "name"), collapse = ",")
    )
  }

  x <- mat[, ids, drop = FALSE]
  n_positive <- Matrix::colSums(x[positive, , drop = FALSE] > 0)
  n_anti <- if (length(anti)) Matrix::colSums(x[anti, , drop = FALSE] > 0) else rep(0, length(ids))
  total_counts <- Matrix::colSums(x)
  positive_score <- Matrix::colSums(x[positive, , drop = FALSE])
  group_counts <- lapply(groups, function(g) Matrix::colSums(x[g$present, , drop = FALSE] > 0))
  group_pass <- if (length(groups)) {
    Reduce(`&`, Map(function(values, g) values >= g$min_detected, group_counts, groups))
  } else {
    rep(TRUE, length(ids))
  }

  min_positive <- if (is.null(spec$min_positive_detected)) 1L else as.integer(spec$min_positive_detected)
  max_anti <- if (is.null(spec$max_anti_detected)) length(anti) else as.integer(spec$max_anti_detected)
  min_total <- if (is.null(spec$min_total_counts)) 0 else as.numeric(spec$min_total_counts)
  max_total <- if (is.null(spec$max_total_counts)) Inf else as.numeric(spec$max_total_counts)
  target_total <- if (is.null(spec$target_total_counts)) NA_real_ else as.numeric(spec$target_total_counts)
  min_selected <- if (is.null(spec$min_selected)) 1L else as.integer(spec$min_selected)

  candidate_metrics <- data.table(
    cell_id = ids,
    n_positive_detected = n_positive,
    n_anti_detected = n_anti,
    total_counts = total_counts,
    positive_score = positive_score,
    required_groups_pass = group_pass
  )
  for (j in seq_along(groups)) {
    candidate_metrics[[paste0("group_", make.names(groups[[j]]$name), "_detected")]] <- group_counts[[j]]
  }
  passing <- candidate_metrics[
    required_groups_pass &
      n_positive_detected >= min_positive &
      n_anti_detected <= max_anti &
      total_counts >= min_total &
      total_counts <= max_total
  ]
  passing[, depth_distance := if (is.finite(target_total)) {
    abs(log1p(total_counts) - log1p(target_total))
  } else 0]
  setorder(passing, depth_distance, -n_positive_detected, n_anti_detected, -positive_score, cell_id)

  n_pass_before_cap <- nrow(passing)
  max_n <- if (is.null(spec$max_n)) n_pass_before_cap else as.integer(spec$max_n)
  if (nrow(passing) > max_n) passing <- passing[seq_len(max_n)]
  if (nrow(passing) < min_selected) {
    stop(
      "insufficient anchors for ", spec$label, ": selected ", nrow(passing),
      ", required ", min_selected
    )
  }

  taken <- c(taken, passing$cell_id)
  passing[, `:=`(
    query_or_anchor = "anchor",
    anchor_label = spec$label,
    source_key = if (is.null(spec$source_key)) paste0("program_anchor|", spec$label) else spec$source_key,
    state_tags = if (is.null(spec$state_tags)) "high_confidence_anchor" else spec$state_tags,
    spatial_tags = "",
    qc_tags = "",
    candidate_lineages = spec$label
  )]
  selected[[i]] <- passing

  group_text <- if (length(groups)) paste(vapply(groups, function(g) {
    paste0(g$name, ">=", g$min_detected, "[", paste(g$present, collapse = ","), "]")
  }, character(1)), collapse = ";") else ""
  summaries[[i]] <- data.table(
    anchor_label = spec$label,
    anchor_tier = if (is.null(spec$anchor_tier)) "strict" else spec$anchor_tier,
    n_candidates_before_overlap = length(ids_all),
    n_candidates_after_overlap = length(ids),
    n_pass_before_cap = n_pass_before_cap,
    n_selected = nrow(passing),
    min_selected = min_selected,
    positive_markers_present = paste(positive, collapse = ","),
    positive_markers_missing = paste(setdiff(requested_positive, positive), collapse = ","),
    anti_markers_present = paste(anti, collapse = ","),
    required_positive_groups = group_text,
    min_positive_detected = min_positive,
    max_anti_detected = max_anti,
    min_total_counts = min_total,
    max_total_counts = max_total,
    target_total_counts = target_total,
    selected_median_total_counts = median(passing$total_counts),
    selected_IQR_total_counts = IQR(passing$total_counts)
  )
}

if (isTRUE(cfg$balance_to_min)) {
  target <- min(vapply(selected, nrow, integer(1)))
  selected <- lapply(selected, function(x) x[seq_len(target)])
}
anchors <- rbindlist(selected, fill = TRUE)
if (anyDuplicated(anchors$cell_id)) stop("anchor specifications overlap; make anchor identities mutually exclusive")
if (allow_overlap) query <- query[!get(cell_col) %in% anchors$cell_id]

query_out <- query[, .(
  cell_id = get(cell_col), query_or_anchor, anchor_label,
  source_key, state_tags, spatial_tags, qc_tags, candidate_lineages
)]
anchor_out <- anchors[, .(
  cell_id, query_or_anchor, anchor_label,
  source_key, state_tags, spatial_tags, qc_tags, candidate_lineages
)]
output <- rbindlist(list(query_out, anchor_out), fill = TRUE)
dir.create(dirname(args$out), recursive = TRUE, showWarnings = FALSE)
fwrite(output, args$out, sep = "\t")
summary_path <- paste0(args$out, ".anchor_summary.tsv")
fwrite(rbindlist(summaries), summary_path, sep = "\t")
manifest <- list(
  status = "PASS",
  n_query = nrow(query_out),
  n_anchors = nrow(anchor_out),
  n_query_promoted_to_anchor = length(intersect(query_ids, anchor_out$cell_id)),
  anchor_counts = as.list(table(anchor_out$anchor_label)),
  balance_to_min = isTRUE(cfg$balance_to_min),
  allow_candidate_overlap_query = allow_overlap,
  membership_sha256 = digest(file = args$out, algo = "sha256"),
  anchor_summary = summary_path
)
write_json(manifest, paste0(args$out, ".manifest.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(manifest, pretty = TRUE, auto_unbox = TRUE), "\n")
