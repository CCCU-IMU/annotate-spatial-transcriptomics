#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(RANN)
})

parse_args <- function(values) {
  result <- list()
  i <- 1L
  while (i <= length(values)) {
    key <- values[[i]]
    if (!startsWith(key, "--") || i == length(values)) {
      stop("arguments must use --name value pairs")
    }
    result[[substring(key, 3L)]] <- values[[i + 1L]]
    i <- i + 2L
  }
  result
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("rds", "profile", "catalog", "partitions", "out", "observation-unit")
missing <- required[!required %in% names(args)]
if (length(missing)) stop("missing required arguments: ", paste(missing, collapse = ", "))
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("cannot resolve scorer script path")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
sha256 <- function(path) {
  output <- system2("sha256sum", shQuote(path), stdout = TRUE)
  strsplit(output[[1L]], "\\s+")[[1L]][[1L]]
}

input_rds <- normalizePath(args$rds, mustWork = TRUE)
profile_path <- normalizePath(args$profile, mustWork = TRUE)
catalog_path <- normalizePath(args$catalog, mustWork = TRUE)
threshold_registry_path <- normalizePath(
  if ("threshold-registry" %in% names(args)) {
    args[["threshold-registry"]]
  } else {
    file.path(dirname(script_path), "..", "references", "controller_thresholds_v2_2.json")
  },
  mustWork = TRUE
)
threshold_registry <- fromJSON(
  threshold_registry_path, simplifyVector = FALSE
)
if (!identical(as.character(threshold_registry$schema_version), "2.2")) {
  stop("controller threshold registry is not schema 2.2")
}
scoring_policy <- threshold_registry$scoring_policy
partition_path <- normalizePath(args$partitions, mustWork = TRUE)
out_dir <- args$out
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

seed <- if ("seed" %in% names(args)) as.integer(args$seed) else 2200L
knn_k <- if ("knn-k" %in% names(args)) as.integer(args[["knn-k"]]) else 31L
workers <- if ("workers" %in% names(args)) as.integer(args$workers) else 1L
if (!is.finite(workers) || workers < 1L) stop("--workers must be >= 1")
if (.Platform$OS.type != "unix") workers <- 1L
observation_unit <- tolower(trimws(args[["observation-unit"]]))
if (!observation_unit %in% c("cell", "nucleus", "cellbin", "spot")) {
  stop("--observation-unit must be cell, nucleus, cellbin or spot")
}
direct_weight <- as.numeric(scoring_policy$direct_weight)
local_weight <- as.numeric(scoring_policy$local_weight)
anti_weight <- as.numeric(scoring_policy$anti_weight)
family_active_threshold <- as.numeric(scoring_policy$family_active_threshold)
local_gene_detection_fraction <- as.numeric(
  scoring_policy$local_gene_detection_fraction
)
if (
  any(!is.finite(c(
    direct_weight, local_weight, anti_weight,
    family_active_threshold, local_gene_detection_fraction
  ))) || abs(direct_weight + local_weight - 1) > 1e-12
) stop("invalid scoring policy in controller threshold registry")
grid_evidence_only <- "grid-evidence-only" %in% names(args) &&
  tolower(args[["grid-evidence-only"]]) %in% c("1", "true", "yes")
set.seed(seed)

read_tsv <- function(path) {
  if (grepl("\\.gz$", path)) {
    fread(cmd = paste("gzip -dc", shQuote(path)), colClasses = "character")
  } else {
    fread(path, colClasses = "character")
  }
}
write_tsv_atomic <- function(value, path) {
  compressed <- grepl("\\.gz$", path)
  temporary <- paste0(
    path, ".tmp-", Sys.getpid(),
    if (compressed) ".gz" else ""
  )
  if (file.exists(temporary)) unlink(temporary)
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  if (compressed && nrow(value) == 0L) {
    # data.table::fwrite() leaves a header-only .gz stream without the gzip
    # end marker.  Write and close that valid empty table explicitly.
    connection <- gzfile(temporary, open = "wt")
    tryCatch(
      writeLines(paste(names(value), collapse = "\t"), connection),
      finally = close(connection)
    )
  } else {
    fwrite(value, temporary, sep = "\t")
  }
  if (compressed) {
    status <- system2("gzip", c("-t", shQuote(temporary)))
    if (!identical(status, 0L)) stop("invalid gzip output: ", temporary)
  }
  if (!file.rename(temporary, path)) {
    stop("could not atomically publish output: ", path)
  }
  invisible(path)
}
get_path <- function(document, path) {
  current <- document
  for (key in strsplit(path, "\\.", fixed = FALSE)[[1L]]) {
    current <- current[[key]]
    if (is.null(current)) return(NULL)
  }
  current
}
as_chars <- function(value) unique(toupper(as.character(unlist(value))))
safe_name <- function(value) gsub("[^A-Za-z0-9]+", "_", tolower(value))
or_value <- function(value, fallback) {
  if (is.null(value) || !length(value) || is.na(value[[1L]])) fallback else value[[1L]]
}
top_two_mean <- function(matrix_value) {
  if (!ncol(matrix_value)) return(rep(0, nrow(matrix_value)))
  if (ncol(matrix_value) == 1L) return(matrix_value[, 1L])
  first <- matrix_value[, 1L]
  second <- rep(-Inf, nrow(matrix_value))
  for (column_index in 2:ncol(matrix_value)) {
    value <- matrix_value[, column_index]
    replace_first <- value > first
    second[replace_first] <- first[replace_first]
    second[!replace_first] <- pmax(
      second[!replace_first], value[!replace_first]
    )
    first[replace_first] <- value[replace_first]
  }
  (first + second) / 2
}
catalog_candidates <- function(catalog) {
  rows <- catalog$candidate_boundaries
  represented <- vapply(rows, function(row) {
    as.character(or_value(row$release_fine_label, ""))
  }, character(1))
  for (items in catalog$machine_actionable_fine_candidate_catalog) {
    for (item in items) {
      fine_label <- as.character(or_value(item$release_label, ""))
      if (!nzchar(fine_label) || fine_label %in% represented) next
      parent <- as.character(or_value(item$parent_release_label, ""))
      context <- unlist(Filter(nzchar, c(
        as.character(or_value(item$context_gate, "")),
        as.character(or_value(item$required_discriminator, ""))
      )))
      item$candidate_role <- "fine"
      item$release_broad_label <- parent
      item$release_fine_label <- fine_label
      item$parent_broad_label <- parent
      item$writeback_strategy <- as.character(or_value(
        item$writeback_strategy, "supported_subset_with_parent_lock"
      ))
      item$specificity_priority <- as.numeric(or_value(item$specificity_priority, 70))
      item$hard_anti_families <- or_value(item$hard_anti_families, list())
      item$soft_anti_families <- or_value(item$soft_anti_families, list())
      item$context_requirements <- or_value(item$context_requirements, context)
      item$review_required <- TRUE
      rows[[length(rows) + 1L]] <- item
      represented <- c(represented, fine_label)
    }
  }
  rows
}

profile <- fromJSON(profile_path, simplifyVector = FALSE)
catalog <- fromJSON(catalog_path, simplifyVector = FALSE)
candidates <- catalog_candidates(catalog)
candidate_ids <- vapply(candidates, `[[`, character(1), "candidate_id")
if (anyDuplicated(candidate_ids) || !length(candidate_ids)) stop("invalid candidate catalog")
catalog_by_id <- setNames(candidates, candidate_ids)

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("--rds must contain a Seurat object")
analysis_membership_path <- NULL
if ("analysis-membership" %in% names(args)) {
  analysis_membership_path <- normalizePath(
    args[["analysis-membership"]], mustWork = TRUE
  )
  analysis_membership <- read_tsv(analysis_membership_path)
  if (!"cell_id" %in% names(analysis_membership)) {
    stop("analysis membership must contain cell_id")
  }
  analysis_ids <- as.character(analysis_membership$cell_id)
  if (
    !length(analysis_ids) || any(!nzchar(analysis_ids)) ||
      anyDuplicated(analysis_ids)
  ) {
    stop("analysis membership contains empty or duplicate observations")
  }
  missing_analysis <- setdiff(analysis_ids, colnames(obj))
  if (length(missing_analysis)) {
    stop("analysis membership contains observations absent from the object")
  }
  obj <- subset(obj, cells = analysis_ids)
}
assay_name <- if ("assay" %in% names(args)) args$assay else {
  assay_candidates <- Assays(obj)
  data_backed <- assay_candidates[vapply(assay_candidates, function(value) {
    value != "SCT" && "data" %in% Layers(obj[[value]])
  }, logical(1))]
  preferred <- c("RNA", "Spatial")
  preferred <- preferred[preferred %in% data_backed]
  if (length(preferred)) preferred[[1L]] else if (length(data_backed)) {
    data_backed[[1L]]
  } else {
    stop("scoring requires a project-local non-SCT full-feature data layer")
  }
}
layer_name <- if ("layer" %in% names(args)) args$layer else "data"
if (!assay_name %in% Assays(obj)) stop("requested assay is absent")
layers <- Layers(obj[[assay_name]])
if (!layer_name %in% layers) {
  stop("full-feature normalized layer is absent; prepare a project-local validation object first")
}
expr <- LayerData(obj[[assay_name]], layer = layer_name)
if (!inherits(expr, "sparseMatrix")) expr <- as(expr, "dgCMatrix")
if (nrow(expr) < 1000L) stop("scoring requires a full-feature expression layer")

partitions <- read_tsv(partition_path)
needed_partition <- c("cell_id", "resolution", "cluster", "resolution_role")
if (!all(needed_partition %in% names(partitions))) {
  stop("partition table must contain: ", paste(needed_partition, collapse = ", "))
}
partitions[, cell_id := as.character(cell_id)]
if (!"boundary_id" %in% names(partitions)) partitions[, boundary_id := "whole_tissue"]
partition_resolutions <- unique(partitions$resolution)
if (length(partition_resolutions) < 3L) {
  stop("selected-and-neighbors scoring requires at least three partition resolutions")
}
if (!all(partitions$cell_id %chin% colnames(obj))) stop("partition contains unknown observations")
partition_roles <- c("selected", "neighbor_1", "neighbor_2")
if (grid_evidence_only) {
  numeric_resolutions <- suppressWarnings(as.numeric(partition_resolutions))
  ordered_resolutions <- if (all(is.finite(numeric_resolutions))) {
    partition_resolutions[order(numeric_resolutions)]
  } else {
    sort(partition_resolutions)
  }
  partitions_by_resolution <- lapply(ordered_resolutions, function(resolution_value) {
    value <- partitions[resolution == resolution_value]
    if (
      nrow(value) != uniqueN(value$cell_id) ||
        nrow(value) != ncol(obj)
    ) {
      stop("grid partition must contain exactly one row per expression observation")
    }
    value <- value[match(colnames(obj), cell_id)]
    if (anyNA(value$cell_id)) {
      stop("grid partition resolution does not cover the expression object")
    }
    value
  })
  # Grid scoring evaluates every supplied resolution. These three views are
  # placeholders for shared downstream code only; their role names are not
  # written back and therefore cannot pre-empt resolution selection.
  partition_by_role <- setNames(partitions_by_resolution[seq_len(3L)], partition_roles)
  selected <- partition_by_role$selected
} else {
  selected <- partitions[resolution_role == "selected"]
  if (nrow(selected) != uniqueN(selected$cell_id)) {
    stop("every observation must have exactly one selected-resolution partition")
  }
  selected <- selected[match(colnames(obj), cell_id)]
  if (anyNA(selected$cell_id)) stop("selected partition does not cover the expression object")
  partition_by_role <- setNames(lapply(partition_roles, function(role) {
    value <- partitions[resolution_role == role]
    if (nrow(value) != uniqueN(value$cell_id)) {
      stop(role, " partition must contain exactly one row per observation")
    }
    value <- value[match(colnames(obj), cell_id)]
    if (anyNA(value$cell_id)) stop(role, " partition does not cover the expression object")
    value
  }), partition_roles)
}

metadata <- as.data.table(obj[[]], keep.rownames = "cell_id")
metadata <- metadata[match(colnames(obj), cell_id)]
x_col <- if ("x-col" %in% names(args)) args[["x-col"]] else "x"
y_col <- if ("y-col" %in% names(args)) args[["y-col"]] else "y"
if (!all(c(x_col, y_col) %in% names(metadata))) {
  stop("metadata lacks requested spatial coordinate columns")
}
coordinates <- as.matrix(metadata[, .(
  x = as.numeric(get(x_col)),
  y = as.numeric(get(y_col))
)])
if (any(!is.finite(coordinates))) stop("spatial coordinates contain non-finite values")

features_upper <- toupper(rownames(expr))
feature_lookup <- setNames(seq_along(features_upper), features_upper)
n_obs <- ncol(expr)
nn_k <- min(max(3L, knn_k), n_obs)
neighbors <- RANN::nn2(coordinates, coordinates, k = nn_k)$nn.idx

positive_families <- list()
soft_anti_genes <- list()
hard_anti_sources <- list()
candidate_meta <- list()
for (candidate in candidates) {
  candidate_id <- candidate$candidate_id
  rule <- get_path(profile, candidate$profile_program)
  if (is.null(rule)) stop("profile program is absent: ", candidate$profile_program)
  families <- rule$positive_families
  inherited_rule <- NULL
  if (is.null(families) && identical(candidate$candidate_role, "state")) {
    families <- candidate$state_positive_families
  }
  if (is.null(families) && identical(candidate$candidate_role, "fine")) {
    parent <- as.character(candidate$parent_broad_label)
    parent_match <- Filter(function(item) {
      identical(as.character(item$candidate_role), "broad") &&
        identical(as.character(item$release_broad_label), parent)
    }, candidates)
    if (!length(parent_match)) {
      stop(candidate_id, " lacks a resolvable broad parent")
    }
    inherited_rule <- get_path(profile, parent_match[[1L]]$profile_program)
    subtype_genes <- as_chars(rule)
    parent_genes <- setdiff(
      as_chars(unlist(inherited_rule$positive_families)), subtype_genes
    )
    configured_fine_families <- candidate$fine_positive_families
    if (!is.null(configured_fine_families)) {
      configured_genes <- as_chars(unlist(configured_fine_families))
      parent_genes <- setdiff(
        as_chars(unlist(inherited_rule$positive_families)), configured_genes
      )
      families <- c(
        list(parent_identity = parent_genes),
        configured_fine_families
      )
    } else {
      families <- list(
        parent_identity = parent_genes,
        fine_discriminator = subtype_genes
      )
    }
    if (!length(candidate$hard_anti_families)) {
      candidate$hard_anti_families <- parent_match[[1L]]$hard_anti_families
    }
  }
  if (is.null(families) || length(families) < 2L) {
    stop(candidate_id, " lacks two positive marker families")
  }
  for (family_name in names(families)) {
    positive_families[[paste(candidate_id, family_name, sep = "::")]] <-
      as_chars(families[[family_name]])
  }
  rule_anti <- if (is.list(rule)) rule$anti_programs else NULL
  inherited_anti <- if (is.list(inherited_rule)) inherited_rule$anti_programs else NULL
  soft_anti_genes[[candidate_id]] <- as_chars(
    if (is.null(rule_anti)) inherited_anti else rule_anti
  )
  if (candidate_id == "oocyte" && is.list(rule) && !is.null(rule$contradictory_somatic)) {
    soft_anti_genes[[candidate_id]] <- as_chars(rule$contradictory_somatic)
  }
  hard_anti <- candidate$hard_anti_families
  unit_specific <- candidate$hard_anti_families_by_observation_unit
  if (is.list(unit_specific) && !is.null(unit_specific[[observation_unit]])) {
    hard_anti <- unit_specific[[observation_unit]]
  }
  hard_anti_sources[[candidate_id]] <- as.character(unlist(hard_anti))
  candidate_meta[[candidate_id]] <- candidate
}

all_marker_genes <- unique(c(unlist(positive_families), unlist(soft_anti_genes)))
available <- intersect(all_marker_genes, names(feature_lookup))
fwrite(
  data.table(gene = all_marker_genes, available = all_marker_genes %chin% available),
  file.path(out_dir, "tables", "marker_availability.tsv"),
  sep = "\t"
)

gene_cache <- new.env(parent = emptyenv())
compute_gene_score <- function(gene) {
  index <- feature_lookup[[gene]]
  if (is.null(index) || is.na(index)) {
    answer <- list(
      direct = numeric(n_obs), local = numeric(n_obs),
      detected = logical(n_obs)
    )
  } else {
    value <- as.numeric(expr[index, ])
    nonzero <- value[value > 0]
    scale <- if (length(nonzero)) {
      as.numeric(quantile(nonzero, 0.95, names = FALSE, type = 8))
    } else {
      1
    }
    if (!is.finite(scale) || scale <= 0) scale <- 1
    direct <- pmin(value / scale, 1)
    local <- rowMeans(matrix(
      direct[neighbors], nrow = nrow(neighbors), ncol = ncol(neighbors)
    ))
    answer <- list(direct = direct, local = local, detected = value > 0)
  }
  answer
}
available_scores <- parallel::mclapply(
  available, compute_gene_score,
  mc.cores = min(workers, max(1L, length(available))),
  mc.preschedule = TRUE
)
for (index in seq_along(available)) {
  assign(available[[index]], available_scores[[index]], envir = gene_cache)
}
rm(available_scores)
gene_score <- function(gene) {
  if (exists(gene, envir = gene_cache, inherits = FALSE)) {
    return(get(gene, envir = gene_cache))
  }
  compute_gene_score(gene)
}

family_ids <- names(positive_families)
family_direct <- matrix(
  0, nrow = n_obs, ncol = length(family_ids),
  dimnames = list(colnames(obj), family_ids)
)
family_local <- family_direct
family_detected <- matrix(
  0L, nrow = n_obs, ncol = length(family_ids),
  dimnames = list(colnames(obj), family_ids)
)
family_local_detected <- family_detected
family_availability <- list()

for (family_index in seq_along(family_ids)) {
  family_id <- family_ids[[family_index]]
  genes <- intersect(positive_families[[family_id]], available)
  family_availability[[family_index]] <- data.table(
    family_id = family_id,
    available_marker_count = length(genes),
    total_marker_count = length(positive_families[[family_id]]),
    available_genes = paste(genes, collapse = ",")
  )
  if (!length(genes)) next
  direct_values <- matrix(0, nrow = n_obs, ncol = length(genes))
  local_values <- direct_values
  direct_detected <- matrix(FALSE, nrow = n_obs, ncol = length(genes))
  local_detected <- direct_detected
  for (gene_index in seq_along(genes)) {
    values <- gene_score(genes[[gene_index]])
    direct_values[, gene_index] <- values$direct
    local_values[, gene_index] <- values$local
    direct_detected[, gene_index] <- values$detected
    local_detected[, gene_index] <- rowMeans(matrix(
      as.numeric(values$detected)[neighbors],
      nrow = nrow(neighbors), ncol = ncol(neighbors)
    )) >= local_gene_detection_fraction
  }
  family_direct[, family_index] <- top_two_mean(direct_values)
  family_local[, family_index] <- top_two_mean(local_values)
  family_detected[, family_index] <- rowSums(direct_detected)
  family_local_detected[, family_index] <- rowSums(local_detected)
}
fwrite(
  rbindlist(family_availability),
  file.path(out_dir, "tables", "positive_family_availability.tsv"),
  sep = "\t"
)

n_candidates <- length(candidate_ids)
candidate_score <- matrix(
  0, nrow = n_obs, ncol = n_candidates,
  dimnames = list(colnames(obj), candidate_ids)
)
candidate_direct <- candidate_score
candidate_local <- candidate_score
candidate_positive_family_count <- matrix(
  0L, nrow = n_obs, ncol = n_candidates,
  dimnames = list(colnames(obj), candidate_ids)
)
candidate_positive_gene_count <- candidate_positive_family_count
candidate_family_coherent <- matrix(
  FALSE, nrow = n_obs, ncol = n_candidates,
  dimnames = list(colnames(obj), candidate_ids)
)
candidate_identity_core_coherent <- candidate_family_coherent
candidate_identity_core_direct <- candidate_family_coherent
candidate_release_family_coherent <- matrix(
  TRUE, nrow = n_obs, ncol = n_candidates,
  dimnames = list(colnames(obj), candidate_ids)
)
candidate_soft_anti <- candidate_score
candidate_direct_anti_count <- candidate_positive_family_count
candidate_direct_anti_family_count <- candidate_positive_family_count
candidate_hard_contradiction <- candidate_family_coherent
family_active <- matrix(
  FALSE, nrow = n_obs, ncol = length(family_ids),
  dimnames = list(colnames(obj), family_ids)
)

for (candidate_index in seq_along(candidate_ids)) {
  candidate_id <- candidate_ids[[candidate_index]]
  family_columns <- grep(
    paste0("^", candidate_id, "::"),
    colnames(family_direct), value = TRUE
  )
  direct_matrix <- family_direct[, family_columns, drop = FALSE]
  local_matrix <- family_local[, family_columns, drop = FALSE]
  combined_matrix <- direct_weight * direct_matrix + local_weight * local_matrix
  coherent <- (
    family_detected[, family_columns, drop = FALSE] >= 2L |
      (
        family_detected[, family_columns, drop = FALSE] >= 1L &
          family_local_detected[, family_columns, drop = FALSE] >= 2L
      ) |
      family_local_detected[, family_columns, drop = FALSE] >= 3L
  ) & combined_matrix >= family_active_threshold
  family_active[, family_columns] <- coherent
  candidate_positive_family_count[, candidate_index] <- rowSums(coherent)
  candidate_positive_gene_count[, candidate_index] <- rowSums(
    family_detected[, family_columns, drop = FALSE]
  )
  candidate_family_coherent[, candidate_index] <- rowSums(coherent) >= 1L
  seed_families <- as.character(unlist(
    catalog_by_id[[candidate_id]]$seed_required_positive_families
  ))
  if (length(seed_families)) {
    seed_columns <- paste(candidate_id, seed_families, sep = "::")
    missing_seed <- setdiff(seed_columns, family_columns)
    if (!length(missing_seed)) {
      candidate_identity_core_coherent[, candidate_index] <- rowSums(
        coherent[, seed_columns, drop = FALSE]
      ) >= 1L
      candidate_identity_core_direct[, candidate_index] <- rowSums(
        coherent[, seed_columns, drop = FALSE] &
          family_detected[, seed_columns, drop = FALSE] >= 1L
      ) >= 1L
    }
  } else {
    # Catalog entries without an explicit identity family remain visible, but
    # their core is still tied to a directly detected coherent family. Local
    # smoothing alone must never create an identity-grade seed.
    directly_supported_core <- rowSums(
      coherent & family_detected[, family_columns, drop = FALSE] >= 1L,
      na.rm = TRUE
    ) >= 1L
    candidate_identity_core_coherent[, candidate_index] <-
      directly_supported_core
    candidate_identity_core_direct[, candidate_index] <-
      directly_supported_core
  }
  required_families <- catalog_by_id[[candidate_id]]$required_positive_families
  if (length(required_families)) {
    required_columns <- paste(candidate_id, unlist(required_families), sep = "::")
    missing_required <- setdiff(required_columns, family_columns)
    if (length(missing_required)) {
      candidate_release_family_coherent[, candidate_index] <- FALSE
    } else {
      candidate_release_family_coherent[, candidate_index] <- rowSums(
        coherent[, required_columns, drop = FALSE]
      ) == length(required_columns)
    }
  }
  candidate_direct[, candidate_index] <- top_two_mean(direct_matrix)
  candidate_local[, candidate_index] <- top_two_mean(local_matrix)
  positive_signal <- top_two_mean(combined_matrix)

  anti_genes <- intersect(soft_anti_genes[[candidate_id]], available)
  anti_signal <- matrix(0, nrow = n_obs, ncol = length(anti_genes))
  if (length(anti_genes)) {
    for (gene_index in seq_along(anti_genes)) {
      values <- gene_score(anti_genes[[gene_index]])
      anti_signal[, gene_index] <- direct_weight * values$direct +
        local_weight * values$local
    }
  }
  soft_anti <- if (length(anti_genes)) top_two_mean(anti_signal) else numeric(n_obs)
  direct_anti_count <- integer(n_obs)
  direct_anti_family_count <- integer(n_obs)
  hard_contradiction <- rep(FALSE, n_obs)
  hard_sources <- hard_anti_sources[[candidate_id]]
  if (length(hard_sources)) {
    for (source in hard_sources) {
      hard_family_ids <- names(positive_families)[
        startsWith(names(positive_families), paste0(source, "::"))
      ]
      source_gene_count <- integer(n_obs)
      source_family_count <- integer(n_obs)
      for (hard_family_id in hard_family_ids) {
        hard_genes <- intersect(positive_families[[hard_family_id]], available)
        family_gene_count <- integer(n_obs)
        for (gene in hard_genes) {
          family_gene_count <- family_gene_count +
            as.integer(gene_score(gene)$direct >= 0.15)
        }
        source_gene_count <- source_gene_count + family_gene_count
        source_family_count <- source_family_count +
          as.integer(family_gene_count >= 1L)
      }
      source_hard <- source_gene_count >= 2L & source_family_count >= 2L
      hard_contradiction <- hard_contradiction | source_hard
      direct_anti_count <- pmax(direct_anti_count, source_gene_count)
      direct_anti_family_count <- pmax(
        direct_anti_family_count, source_family_count
      )
    }
  }
  candidate_soft_anti[, candidate_index] <- soft_anti
  candidate_direct_anti_count[, candidate_index] <- direct_anti_count
  candidate_direct_anti_family_count[, candidate_index] <-
    direct_anti_family_count
  candidate_hard_contradiction[, candidate_index] <- hard_contradiction
  candidate_score[, candidate_index] <- positive_signal - anti_weight * soft_anti
}

# Candidate-local seeds are computed independently; aggregate winners never
# suppress a coherent alternative program.
candidate_seed <- candidate_identity_core_coherent &
  candidate_score >= 0.04 & !candidate_hard_contradiction
candidate_local_seed_fraction <- matrix(
  0, nrow = n_obs, ncol = n_candidates,
  dimnames = list(colnames(obj), candidate_ids)
)
for (candidate_index in seq_along(candidate_ids)) {
  candidate_local_seed_fraction[, candidate_index] <- rowMeans(matrix(
    as.numeric(candidate_seed[, candidate_index])[neighbors],
    nrow = nrow(neighbors), ncol = ncol(neighbors)
  ))
}

# Candidate-specific stability is assessed at the selected resolution and its
# two contract-selected nearest neighbours. A weak stable program remains visible.
partition_list <- split(partitions, partitions$resolution)
cross_resolution_count <- matrix(
  0L, nrow = n_obs, ncol = n_candidates,
  dimnames = list(colnames(obj), candidate_ids)
)
stability_rows <- list()
stability_i <- 1L
for (resolution in sort(names(partition_list))) {
  partition <- partition_list[[resolution]]
  partition <- partition[match(colnames(obj), cell_id)]
  if (anyNA(partition$cell_id)) stop("partition resolution does not cover all observations")
  group <- paste(partition$boundary_id, partition$cluster, sep = "::")
  for (candidate_index in seq_along(candidate_ids)) {
    dt <- data.table(
      group = group,
      seed = candidate_seed[, candidate_index],
      score = candidate_score[, candidate_index]
    )
    summary <- dt[, .(
      n = .N,
      seed_fraction = mean(seed),
      mean_score = mean(score)
    ), by = group]
    mapped <- summary[match(group, summary$group)]
    group_stable <- mapped$n >= 20L & mapped$seed_fraction >= 0.03 &
      mapped$mean_score >= 0.02
    rare_spatial_stable <- candidate_seed[, candidate_index] &
      candidate_local_seed_fraction[, candidate_index] >= 0.10
    stable <- group_stable | rare_spatial_stable
    cross_resolution_count[, candidate_index] <-
      cross_resolution_count[, candidate_index] + as.integer(stable)
    summary[, `:=`(
      resolution = resolution,
      candidate_id = candidate_ids[[candidate_index]]
    )]
    stability_rows[[stability_i]] <- summary
    stability_i <- stability_i + 1L
  }
}
write_tsv_atomic(
  rbindlist(stability_rows),
  file.path(out_dir, "tables", "cross_resolution_candidate_stability.tsv.gz")
)

# Normalize within candidate only after absolute evidence and contradictions
# are frozen. This is used for overlap discrimination, never presence.
candidate_normalized <- candidate_score
for (candidate_index in seq_along(candidate_ids)) {
  positive <- candidate_score[, candidate_index][candidate_score[, candidate_index] > 0]
  scale <- if (length(positive)) {
    as.numeric(quantile(positive, 0.95, names = FALSE, type = 8))
  } else {
    1
  }
  if (!is.finite(scale) || scale <= 0) scale <- 1
  candidate_normalized[, candidate_index] <- candidate_score[, candidate_index] / scale
}

n_feature_columns <- grep("^nFeature_", names(metadata), value = TRUE)
technical_flag <- rep(FALSE, n_obs)
if (length(n_feature_columns)) {
  features_per_observation <- as.numeric(metadata[[n_feature_columns[[1L]]]])
  threshold <- quantile(features_per_observation, 0.01, na.rm = TRUE, names = FALSE)
  technical_flag <- !is.finite(features_per_observation) |
    features_per_observation <= threshold
}

score_path <- file.path(out_dir, "tables", "observation_lineage_scores.tsv.gz")
if (!grid_evidence_only) {
  score_rows <- vector("list", n_candidates)
  for (candidate_index in seq_along(candidate_ids)) {
  candidate_id <- candidate_ids[[candidate_index]]
  candidate <- candidate_meta[[candidate_id]]
  family_columns <- grep(
    paste0("^", candidate_id, "::"),
    colnames(family_direct), value = TRUE
  )
  positive_names <- apply(family_active[, family_columns, drop = FALSE], 1L, function(value) {
    paste(sub("^[^:]+::", "", family_columns[value]), collapse = ";")
  })
  score_rows[[candidate_index]] <- data.table(
    cell_id = colnames(obj),
    source_boundary = selected$boundary_id,
    source_cluster = selected$cluster,
    neighbor_1_boundary = partition_by_role$neighbor_1$boundary_id,
    neighbor_1_cluster = partition_by_role$neighbor_1$cluster,
    neighbor_2_boundary = partition_by_role$neighbor_2$boundary_id,
    neighbor_2_cluster = partition_by_role$neighbor_2$cluster,
    candidate_id = candidate_id,
    candidate_role = or_value(candidate$candidate_role, "exploratory"),
    release_broad_label = or_value(candidate$release_broad_label, ""),
    release_fine_label = or_value(candidate$release_fine_label, ""),
    specificity_priority = or_value(candidate$specificity_priority, 0),
    direct_signal = candidate_direct[, candidate_index],
    local_signal = candidate_local[, candidate_index],
    ambient_suspect = candidate_direct[, candidate_index] < 0.03 &
      candidate_local[, candidate_index] >= 0.12,
    program_score = candidate_score[, candidate_index],
    normalized_evidence = candidate_normalized[, candidate_index],
    positive_family_count = candidate_positive_family_count[, candidate_index],
    positive_families = positive_names,
    positive_gene_count = candidate_positive_gene_count[, candidate_index],
    family_coherent = candidate_family_coherent[, candidate_index],
    identity_core_coherent =
      candidate_identity_core_coherent[, candidate_index],
    identity_core_direct = candidate_identity_core_direct[, candidate_index],
    release_family_coherent =
      candidate_release_family_coherent[, candidate_index],
    soft_anti_score = candidate_soft_anti[, candidate_index],
    direct_anti_gene_count = candidate_direct_anti_count[, candidate_index],
    direct_anti_family_count =
      candidate_direct_anti_family_count[, candidate_index],
    hard_contradiction = candidate_hard_contradiction[, candidate_index],
    candidate_seed = candidate_seed[, candidate_index],
    local_seed_fraction = candidate_local_seed_fraction[, candidate_index],
    cross_resolution_support_count = cross_resolution_count[, candidate_index],
    whole_subcluster_inherit = FALSE,
    technical_flag = technical_flag,
    x = coordinates[, 1L],
    y = coordinates[, 2L]
  )
  }
  scores <- rbindlist(score_rows, use.names = TRUE)
  setorder(scores, cell_id, candidate_id)
  write_tsv_atomic(scores, score_path)
}

# Group-level DEG/pseudobulk and spatial summaries supplement, but never
# replace, observation evidence.
group_rows <- list()
group_i <- 1L
for (resolution in sort(names(partition_list))) {
  partition <- partition_list[[resolution]]
  partition <- partition[match(colnames(obj), cell_id)]
  source_groups <- paste(partition$boundary_id, partition$cluster, sep = "::")
  resolution_role <- paste(unique(partition$resolution_role), collapse = ";")
  group_levels <- unique(source_groups)
  group_factor <- factor(source_groups, levels = group_levels)
  group_codes <- as.integer(group_factor)
  n_groups <- length(group_levels)
  group_n <- tabulate(group_codes, nbins = n_groups)
  neighbor_group_codes <- matrix(
    group_codes[neighbors], nrow = nrow(neighbors), ncol = ncol(neighbors)
  )
  same_group_neighbor_fraction <- rowMeans(
    neighbor_group_codes == group_codes
  )
  group_spatial_connectivity <- as.numeric(rowsum(
    as.numeric(same_group_neighbor_fraction >= 0.50),
    group_factor, reorder = FALSE
  )) / group_n
  first_member <- match(group_levels, source_groups)
  membership <- sparseMatrix(
    i = seq_len(n_obs), j = group_codes, x = 1,
    dims = c(n_obs, n_groups)
  )
  group_mean_matrix <- function(value) {
    sweep(rowsum(value, group_factor, reorder = FALSE), 1L, group_n, "/")
  }
  seed_fraction <- group_mean_matrix(candidate_seed * 1)
  identity_core_fraction <- group_mean_matrix(
    candidate_identity_core_coherent * 1
  )
  identity_core_direct_fraction <- group_mean_matrix(
    candidate_identity_core_direct * 1
  )
  coherent_fraction <- group_mean_matrix(candidate_family_coherent * 1)
  release_coherent_fraction <- group_mean_matrix(
    (candidate_family_coherent & candidate_release_family_coherent) * 1
  )
  contradiction_fraction <- group_mean_matrix(candidate_hard_contradiction * 1)
  mean_program_score <- group_mean_matrix(candidate_score)
  spatial_support_fraction <- group_mean_matrix(
    (candidate_local_seed_fraction >= 0.03) * 1
  )
  stable_fraction <- group_mean_matrix(
    (cross_resolution_count >= 2L) * 1
  )
  summarize_gene_block <- function(gene_index) {
    if (!length(gene_index)) {
      return(list(
        deg = rep(0, n_groups),
        pseudobulk = rep(0, n_groups),
        detection = rep(0, n_groups)
      ))
    }
    block <- expr[gene_index, , drop = FALSE]
    group_sums <- as.matrix(block %*% membership)
    total_sums <- Matrix::rowSums(block)
    group_means <- sweep(group_sums, 2L, group_n, "/")
    rest_n <- n_obs - group_n
    rest_sums <- outer(total_sums, rep(1, n_groups)) - group_sums
    rest_means <- sweep(rest_sums, 2L, pmax(rest_n, 1L), "/")
    if (any(rest_n == 0L)) rest_means[, rest_n == 0L] <- 0
    detected <- as.numeric(Matrix::colSums(block > 0) > 0)
    deg <- colMeans(log2((group_means + 0.01) / (rest_means + 0.01)))
    # A one-cluster partition has no biological one-vs-rest contrast. Preserve
    # its absolute pseudobulk/detection evidence and mark DEG as neutral rather
    # than manufacturing enrichment against an empty background.
    deg[rest_n < 1L] <- 0
    list(
      deg = deg,
      pseudobulk = colSums(group_sums),
      detection = as.numeric(crossprod(membership, detected)) / group_n
    )
  }
  for (candidate_index in seq_along(candidate_ids)) {
    candidate_id <- candidate_ids[[candidate_index]]
    positive_genes <- unique(unlist(
      positive_families[startsWith(names(positive_families), paste0(candidate_id, "::"))]
    ))
    positive_genes <- intersect(positive_genes, available)
    positive_family_keys <- names(positive_families)[
      startsWith(names(positive_families), paste0(candidate_id, "::"))
    ]
    available_positive_family_count <- sum(vapply(
      positive_families[positive_family_keys],
      function(genes) any(genes %chin% available), logical(1)
    ))
    family_prevalence <- group_mean_matrix(
      family_active[, positive_family_keys, drop = FALSE] * 1
    )
    identity_floor <- max(
      0.03,
      as.numeric(or_value(
        candidate_meta[[candidate_id]]$minimum_identity_core_fraction,
        0.03
      )) / 2
    )
    group_positive_family_supported_count <- rowSums(
      family_prevalence >= identity_floor
    )
    group_positive_family_mean_fraction <- rowMeans(family_prevalence)
    required_group_families <- as.character(unlist(
      candidate_meta[[candidate_id]]$required_positive_families
    ))
    group_required_positive_families_pass <- rep(TRUE, n_groups)
    if (length(required_group_families)) {
      required_group_columns <- paste(
        candidate_id, required_group_families, sep = "::"
      )
      if (length(setdiff(required_group_columns, colnames(family_prevalence)))) {
        group_required_positive_families_pass[] <- FALSE
      } else {
        group_required_positive_families_pass <- rowSums(
          family_prevalence[, required_group_columns, drop = FALSE] >=
            identity_floor
        ) == length(required_group_columns)
      }
    }
    positive_index <- unname(feature_lookup[positive_genes])
    positive_summary <- summarize_gene_block(positive_index)
    anti_genes <- intersect(soft_anti_genes[[candidate_id]], available)
    anti_index <- unname(feature_lookup[anti_genes])
    anti_summary <- summarize_gene_block(anti_index)
    group_rows[[group_i]] <- data.table(
      resolution = resolution,
      resolution_role = resolution_role,
      source_boundary = partition$boundary_id[first_member],
      source_cluster = partition$cluster[first_member],
      candidate_id = candidate_id,
      n_observations = group_n,
      observation_seed_fraction = seed_fraction[, candidate_index],
      observation_identity_core_fraction =
        identity_core_fraction[, candidate_index],
      observation_identity_core_direct_fraction =
        identity_core_direct_fraction[, candidate_index],
      observation_coherent_fraction = coherent_fraction[, candidate_index],
      observation_release_family_coherent_fraction =
        release_coherent_fraction[, candidate_index],
      hard_contradiction_fraction = contradiction_fraction[, candidate_index],
      mean_program_score = mean_program_score[, candidate_index],
      available_positive_gene_count = length(positive_genes),
      available_positive_family_count = available_positive_family_count,
      group_positive_family_supported_count =
        group_positive_family_supported_count,
      group_positive_family_mean_fraction =
        group_positive_family_mean_fraction,
      group_required_positive_families_pass =
        group_required_positive_families_pass,
      positive_marker_detection_fraction = positive_summary$detection,
      positive_marker_pseudobulk_sum = positive_summary$pseudobulk,
      marker_deg_log2fc_mean = positive_summary$deg,
      anti_marker_detection_fraction = anti_summary$detection,
      anti_marker_pseudobulk_sum = anti_summary$pseudobulk,
      anti_marker_deg_log2fc_mean = anti_summary$deg,
      spatial_local_support_fraction = spatial_support_fraction[, candidate_index],
      spatial_group_connectivity_fraction = group_spatial_connectivity,
      cross_resolution_stable_fraction = stable_fraction[, candidate_index]
    )
    group_i <- group_i + 1L
  }
}
write_tsv_atomic(
  rbindlist(group_rows),
  file.path(out_dir, "tables", "cluster_candidate_multichannel_evidence.tsv.gz")
)

# Extract stable-program candidates without a fixed lineage whitelist. This
# uses sparse pseudobulk/detection aggregates at the selected and neighboring
# resolutions, then removes technical/state-only genes before cross-resolution
# matching by the controller.
technical_genes <- toupper(c(
  "FOS", "JUN", "JUNB", "DUSP1", "HSPA1A", "HSPA1B", "HSP90AA1",
  "HIF1A", "VEGFA", "BNIP3", "CA9", "EGLN3",
  "MKI67", "TOP2A", "UBE2C", "CENPF", "PCNA",
  "DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2"
))
catalog_markers <- unique(unlist(positive_families))
program_rows <- list()
program_i <- 1L
detection_matrix <- expr > 0
total_sum <- Matrix::rowSums(expr)
total_detection <- Matrix::rowSums(detection_matrix)
for (resolution in sort(names(partition_list))) {
  partition <- partition_list[[resolution]]
  partition <- partition[match(colnames(obj), cell_id)]
  group <- factor(paste(partition$boundary_id, partition$cluster, sep = "::"))
  # Unmodeled-program discovery is differential by definition. A one-cluster
  # partition remains valid for broad absolute evidence but cannot generate a
  # cluster-vs-rest program.
  if (nlevels(group) < 2L) next
  group_code <- as.integer(group)
  design <- Matrix::sparse.model.matrix(~ 0 + group)
  group_sum <- expr %*% design
  group_detection <- detection_matrix %*% design
  group_n <- as.numeric(table(group))
  for (group_index in seq_len(ncol(group_sum))) {
    n_group <- group_n[[group_index]]
    n_rest <- n_obs - n_group
    if (n_group < 20L || n_rest < 20L) next
    mean_group <- as.numeric(group_sum[, group_index]) / n_group
    mean_rest <- (total_sum - as.numeric(group_sum[, group_index])) / n_rest
    pct_group <- as.numeric(group_detection[, group_index]) / n_group
    pct_rest <- (total_detection - as.numeric(group_detection[, group_index])) / n_rest
    logfc <- log2((mean_group + 0.01) / (mean_rest + 0.01))
    allowed <- !grepl("^(RPL|RPS|MRPL|MRPS|MT-)", features_upper) &
      !features_upper %chin% technical_genes &
      !grepl("^COL[0-9]", features_upper) &
      pct_group >= 0.05 & logfc >= 0.25
    ordered <- order(logfc, pct_group - pct_rest, decreasing = TRUE)
    ordered <- ordered[allowed[ordered]]
    top <- head(ordered, 30L)
    if (length(top) < 2L) next
    genes <- features_upper[top]
    catalog_overlap <- mean(genes %chin% catalog_markers)
    catalog_matches <- lapply(candidate_ids, function(candidate_id) {
      family_keys <- names(positive_families)[
        startsWith(names(positive_families), paste0(candidate_id, "::"))
      ]
      overlap_by_family <- vapply(
        positive_families[family_keys],
        function(markers) length(intersect(genes, markers)),
        integer(1)
      )
      list(
        candidate_id = candidate_id,
        family_count = sum(overlap_by_family > 0L),
        gene_count = length(intersect(
          genes, unique(unlist(positive_families[family_keys]))
        ))
      )
    })
    catalog_match_order <- order(
      -vapply(catalog_matches, `[[`, integer(1), "family_count"),
      -vapply(catalog_matches, `[[`, integer(1), "gene_count"),
      vapply(catalog_matches, `[[`, character(1), "candidate_id")
    )
    best_catalog_match <- catalog_matches[[catalog_match_order[[1L]]]]
    group_members <- which(group == levels(group)[[group_index]])
    neighbor_group <- matrix(
      group_code[neighbors[group_members, , drop = FALSE]],
      nrow = length(group_members), ncol = ncol(neighbors)
    )
    same_neighbor <- rowMeans(neighbor_group == group_code[group_members])
    spatial_coherent <- mean(same_neighbor >= 0.50) >= 0.70
    group_parts <- strsplit(levels(group)[[group_index]], "::", fixed = TRUE)[[1L]]
    program_rows[[program_i]] <- data.table(
      program_id = paste0("res", resolution, "__", safe_name(levels(group)[[group_index]])),
      resolution = resolution,
      source_boundary = group_parts[[1L]],
      source_cluster = group_parts[[2L]],
      n_observations = n_group,
      genes = paste(genes, collapse = ";"),
      coexpressed_gene_count = length(genes),
      mean_top_log2fc = mean(logfc[top]),
      mean_detection_difference = mean(pct_group[top] - pct_rest[top]),
      catalog_marker_overlap_fraction = catalog_overlap,
      best_catalog_candidate_id = best_catalog_match$candidate_id,
      best_catalog_overlap_gene_count = best_catalog_match$gene_count,
      best_catalog_overlap_family_count = best_catalog_match$family_count,
      spatially_coherent = spatial_coherent,
      excluded_program_classes = "",
      candidate_status = ifelse(
        catalog_overlap < 0.30 && spatial_coherent,
        "unmodeled_program_seed", "modeled_or_incoherent"
      )
    )
    program_i <- program_i + 1L
  }
}
program_table <- if (length(program_rows)) rbindlist(program_rows, fill = TRUE) else data.table(
  program_id = character(), resolution = character(), source_boundary = character(),
  source_cluster = character(), n_observations = integer(), genes = character(),
  coexpressed_gene_count = integer(), mean_top_log2fc = numeric(),
  mean_detection_difference = numeric(), catalog_marker_overlap_fraction = numeric(),
  best_catalog_candidate_id = character(),
  best_catalog_overlap_gene_count = integer(),
  best_catalog_overlap_family_count = integer(),
  spatially_coherent = logical(), excluded_program_classes = character(),
  candidate_status = character()
)
write_tsv_atomic(
  program_table,
  file.path(out_dir, "tables", "resolution_deg_coexpression_programs.tsv.gz")
)

manifest <- list(
  status = "PASS",
  schema_version = "2.2",
  controller_version = "2.2.0",
  stage = "label_blind_observation_lineage_scoring",
  scorer = list(path = script_path, sha256 = sha256(script_path)),
  input_rds = input_rds,
  profile = profile_path,
  catalog = catalog_path,
  threshold_registry = list(
    path = threshold_registry_path,
    sha256 = sha256(threshold_registry_path)
  ),
  partitions = partition_path,
  analysis_membership = analysis_membership_path,
  n_observations = n_obs,
  candidate_universe = candidate_ids,
  historical_labels_read = FALSE,
  assay = assay_name,
  layer = layer_name,
  parameters = list(
    seed = seed,
    workers = workers,
    knn_k = nn_k,
    gene_scale = "query_nonzero_q95_capped",
    family_aggregation = "mean_top_two_available_genes",
    direct_weight = direct_weight,
    local_weight = local_weight,
    anti_weight = anti_weight,
    family_active_threshold = family_active_threshold,
    local_gene_detection_fraction = local_gene_detection_fraction,
    family_coherence = "direct>=2 OR direct>=1+local>=2 OR local>=3",
    anti_policy = paste(
      "single_or_local_anti_soft_penalty;",
      "direct_multigene_multifamily_hard_contradiction",
      sep = ""
    ),
    candidate_seed_policy = paste(
      "declared_identity_core_families;",
      "per_observation_seed_any_core;group_required_families_all;",
      "candidate_local_independent_of_aggregate_winner",
      sep = ""
    )
    ,
    grid_evidence_only = grid_evidence_only
  ),
  outputs = list(
    observation_scores = if (file.exists(score_path)) normalizePath(score_path) else NULL,
    cluster_multichannel_evidence = normalizePath(file.path(
      out_dir, "tables", "cluster_candidate_multichannel_evidence.tsv.gz"
    )),
    resolution_deg_coexpression_programs = normalizePath(file.path(
      out_dir, "tables", "resolution_deg_coexpression_programs.tsv.gz"
    ))
  )
)
write_json(
  manifest, file.path(out_dir, "observation_scoring_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE
)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
writeLines("status\tPASS", file.path(out_dir, "RUN_COMPLETE.tsv"))
