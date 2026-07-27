#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(RANN)
  library(igraph)
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
required <- c("scores", "cluster-evidence", "catalog", "out")
missing <- required[!required %in% names(args)]
if (length(missing)) stop("missing required arguments: ", paste(missing, collapse = ", "))

score_path <- normalizePath(args$scores, mustWork = TRUE)
cluster_evidence_path <- normalizePath(args[["cluster-evidence"]], mustWork = TRUE)
catalog_path <- normalizePath(args$catalog, mustWork = TRUE)
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("cannot resolve local-subset script path")
script_path <- normalizePath(
  sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE
)
sha256 <- function(path) {
  output <- system2("sha256sum", shQuote(path), stdout = TRUE)
  strsplit(output[[1L]], "\\s+")[[1L]][[1L]]
}
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
writeback_policy <- threshold_registry$observation_writeback_policy
subset_support_threshold <- as.numeric(
  writeback_policy$supported_subset_min_lineage_supported_fraction
)
subset_margin_threshold <- as.numeric(
  writeback_policy$supported_subset_min_purity_margin
)
maximum_contradiction_threshold <- as.numeric(
  writeback_policy$maximum_contradiction_fraction
)
parent_support_threshold <- as.numeric(
  writeback_policy$whole_subcluster_min_lineage_supported_fraction
)
minimum_component_n <- as.integer(
  threshold_registry$local_subset_policy$minimum_component_members
)
out_dir <- args$out
release_level <- if ("release-level" %in% names(args)) args[["release-level"]] else "broad"
if (!release_level %in% c("broad", "fine", "all")) stop("invalid --release-level")
workers <- if ("workers" %in% names(args)) {
  as.integer(args$workers)
} else {
  as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "1"))
}
if (is.na(workers) || workers < 1L) stop("--workers must be a positive integer")
setDTthreads(max(1L, workers))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  if (grepl("\\.gz$", path)) {
    fread(cmd = paste("gzip -dc", shQuote(path)), colClasses = "character")
  } else {
    fread(path, colClasses = "character")
  }
}
truth <- function(value) tolower(as.character(value)) %in% c("1", "true", "yes", "pass")
safe <- function(value) gsub("[^A-Za-z0-9]+", "_", value)
or_value <- function(value, fallback) {
  if (is.null(value) || !length(value) || is.na(value[[1L]])) fallback else value[[1L]]
}
or_values <- function(value, fallback = list()) {
  if (is.null(value) || !length(value) || all(is.na(value))) fallback else value
}
catalog_candidates <- function(catalog) {
  rows <- catalog$candidate_boundaries
  broad_by_label <- list()
  for (row in rows) {
    if (
      identical(tolower(as.character(or_value(row$candidate_role, ""))), "broad") &&
        nzchar(as.character(or_value(row$release_broad_label, "")))
    ) {
      broad_by_label[[as.character(row$release_broad_label)]] <- row
    }
  }
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
      parent_candidate <- broad_by_label[[parent]]
      parent_candidate_id <- if (is.null(parent_candidate)) "" else {
        as.character(or_value(parent_candidate$candidate_id, ""))
      }
      item$candidate_role <- "fine"
      item$release_broad_label <- parent
      item$release_fine_label <- fine_label
      item$parent_broad_label <- parent
      item$writeback_strategy <- as.character(or_value(
        item$writeback_strategy, "supported_subset_with_parent_lock"
      ))
      item$specificity_priority <- as.numeric(or_value(item$specificity_priority, 70))
      item$hard_anti_families <- or_values(item$hard_anti_families)
      item$soft_anti_families <- or_values(item$soft_anti_families)
      item$context_requirements <- or_values(item$context_requirements, context)
      item$required_positive_families <- or_values(
        item$required_positive_families,
        list("parent_identity", "fine_discriminator")
      )
      item$seed_required_positive_families <- or_values(
        item$seed_required_positive_families,
        list("fine_discriminator")
      )
      item$formal_context_evidence_required <- isTRUE(
        item$formal_context_evidence_required
      ) || nzchar(as.character(or_value(item$context_gate, ""))) || (
        !is.null(parent_candidate) &&
          isTRUE(parent_candidate$formal_context_evidence_required)
      )
      item$context_evidence_candidate_id <- as.character(or_value(
        item$context_evidence_candidate_id, parent_candidate_id
      ))
      item$review_required <- TRUE
      item$parent_broad_writeback_strategy <- as.character(or_value(
        parent_candidate$writeback_strategy, ""
      ))
      rows[[length(rows) + 1L]] <- item
      represented <- c(represented, fine_label)
    }
  }
  rows
}

scores <- read_tsv(score_path)
if (!("source-boundary" %in% names(args)) || !("source-cluster" %in% names(args))) {
  stop("local subset derivation requires --source-boundary and --source-cluster")
}
source_boundary_filter <- as.character(args[["source-boundary"]])
source_cluster_filter <- as.character(args[["source-cluster"]])
numeric_columns <- c(
  "specificity_priority", "direct_signal", "local_signal",
  "program_score", "normalized_evidence",
  "positive_family_count", "positive_gene_count", "soft_anti_score",
  "direct_anti_gene_count", "direct_anti_family_count", "local_seed_fraction",
  "cross_resolution_support_count", "x", "y"
)
for (column in intersect(numeric_columns, names(scores))) {
  set(scores, j = column, value = as.numeric(scores[[column]]))
}
required_score <- c(
  "cell_id", "source_boundary", "source_cluster", "candidate_id",
  "neighbor_1_boundary", "neighbor_1_cluster",
  "neighbor_2_boundary", "neighbor_2_cluster",
  "candidate_role", "release_broad_label", "program_score",
  "normalized_evidence", "positive_family_count", "family_coherent",
  "positive_families",
  "release_family_coherent", "direct_anti_gene_count",
  "direct_anti_family_count",
  "hard_contradiction", "identity_core_coherent",
  "identity_core_direct", "candidate_seed",
  "local_seed_fraction", "cross_resolution_support_count", "x", "y"
)
if (!all(required_score %in% names(scores))) {
  stop("score table lacks: ", paste(setdiff(required_score, names(scores)), collapse = ", "))
}
scores <- scores[
  source_boundary == source_boundary_filter &
    source_cluster == source_cluster_filter
]
if (!nrow(scores)) stop("requested local mixed subcluster is absent from scores")
if (scores[, uniqueN(candidate_id), by = cell_id][, uniqueN(V1)] != 1L) {
  stop("score table is not a complete observation x candidate product")
}
cluster_evidence <- read_tsv(cluster_evidence_path)
cluster_evidence <- cluster_evidence[
  source_boundary == source_boundary_filter &
    source_cluster == source_cluster_filter
]
if (!nrow(cluster_evidence)) {
  stop("requested local mixed subcluster is absent from cluster evidence")
}
cluster_numeric <- c(
  "n_observations", "observation_seed_fraction",
  "observation_identity_core_fraction",
  "observation_identity_core_direct_fraction",
  "observation_coherent_fraction",
  "observation_release_family_coherent_fraction",
  "available_positive_family_count",
  "group_positive_family_supported_count",
  "hard_contradiction_fraction",
  "mean_program_score", "positive_marker_detection_fraction",
  "positive_marker_pseudobulk_sum", "marker_deg_log2fc_mean",
  "anti_marker_detection_fraction", "anti_marker_pseudobulk_sum",
  "anti_marker_deg_log2fc_mean", "spatial_local_support_fraction",
  "cross_resolution_stable_fraction"
)
for (column in intersect(cluster_numeric, names(cluster_evidence))) {
  set(cluster_evidence, j = column, value = as.numeric(cluster_evidence[[column]]))
}
required_cluster <- c(
  "resolution_role", "source_boundary", "source_cluster", "candidate_id",
  "n_observations", "observation_seed_fraction",
  "observation_identity_core_fraction",
  "observation_identity_core_direct_fraction",
  "observation_coherent_fraction", "hard_contradiction_fraction",
  "mean_program_score", "positive_marker_detection_fraction",
  "marker_deg_log2fc_mean", "anti_marker_deg_log2fc_mean",
  "spatial_local_support_fraction", "cross_resolution_stable_fraction"
)
if (!all(required_cluster %in% names(cluster_evidence))) {
  stop(
    "cluster evidence lacks: ",
    paste(setdiff(required_cluster, names(cluster_evidence)), collapse = ", ")
  )
}

catalog <- fromJSON(catalog_path, simplifyVector = FALSE)
catalog_candidate_rows <- catalog_candidates(catalog)
candidate_ids <- vapply(catalog_candidate_rows, `[[`, character(1), "candidate_id")
if (!setequal(candidate_ids, unique(scores$candidate_id))) {
  stop("score candidate universe differs from the catalog")
}
catalog_by_id <- setNames(catalog_candidate_rows, candidate_ids)
effective_broad_writeback_strategy <- function(candidate_id) {
  candidate <- catalog_by_id[[candidate_id]]
  inherited <- as.character(or_value(
    candidate$parent_broad_writeback_strategy, ""
  ))
  if (nzchar(inherited)) return(inherited)
  if (identical(as.character(candidate$candidate_role), "fine")) {
    parent <- as.character(or_value(candidate$release_broad_label, ""))
    matches <- Filter(function(item) {
      identical(as.character(or_value(item$candidate_role, "")), "broad") &&
        identical(as.character(or_value(item$release_broad_label, "")), parent)
    }, catalog_candidate_rows)
    if (length(matches)) {
      strict_parent <- Filter(function(item) {
        !nzchar(as.character(or_value(item$release_fine_label, "")))
      }, matches)
      chosen <- if (length(strict_parent)) strict_parent[[1L]] else matches[[1L]]
      return(as.character(or_value(chosen$writeback_strategy, "")))
    }
  }
  as.character(or_value(candidate$writeback_strategy, ""))
}
context_ok <- setNames(rep(FALSE, length(candidate_ids)), candidate_ids)
accepted_context <- character()
if ("context-evidence" %in% names(args)) {
  context_table <- read_tsv(normalizePath(args[["context-evidence"]], mustWork = TRUE))
  if (!all(c("candidate_id", "status") %in% names(context_table))) {
    stop("context evidence must contain candidate_id and status")
  }
  accepted_context <- context_table[
    tolower(status) %chin% c("supported", "pass"), candidate_id
  ]
}
for (candidate_id in candidate_ids) {
  evidence_id <- as.character(or_value(
    catalog_by_id[[candidate_id]]$context_evidence_candidate_id,
    candidate_id
  ))
  context_ok[[candidate_id]] <- candidate_id %chin% accepted_context ||
    evidence_id %chin% accepted_context
}
candidate_broad <- setNames(vapply(catalog_candidate_rows, function(candidate) {
  as.character(or_value(candidate$release_broad_label, ""))
}, character(1)), candidate_ids)
candidate_role <- setNames(vapply(catalog_candidate_rows, function(candidate) {
  as.character(or_value(candidate$candidate_role, "exploratory"))
}, character(1)), candidate_ids)
scores[, effective_broad_label := unname(candidate_broad[candidate_id])]
scores[, effective_candidate_role := unname(candidate_role[candidate_id])]
release_eligible <- vapply(candidate_ids, function(candidate_id) {
  candidate <- catalog_by_id[[candidate_id]]
  role <- as.character(candidate$candidate_role)
  broad <- as.character(candidate$release_broad_label)
  requires_context <- isTRUE(candidate$formal_context_evidence_required)
  role %in% c("broad", "fine") && length(broad) && nzchar(broad) &&
    (!requires_context || context_ok[[candidate_id]])
}, logical(1))

identity_core_mask <- function(candidate_rows, candidate_id) {
  base_core <- truth(candidate_rows$identity_core_coherent)
  activated <- truth(candidate_rows$candidate_seed) |
    candidate_rows$program_score >= 0.02 |
    candidate_rows$direct_signal > 0
  base_core & activated
}

minimum_identity_core_fraction <- function(candidate_id) {
  configured <- as.numeric(or_value(
    catalog_by_id[[candidate_id]]$minimum_identity_core_fraction,
    0.03
  ))
  if (!is.finite(configured)) configured <- 0.03
  max(0.03, min(0.50, configured))
}

# The strict subset thresholds validate a generated group, never one cellbin.
# Candidate-local admission begins from a declared identity family. Generic
# support genes are not allowed to connect identity cores transitively. Sparse
# non-core observations may inherit a label only from an independently
# validated whole or neighboring-resolution expression subcluster.
scores[, hard_block := (
  truth(hard_contradiction) & direct_anti_gene_count >= 2 &
    direct_anti_family_count >= 2
)]
scores[, anchor := truth(family_coherent) & positive_family_count >= 1 &
         !hard_block]
scores[, seed := truth(candidate_seed) & anchor &
         program_score >= 0.04]
scores[, strict_support := seed & positive_family_count >= 2 &
         truth(release_family_coherent)]
scores[, expanded_seed := anchor & (
  seed | program_score >= -0.02 | local_seed_fraction >= 0.03 |
    cross_resolution_support_count >= 1
)]

group_family_evidence <- function(candidate_rows, candidate_id) {
  n <- nrow(candidate_rows)
  if (!n) {
    return(list(
      pass = FALSE, supported_families = "", family_prevalence = ""
    ))
  }
  family_values <- if ("positive_families" %in% names(candidate_rows)) {
    candidate_rows$positive_families
  } else {
    rep("", n)
  }
  per_row <- strsplit(
    ifelse(is.na(family_values), "", as.character(family_values)),
    ";", fixed = TRUE
  )
  family_counts <- table(unlist(lapply(per_row, function(values) {
    unique(values[nzchar(values)])
  }), use.names = FALSE))
  prevalence <- family_counts / n
  supported <- sort(names(prevalence)[prevalence >= 0.03])
  required <- as.character(unlist(or_values(
    catalog_by_id[[candidate_id]]$required_positive_families
  ), use.names = FALSE))
  required <- unique(required[nzchar(required)])
  prevalence_text <- if (length(prevalence)) {
    paste(
      paste0(
        names(prevalence), "=",
        formatC(as.numeric(prevalence), digits = 6L, format = "f")
      ),
      collapse = ";"
    )
  } else {
    ""
  }
  list(
    pass = length(supported) >= 2L && all(required %chin% supported),
    supported_families = paste(supported, collapse = ";"),
    family_prevalence = prevalence_text
  )
}

aggregate_row_for <- function(
  resolution_role, source_boundary, source_cluster, candidate_id
) {
  role_value <- resolution_role
  boundary_value <- source_boundary
  cluster_value <- source_cluster
  candidate_value <- candidate_id
  value <- cluster_evidence[
    resolution_role == role_value &
      source_boundary == boundary_value &
      source_cluster == cluster_value &
      candidate_id == candidate_value
  ]
  if (nrow(value) != 1L) return(NULL)
  value[1L]
}

aggregate_score <- function(row) {
  if (is.null(row) || !nrow(row)) return(-Inf)
  row$marker_deg_log2fc_mean +
    0.25 * row$observation_coherent_fraction +
    0.25 * row$observation_seed_fraction +
    0.10 * row$mean_program_score -
    0.25 * max(0, row$anti_marker_deg_log2fc_mean)
}

canonical_cluster_challenger_supported <- function(row, candidate_id) {
  if (is.null(row) || !nrow(row)) return(FALSE)
  candidate <- catalog_by_id[[candidate_id]]
  strategy <- as.character(or_value(candidate$writeback_strategy, ""))
  if (!identical(strategy, "canonical_cluster_membership")) return(FALSE)
  anti_compatible <- row$anti_marker_deg_log2fc_mean <= 0 ||
    row$marker_deg_log2fc_mean - row$anti_marker_deg_log2fc_mean >= 1.00
  row$available_positive_family_count >= 2 &&
    row$group_positive_family_supported_count >= 2 &&
    row$observation_identity_core_fraction >= max(
      0.10, minimum_identity_core_fraction(candidate_id)
    ) &&
    row$observation_identity_core_direct_fraction >= 0.10 &&
    row$observation_release_family_coherent_fraction >= 0.02 &&
    row$positive_marker_detection_fraction >= 0.25 &&
    row$marker_deg_log2fc_mean >= 1.50 && anti_compatible
}

aggregate_program_supported <- function(row, candidate_id, family_evidence) {
  if (canonical_cluster_challenger_supported(row, candidate_id)) {
    return(TRUE)
  }
  if (
    is.null(row) || !nrow(row) || !family_evidence$pass ||
      row$positive_marker_detection_fraction < 0.05 ||
      row$mean_program_score < 0.02
  ) {
    return(FALSE)
  }
  strategy <- as.character(or_value(
    catalog_by_id[[candidate_id]]$writeback_strategy, ""
  ))
  canonical <- identical(strategy, "canonical_cluster_membership") &&
    row$marker_deg_log2fc_mean >= 1.50 &&
    row$observation_coherent_fraction >= 0.25
  common <- row$marker_deg_log2fc_mean >= 0.50 &&
    row$observation_coherent_fraction >= 0.05
  seeded <- row$marker_deg_log2fc_mean >= 0.25 &&
    row$observation_seed_fraction >= 0.05 &&
    row$observation_coherent_fraction >= 0.10
  anti_compatible <- row$anti_marker_deg_log2fc_mean <= 0 ||
    row$marker_deg_log2fc_mean -
      row$anti_marker_deg_log2fc_mean >= 0.50
  (canonical || common || seeded) && anti_compatible
}

component_enrichment <- function(members, target_id, group_scores) {
  target_all <- group_scores[candidate_id == target_id]
  member_rows <- target_all[cell_id %chin% members]
  background <- target_all[!cell_id %chin% members]
  if (!nrow(member_rows)) {
    return(list(
      pass = FALSE, program_delta = -Inf, direct_delta = -Inf,
      member_mean_program = -Inf, member_mean_direct = -Inf
    ))
  }
  member_program <- mean(member_rows$program_score)
  member_direct <- mean(member_rows$direct_signal)
  if (nrow(background) < 20L) {
    return(list(
      pass = member_program >= 0.04 || member_direct >= 0.08,
      program_delta = 0, direct_delta = 0,
      member_mean_program = member_program,
      member_mean_direct = member_direct
    ))
  }
  program_delta <- member_program - mean(background$program_score)
  direct_delta <- member_direct - mean(background$direct_signal)
  list(
    pass = member_program >= 0.02 &&
      (program_delta >= 0.02 || direct_delta >= 0.03),
    program_delta = program_delta,
    direct_delta = direct_delta,
    member_mean_program = member_program,
    member_mean_direct = member_direct
  )
}

cell_table <- unique(scores[, .(
  cell_id, source_boundary, source_cluster, x, y
)])
setorder(cell_table, cell_id)
coordinates <- as.matrix(cell_table[, .(x, y)])
nn_k <- min(12L, nrow(cell_table))
nn <- RANN::nn2(coordinates, coordinates, k = nn_k)
positive_distances <- nn$nn.dists[, 2L][nn$nn.dists[, 2L] > 0]
spatial_radius <- 3 * median(positive_distances, na.rm = TRUE)
if (!is.finite(spatial_radius) || spatial_radius <= 0) {
  stop("cannot derive spatial component radius")
}
cell_index <- setNames(seq_len(nrow(cell_table)), cell_table$cell_id)

subset_rows <- list()
membership_rows <- list()
watch_rows <- list()
subset_i <- 1L
membership_i <- 1L
watch_i <- 1L

validate_group <- function(
  members, target_id, group_scores, aggregate_row = NULL,
  minimum_support = subset_support_threshold,
  require_component_enrichment = TRUE
) {
  target <- group_scores[candidate_id == target_id & cell_id %chin% members]
  target_broad <- candidate_broad[[target_id]]
  target_role <- candidate_role[[target_id]]
  target_family <- group_family_evidence(target, target_id)
  enrichment <- component_enrichment(members, target_id, group_scores)
  parent_candidate_id <- ""
  parent_family_status <- "NOT_REQUIRED"
  parent_supported_families <- ""
  parent_support <- 1
  parent_contradiction <- 0
  parent_enrichment_status <- "NOT_REQUIRED"
  parent_identity_status <- "NOT_REQUIRED"
  if (identical(target_role, "fine")) {
    proposed_parent <- as.character(or_value(
      catalog_by_id[[target_id]]$context_evidence_candidate_id, ""
    ))
    parent_candidate <- catalog_by_id[[proposed_parent]]
    if (
      nzchar(proposed_parent) && !is.null(parent_candidate) &&
        identical(
          as.character(or_value(parent_candidate$candidate_role, "")),
          "broad"
        ) &&
        identical(
          as.character(or_value(parent_candidate$release_broad_label, "")),
          target_broad
        )
    ) {
      parent_candidate_id <- proposed_parent
      parent_target <- group_scores[
        candidate_id == proposed_parent & cell_id %chin% members
      ]
      parent_family <- group_family_evidence(
        parent_target, proposed_parent
      )
      parent_enrichment <- component_enrichment(
        members, proposed_parent, group_scores
      )
      parent_family_status <- if (parent_family$pass) "PASS" else "FAIL"
      parent_supported_families <- parent_family$supported_families
      parent_support <- mean(parent_target$anchor)
      parent_contradiction <- mean(parent_target$hard_block)
      parent_enrichment_status <- if (
        parent_enrichment$pass
      ) "PASS" else "FAIL"
      parent_required_families <- as.character(unlist(
        parent_candidate$required_positive_families
      ))
      parent_family_floor_pass <- if (length(parent_required_families)) {
        parent_family$pass
      } else {
        nzchar(parent_supported_families)
      }
      parent_identity_status <- if (
        parent_family_floor_pass && parent_support >= parent_support_threshold &&
          parent_contradiction <= maximum_contradiction_threshold
      ) "PASS" else "FAIL"
    }
  }
  aggregate_override <- !is.null(aggregate_row) && nrow(aggregate_row)
  support <- if (aggregate_override) {
    strategy <- as.character(or_value(
      catalog_by_id[[target_id]]$writeback_strategy, ""
    ))
    if (identical(strategy, "canonical_cluster_membership")) {
      aggregate_row$observation_coherent_fraction
    } else {
      aggregate_row$observation_identity_core_fraction
    }
  } else {
    mean(target$anchor)
  }
  contradiction <- if (
    aggregate_override &&
      aggregate_row$marker_deg_log2fc_mean >= 1.50 &&
      aggregate_row$anti_marker_deg_log2fc_mean <= 0
  ) {
    0
  } else if (aggregate_override) {
    aggregate_row$hard_contradiction_fraction
  } else {
    mean(target$hard_block)
  }
  competitor_rows <- group_scores[
    cell_id %chin% members & effective_broad_label != target_broad
  ]
  if (release_level != "broad" && target_role == "fine") {
    competitor_rows <- rbind(
      competitor_rows,
      group_scores[
        cell_id %chin% members &
          effective_broad_label == target_broad &
          effective_candidate_role == "fine" &
          candidate_id != target_id
      ],
      use.names = TRUE
    )
  }
  target_normalized <- setNames(target$normalized_evidence, target$cell_id)
  competitors <- competitor_rows[candidate_id != target_id, {
    evidence <- group_family_evidence(.SD, candidate_id[[1L]])
    target_value <- unname(target_normalized[cell_id])
    .(
      fraction = if (evidence$pass) {
        mean(
          anchor & is.finite(target_value) &
            normalized_evidence >= target_value + 0.05
        )
      } else {
        0
      },
      family_support_status = if (evidence$pass) "PASS" else "FAIL"
    )
  }, by = candidate_id][order(-fraction, candidate_id)]
  competitor_fraction <- if (nrow(competitors)) competitors$fraction[[1L]] else 0
  competitor_id <- if (nrow(competitors)) competitors$candidate_id[[1L]] else ""
  list(
    support = support,
    contradiction = contradiction,
    competitor_fraction = competitor_fraction,
    competitor_id = competitor_id,
    margin = support - competitor_fraction,
    family_support_status = if (target_family$pass) "PASS" else "FAIL",
    supported_families = target_family$supported_families,
    family_prevalence = target_family$family_prevalence,
    component_enrichment_status = if (enrichment$pass) "PASS" else "FAIL",
    parent_candidate_id = parent_candidate_id,
    parent_identity_status = parent_identity_status,
    parent_lineage_supported_fraction = parent_support,
    parent_contradiction_fraction = parent_contradiction,
    parent_family_support_status = parent_family_status,
    parent_supported_families = parent_supported_families,
    parent_component_enrichment_status = parent_enrichment_status,
    program_score_delta = enrichment$program_delta,
    direct_signal_delta = enrichment$direct_delta,
    pass = support >= minimum_support &&
      support - competitor_fraction >= subset_margin_threshold &&
      contradiction <= maximum_contradiction_threshold && target_family$pass &&
      parent_identity_status != "FAIL" &&
      (!require_component_enrichment || enrichment$pass)
  )
}

for (boundary in unique(scores$source_boundary)) {
  boundary_scores <- scores[source_boundary == boundary]
  for (cluster in unique(boundary_scores$source_cluster)) {
    group_scores <- boundary_scores[source_cluster == cluster]
    group_cells <- unique(group_scores$cell_id)
    group_n <- length(group_cells)
    group_summary <- group_scores[, {
      family_evidence <- group_family_evidence(.SD, candidate_id[[1L]])
      .(
        seed_fraction = mean(anchor),
        contradiction_fraction = mean(hard_block),
        mean_program_score = mean(program_score),
        cross_resolution_fraction = mean(cross_resolution_support_count >= 2),
        spatial_support_fraction = mean(local_seed_fraction >= 0.03),
        family_support = family_evidence$pass,
        supported_families = family_evidence$supported_families
      )
    }, by = candidate_id][order(-seed_fraction, -mean_program_score, candidate_id)]

    # Only a single, aggregate-supported specific broad program may return a
    # selected-resolution cluster wholesale. Generic Stromal is deliberately
    # deferred until all embedded specific programs and exact remainders close.
    selected_aggregate <- list()
    selected_i <- 1L
    for (candidate_id in candidate_ids) {
      candidate_value <- candidate_id
      candidate_rows <- group_scores[candidate_id == candidate_value]
      family_evidence <- group_family_evidence(candidate_rows, candidate_id)
      aggregate_row <- aggregate_row_for(
        "selected", boundary, cluster, candidate_id
      )
      identity_core_fraction <- mean(identity_core_mask(
        candidate_rows, candidate_id
      ))
      if (aggregate_program_supported(
        aggregate_row, candidate_id, family_evidence
      ) && identity_core_fraction >= minimum_identity_core_fraction(
        candidate_id
      )) {
        selected_aggregate[[selected_i]] <- data.table(
          candidate_id = candidate_id,
          candidate_role = candidate_role[[candidate_id]],
          broad_label = candidate_broad[[candidate_id]],
          identity_token = ifelse(
            nzchar(candidate_broad[[candidate_id]]),
            candidate_broad[[candidate_id]],
            paste0("exploratory::", candidate_id)
          ),
          aggregate_score = aggregate_score(aggregate_row),
          identity_core_fraction = identity_core_fraction
        )
        selected_i <- selected_i + 1L
      }
    }
    selected_aggregate <- if (length(selected_aggregate)) {
      rbindlist(selected_aggregate)
    } else {
      data.table()
    }
    whole_candidate <- ""
    if (nrow(selected_aggregate)) {
      supported_identities <- unique(selected_aggregate$identity_token)
      releasable_broad <- selected_aggregate[
        candidate_id %chin% candidate_ids[release_eligible] &
          candidate_role == "broad"
      ]
      if (length(supported_identities) == 1L && nrow(releasable_broad)) {
        releasable_broad <- releasable_broad[
          order(-aggregate_score, candidate_id)
        ]
        proposed <- releasable_broad$candidate_id[[1L]]
        strategy <- as.character(or_value(
          effective_broad_writeback_strategy(proposed), ""
        ))
        if (!identical(
          strategy, "generic_exact_remainder_after_specific_lineages"
        )) {
          whole_candidate <- proposed
        }
      }
    }
    if (nzchar(whole_candidate)) {
      target <- group_scores[candidate_id == whole_candidate]
      identity_core_fraction <- mean(identity_core_mask(
        target, whole_candidate
      ))
      aggregate_row <- aggregate_row_for(
        "selected", boundary, cluster, whole_candidate
      )
      canonical_strong <- aggregate_row$marker_deg_log2fc_mean >= 1.50 &&
        aggregate_row$anti_marker_deg_log2fc_mean <= 0
      members <- target[
        !truth(technical_flag) & (canonical_strong | !hard_block),
        cell_id
      ]
      evidence <- validate_group(
        members, whole_candidate, group_scores,
        aggregate_row = aggregate_row,
        minimum_support = subset_support_threshold,
        require_component_enrichment = FALSE
      )
      subset_id <- paste(
        safe(boundary), safe(cluster), whole_candidate, "whole_subcluster",
        sep = "__"
      )
      subset_rows[[subset_i]] <- data.table(
        subset_id = subset_id,
        source_boundary = boundary,
        source_cluster = cluster,
        candidate_id = whole_candidate,
        proposal_scope = "whole_subcluster",
        n_observations = length(members),
        lineage_supported_fraction = evidence$support,
        strongest_competing_candidate = evidence$competitor_id,
        strongest_competing_fraction = evidence$competitor_fraction,
        support_margin = evidence$margin,
        contradiction_fraction = evidence$contradiction,
        family_support_status = evidence$family_support_status,
        supported_families = evidence$supported_families,
        family_prevalence = evidence$family_prevalence,
        component_enrichment_status = evidence$component_enrichment_status,
        program_score_delta = evidence$program_score_delta,
        direct_signal_delta = evidence$direct_signal_delta,
        cross_resolution_supported_fraction = mean(
          target[cell_id %chin% members]$cross_resolution_support_count >= 2
        ),
        spatially_supported_fraction = mean(
          target[cell_id %chin% members]$local_seed_fraction >= 0.03
        ),
        identity_core_fraction = identity_core_fraction,
        component_geodesic_policy =
          "whole_expression_subcluster_sparse_tail_inheritance",
        maximum_expansion_hops = 0L,
        status = ifelse(evidence$pass, "PASS", "FAIL")
      )
      membership_rows[[membership_i]] <- data.table(
        subset_id = subset_id, cell_id = sort(members)
      )
      subset_i <- subset_i + 1L
      membership_i <- membership_i + 1L
      if (evidence$pass) next
    }

    # Every candidate proposes independently. Multiple candidates may emerge
    # from one mixed subcluster; overlap is resolved later from normalized
    # evidence and pairwise discriminators.
    build_candidate_result <- function(candidate_id) {
      data.table::setDTthreads(1L)
      target_id <- candidate_id
      candidate_scores <- group_scores[candidate_id == target_id]
      candidate <- catalog_by_id[[candidate_id]]
      if (!release_eligible[[candidate_id]]) {
        if (mean(candidate_scores$seed) >= 0.03) {
          return(list(
            subsets = list(),
            memberships = list(),
            watches = list(data.table(
              source_boundary = boundary,
              source_cluster = cluster,
              candidate_id = candidate_id,
              candidate_role = as.character(candidate$candidate_role),
              seed_fraction = mean(candidate_scores$seed),
              mean_program_score = mean(candidate_scores$program_score),
              status = "watch"
            ))
          ))
        }
        return(list(subsets = list(), memberships = list(), watches = list()))
      }
      if (identical(
        as.character(or_value(candidate$writeback_strategy, "")),
        "generic_exact_remainder_after_specific_lineages"
      )) {
        return(list(subsets = list(), memberships = list(), watches = list()))
      }
      local_subsets <- list()
      local_memberships <- list()
      effective_strategy <- effective_broad_writeback_strategy(candidate_id)

      # A neighboring-resolution expression subcluster is an independent
      # proposal source. Strong aggregate DEG/pseudobulk evidence may return its sparse
      # noncontradictory tail even when few individual cellbins carry two
      # marker families.
      for (resolution_role in if (
        identical(
          effective_strategy,
          "candidate_local_component_never_parent_expansion"
        )
      ) character() else c("neighbor_1", "neighbor_2")) {
        boundary_column <- paste0(resolution_role, "_boundary")
        cluster_column <- paste0(resolution_role, "_cluster")
        subclusters <- unique(candidate_scores[, .(
          neighbor_boundary = get(boundary_column),
          neighbor_cluster = get(cluster_column)
        )])
        for (subcluster_i in seq_len(nrow(subclusters))) {
          neighbor_boundary <- subclusters$neighbor_boundary[[subcluster_i]]
          neighbor_cluster <- subclusters$neighbor_cluster[[subcluster_i]]
          neighbor_target <- candidate_scores[
            get(boundary_column) == neighbor_boundary &
              get(cluster_column) == neighbor_cluster
          ]
          if (nrow(neighbor_target) < 5L) next
          neighbor_members <- neighbor_target$cell_id
          family_evidence <- group_family_evidence(
            neighbor_target, candidate_id
          )
          target_aggregate <- aggregate_row_for(
            resolution_role, neighbor_boundary, neighbor_cluster, candidate_id
          )
          identity_core_fraction <- mean(identity_core_mask(
            neighbor_target, candidate_id
          ))
          if (!aggregate_program_supported(
            target_aggregate, candidate_id, family_evidence
          ) || identity_core_fraction < minimum_identity_core_fraction(
            candidate_id
          )) {
            next
          }
          target_aggregate_score <- aggregate_score(target_aggregate)
          aggregate_competitors <- list()
          aggregate_competitor_i <- 1L
          for (other_id in candidate_ids[candidate_ids != candidate_id]) {
            other_broad <- candidate_broad[[other_id]]
            other_role <- candidate_role[[other_id]]
            same_broad <- identical(other_broad, candidate_broad[[candidate_id]])
            if (
              same_broad &&
                (
                  release_level == "broad" ||
                    candidate_role[[candidate_id]] != "fine" ||
                    other_role != "fine"
                )
            ) {
              next
            }
            other_rows <- group_scores[
              candidate_id == other_id & cell_id %chin% neighbor_members
            ]
            other_family <- group_family_evidence(other_rows, other_id)
            other_aggregate <- aggregate_row_for(
              resolution_role, neighbor_boundary, neighbor_cluster, other_id
            )
            if (aggregate_program_supported(
              other_aggregate, other_id, other_family
            )) {
              aggregate_competitors[[aggregate_competitor_i]] <- data.table(
                candidate_id = other_id,
                aggregate_score = aggregate_score(other_aggregate)
              )
              aggregate_competitor_i <- aggregate_competitor_i + 1L
            }
          }
          aggregate_competitors <- if (length(aggregate_competitors)) {
            rbindlist(aggregate_competitors)[
              order(-aggregate_score, candidate_id)
            ]
          } else {
            data.table()
          }
          competitor_id <- if (nrow(aggregate_competitors)) {
            aggregate_competitors$candidate_id[[1L]]
          } else {
            ""
          }
          competitor_score <- if (nrow(aggregate_competitors)) {
            aggregate_competitors$aggregate_score[[1L]]
          } else {
            -Inf
          }
          aggregate_margin <- target_aggregate_score - competitor_score
          strong_cluster_program <-
            target_aggregate$marker_deg_log2fc_mean >= 1.50 &&
            target_aggregate$anti_marker_deg_log2fc_mean <= 0
          members <- neighbor_target[
            !truth(technical_flag) &
              (strong_cluster_program | !hard_block),
            cell_id
          ]
          evidence <- validate_group(
            members, candidate_id, group_scores,
            aggregate_row = target_aggregate,
            minimum_support = subset_support_threshold,
            require_component_enrichment = FALSE
          )
          aggregate_pass <- evidence$family_support_status == "PASS" &&
            evidence$support >= subset_support_threshold &&
            evidence$contradiction <= maximum_contradiction_threshold &&
            aggregate_margin >= 0.25
          subset_id <- paste(
            safe(boundary), safe(cluster), candidate_id,
            paste0(
              resolution_role, "_", safe(neighbor_boundary), "_",
              safe(neighbor_cluster)
            ),
            sep = "__"
          )
          local_subsets[[length(local_subsets) + 1L]] <- data.table(
            subset_id = subset_id,
            source_boundary = boundary,
            source_cluster = cluster,
            candidate_id = candidate_id,
            proposal_scope = "neighboring_resolution_expression_subcluster",
            n_observations = length(members),
            lineage_supported_fraction = evidence$support,
            strongest_competing_candidate = competitor_id,
            strongest_competing_fraction =
              evidence$competitor_fraction,
            support_margin = evidence$margin,
            contradiction_fraction = evidence$contradiction,
            family_support_status = evidence$family_support_status,
            supported_families = evidence$supported_families,
            family_prevalence = evidence$family_prevalence,
            component_enrichment_status = "aggregate_multichannel",
            program_score_delta = evidence$program_score_delta,
            direct_signal_delta = evidence$direct_signal_delta,
            aggregate_score = target_aggregate_score,
            aggregate_competitor_score = competitor_score,
            aggregate_margin = aggregate_margin,
            source_resolution_role = resolution_role,
            source_resolution_cluster = neighbor_cluster,
            cross_resolution_supported_fraction =
              target_aggregate$cross_resolution_stable_fraction,
            spatially_supported_fraction =
              target_aggregate$spatial_local_support_fraction,
            identity_core_fraction = identity_core_fraction,
            component_geodesic_policy =
              "expression_subcluster_sparse_tail_inheritance",
            maximum_expansion_hops = 0L,
            status = ifelse(aggregate_pass, "PASS", "FAIL")
          )
          local_memberships[[length(local_memberships) + 1L]] <- data.table(
            subset_id = subset_id, cell_id = sort(members)
          )
        }
      }
      # Components are made only from candidate identity cores.  In v2.2.0
      # generic parent/support anchors must never form a transitive bridge
      # between distant identity-positive foci.  Dropout tails are handled by
      # validated expression subclusters above, not by spatial propagation.
      # A directly contradicted observation cannot seed or connect a local
      # lineage component.  Keeping it in the seed graph made a coherent
      # Smooth-muscle/Granulosa core fail only because contradictory members
      # inflated the component-wide hard-anti fraction.  Those observations
      # remain in the exact remainder and are reconsidered for a coarse parent;
      # the final subset thresholds themselves stay unchanged.
      core_mask <- identity_core_mask(candidate_scores, candidate_id) &
        !candidate_scores$hard_block
      seeds <- candidate_scores[core_mask == TRUE, cell_id]
      if (!length(seeds)) {
        return(list(
          subsets = local_subsets,
          memberships = local_memberships,
          watches = list()
        ))
      }
      global_rows <- unname(cell_index[seeds])
      core_coordinates <- coordinates[global_rows, , drop = FALSE]
      core_k <- min(12L, nrow(core_coordinates))
      core_nn <- RANN::nn2(
        core_coordinates, core_coordinates, k = core_k
      )
      neighbor_matrix <- core_nn$nn.idx[, -1L, drop = FALSE]
      distance_matrix <- core_nn$nn.dists[, -1L, drop = FALSE]
      edge_from_all <- rep(
        seq_len(nrow(core_coordinates)), each = ncol(neighbor_matrix)
      )
      edge_to_all <- as.vector(t(neighbor_matrix))
      neighbor_distance <- as.vector(t(distance_matrix))
      keep <- neighbor_distance <= 2 * spatial_radius
      edge_from <- edge_from_all[keep]
      edge_to <- edge_to_all[keep]
      graph <- make_empty_graph(n = length(global_rows), directed = FALSE)
      if (length(edge_from)) {
        graph <- add_edges(graph, as.vector(rbind(edge_from, edge_to)))
      }
      components <- igraph::components(graph)$membership
      for (component in sort(unique(components))) {
        members <- seeds[components == component]
        if (length(members) < minimum_component_n) next
        evidence <- validate_group(members, candidate_id, group_scores)
        target <- candidate_scores[cell_id %chin% members]
        selected_aggregate_row <- aggregate_row_for(
          "selected", boundary, cluster, candidate_id
        )
        spatial_fraction <- mean(target$local_seed_fraction >= 0.03)
        canonical_component_pass <-
          canonical_cluster_challenger_supported(
            selected_aggregate_row, candidate_id
          ) &&
          evidence$family_support_status == "PASS" &&
          evidence$component_enrichment_status == "PASS" &&
          evidence$program_score_delta >= 0.05 &&
          evidence$direct_signal_delta >= 0.05 &&
          spatial_fraction >= 0.40
        ordinary_component_pass <- evidence$pass &&
          (
            mean(target$cross_resolution_support_count >= 2) >= 0.25 ||
              evidence$program_score_delta >= 0.05 ||
              evidence$direct_signal_delta >= 0.05
          )
        subset_id <- paste(
          safe(boundary), safe(cluster), candidate_id,
          paste0("component", component), sep = "__"
        )
        local_subsets[[length(local_subsets) + 1L]] <- data.table(
          subset_id = subset_id,
          source_boundary = boundary,
          source_cluster = cluster,
          candidate_id = candidate_id,
          proposal_scope = ifelse(
            canonical_component_pass,
            "canonical_identity_component",
            "candidate_local_spatial_component"
          ),
          n_observations = length(members),
          lineage_supported_fraction = evidence$support,
          strongest_competing_candidate = evidence$competitor_id,
          strongest_competing_fraction = evidence$competitor_fraction,
          support_margin = evidence$margin,
          contradiction_fraction = evidence$contradiction,
          family_support_status = evidence$family_support_status,
          supported_families = evidence$supported_families,
          family_prevalence = evidence$family_prevalence,
          component_enrichment_status = evidence$component_enrichment_status,
          program_score_delta = evidence$program_score_delta,
          direct_signal_delta = evidence$direct_signal_delta,
          aggregate_score = NA_real_,
          aggregate_competitor_score = NA_real_,
          aggregate_margin = NA_real_,
          source_resolution_role = "",
          source_resolution_cluster = "",
          cross_resolution_supported_fraction = mean(
            target$cross_resolution_support_count >= 2
          ),
          spatially_supported_fraction = spatial_fraction,
          identity_core_fraction = 1,
          component_geodesic_policy =
            "identity_core_only_no_support_cell_bridges",
          maximum_expansion_hops = 0L,
          status = ifelse(
            canonical_component_pass || ordinary_component_pass,
            "PASS", "FAIL"
          )
        )
        local_memberships[[length(local_memberships) + 1L]] <- data.table(
          subset_id = subset_id, cell_id = sort(members)
        )
      }
      list(
        subsets = local_subsets,
        memberships = local_memberships,
        watches = list()
      )
    }
    candidate_workers <- min(workers, length(candidate_ids))
    candidate_results <- if (candidate_workers > 1L) {
      parallel::mclapply(
        candidate_ids,
        build_candidate_result,
        mc.cores = candidate_workers,
        mc.preschedule = TRUE
      )
    } else {
      lapply(candidate_ids, build_candidate_result)
    }
    setDTthreads(max(1L, workers))
    for (result in candidate_results) {
      for (row in result$subsets) {
        subset_rows[[subset_i]] <- row
        subset_i <- subset_i + 1L
      }
      for (row in result$memberships) {
        membership_rows[[membership_i]] <- row
        membership_i <- membership_i + 1L
      }
      for (row in result$watches) {
        watch_rows[[watch_i]] <- row
        watch_i <- watch_i + 1L
      }
    }
  }
}

subset_table <- if (length(subset_rows)) rbindlist(subset_rows, fill = TRUE) else data.table()
membership_table <- if (length(membership_rows)) rbindlist(membership_rows, fill = TRUE) else data.table()
watch_table <- if (length(watch_rows)) rbindlist(watch_rows, fill = TRUE) else data.table()
if (!nrow(subset_table)) {
  subset_table <- data.table(
    subset_id = character(), source_boundary = character(),
    source_cluster = character(), candidate_id = character(),
    proposal_scope = character(), n_observations = integer(),
    lineage_supported_fraction = numeric(),
    strongest_competing_candidate = character(),
    strongest_competing_fraction = numeric(), support_margin = numeric(),
    contradiction_fraction = numeric(),
    family_support_status = character(), supported_families = character(),
    family_prevalence = character(),
    component_enrichment_status = character(),
    program_score_delta = numeric(), direct_signal_delta = numeric(),
    aggregate_score = numeric(), aggregate_competitor_score = numeric(),
    aggregate_margin = numeric(), source_resolution_role = character(),
    source_resolution_cluster = character(),
    cross_resolution_supported_fraction = numeric(),
    spatially_supported_fraction = numeric(),
    identity_core_fraction = numeric(),
    component_geodesic_policy = character(),
    maximum_expansion_hops = integer(), status = character()
  )
}
if (!nrow(membership_table)) {
  membership_table <- data.table(subset_id = character(), cell_id = character())
}
if (!nrow(watch_table)) {
  watch_table <- data.table(
    source_boundary = character(), source_cluster = character(),
    candidate_id = character(), candidate_role = character(),
    seed_fraction = numeric(), mean_program_score = numeric(),
    status = character()
  )
}
setorder(subset_table, source_boundary, source_cluster, candidate_id, subset_id)
setorder(membership_table, subset_id, cell_id)
fwrite(subset_table, file.path(out_dir, "candidate_subset_evidence.tsv"), sep = "\t")
fwrite(
  membership_table,
  file.path(out_dir, "candidate_subset_membership.tsv.gz"),
  sep = "\t"
)
fwrite(watch_table, file.path(out_dir, "candidate_watch_ledger.tsv"), sep = "\t")

manifest <- list(
  status = "PASS",
  schema_version = "2.2",
  controller_version = "2.2.0",
  stage = "candidate_local_subset_derivation",
  activation_scope = "second_round_mixed_subcluster_only",
  source_boundary = source_boundary_filter,
  source_cluster = source_cluster_filter,
  scores = score_path,
  cluster_evidence = cluster_evidence_path,
  catalog = catalog_path,
  threshold_registry = list(
    path = threshold_registry_path,
    sha256 = sha256(threshold_registry_path)
  ),
  release_level = release_level,
  requested_workers = workers,
  effective_candidate_workers = min(workers, length(candidate_ids)),
  parallel_backend = "parallel_mclapply_candidate_components",
  parallel_unit = "candidate_within_source_cluster",
  n_subsets_proposed = nrow(subset_table),
  n_subsets_preliminarily_passing = sum(subset_table$status == "PASS"),
  n_watch_programs = nrow(watch_table),
  spatial_radius = spatial_radius,
  policy = list(
    aggregate_winner_can_veto_candidate = FALSE,
    candidates_proposed_independently = TRUE,
    overlap_assignment_by_catalog_order = FALSE,
    generic_support_cells_can_bridge_identity_cores = FALSE,
    spatial_component_membership =
      "declared identity-core observations only",
    group_validation = paste0(
      "support>=", subset_support_threshold,
      ";margin>=", subset_margin_threshold,
      ";contradiction<=", maximum_contradiction_threshold
    ),
    validation_threshold_used_as_cell_admission = FALSE,
    high_purity_sparse_tail_inherits_broad = TRUE
  )
)
write_json(
  manifest, file.path(out_dir, "candidate_subset_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE
)
writeLines("status\tPASS", file.path(out_dir, "RUN_COMPLETE.tsv"))
