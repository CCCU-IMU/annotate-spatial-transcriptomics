#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(data.table)
  library(jsonlite); library(RANN); library(ggplot2)
})

parse_args <- function(x) {
  o <- list(); i <- 1L
  while (i <= length(x)) {
    k <- sub("^--", "", x[i])
    if (i == length(x) || startsWith(x[i + 1L], "--")) { o[[k]] <- TRUE; i <- i + 1L }
    else { o[[k]] <- x[i + 1L]; i <- i + 2L }
  }
  o
}
read_any <- function(p) if (grepl("\\.gz$", p)) fread(cmd = paste("gzip -dc", shQuote(p))) else fread(p)
canon <- function(x) {
  x <- as.character(x); out <- x
  out[grepl("stromal", x, ignore.case = TRUE)] <- "Stromal"
  out[grepl("immune|plasma|macroph|neutroph", x, ignore.case = TRUE)] <- "Immune/plasma"
  out[grepl("epithelial|mesothelial", x, ignore.case = TRUE)] <- "Epithelial/mesothelial"
  out[grepl("granulosa|follicular somatic", x, ignore.case = TRUE)] <- "Granulosa"
  out[grepl("theca", x, ignore.case = TRUE)] <- "Theca"
  out[grepl("vascular|endothelial|mural|perivascular", x, ignore.case = TRUE)] <- "Vascular/perivascular"
  out[grepl("oocyte", x, ignore.case = TRUE)] <- "Oocyte"
  out[grepl("stromal", x, ignore.case = TRUE) & grepl("vascular|perivascular", x, ignore.case = TRUE)] <- "Stromal"
  out
}

a <- parse_args(commandArgs(trailingOnly = TRUE))
need <- c("rds", "calibrated-mapping", "clusters", "cell-ledger", "out", "cell-id-col", "cluster-col")
miss <- need[!need %in% names(a)]; if (length(miss)) stop("Missing: ", paste(miss, collapse = ", "))
cc <- a$`cell-id-col`; gc <- a$`cluster-col`; outdir <- a$out
dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "provenance"), recursive = TRUE, showWarnings = FALSE)

min_program_hits <- as.integer(ifelse(is.null(a$`min-program-hits`), 2, a$`min-program-hits`))
min_cluster_any <- as.numeric(ifelse(is.null(a$`min-cluster-any-fraction`), 0.15, a$`min-cluster-any-fraction`))
min_cluster_strong <- as.numeric(ifelse(is.null(a$`min-cluster-strong-fraction`), 0.03, a$`min-cluster-strong-fraction`))
spatial_fraction <- as.numeric(ifelse(is.null(a$`min-spatial-fraction`), 0.35, a$`min-spatial-fraction`))
radius_multiplier <- as.numeric(ifelse(is.null(a$`radius-multiplier`), 10, a$`radius-multiplier`))
k <- as.integer(ifelse(is.null(a$`spatial-neighbors`), 20, a$`spatial-neighbors`))
seed <- as.integer(ifelse(is.null(a$seed), 20260713, a$seed)); set.seed(seed)

obj <- readRDS(a$rds); mp <- read_any(a$`calibrated-mapping`); cl <- read_any(a$clusters); led <- read_any(a$`cell-ledger`)
required_mapping_cols <- c("mapping_tier", "meets_moderate_or_higher", "fine_anchor_eligible")
if (!all(required_mapping_cols %in% names(mp))) {
  stop("Legacy/combined mapping tiers are forbidden. Use calibrate_tiered_mapping_thresholds.py and provide high/moderate/low tiers.")
}
if (any(tolower(as.character(mp$fine_anchor_eligible)) %in% c("true", "1", "yes"))) stop("Atlas mapping cannot create fine anchors")
if (any(as.character(mp$mapping_tier) %in% c("high", "moderate") & !tolower(as.character(mp$meets_moderate_or_higher)) %in% c("true", "1", "yes"))) stop("accepted tier violates moderate-or-higher invariant")
mp[[cc]] <- as.character(mp[[cc]]); cl[[cc]] <- as.character(cl[[cc]])
if (uniqueN(mp[[cc]]) != nrow(mp) || uniqueN(cl[[cc]]) != nrow(cl)) stop("duplicate query IDs")
if (!all(mp[[cc]] %in% cl[[cc]]) || !all(mp[[cc]] %in% colnames(obj))) stop("query boundary mismatch")
mp <- merge(mp, cl[, .(cell_id = get(cc), cluster = as.character(get(gc)))], by.x = cc, by.y = "cell_id", all.x = TRUE, sort = FALSE)
setnames(mp, cc, "cell_id"); mp[, predicted_label_canonical := canon(predicted_label)]

programs <- list(
  Stromal = c("DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "COL6A1", "PDGFRA", "TCF21", "NR2F2", "OGN", "SPARC", "MGP", "CXCL14", "IGFBP5"),
  Granulosa = c("FOXL2", "FSHR", "FST", "AMH", "CDH2", "INHA", "INHBA", "INHBB", "GJA1", "HSD17B1", "CYP19", "CYP19A1", "SERPINE2"),
  Theca = c("LHCGR", "CYP17A1", "STAR", "CYP11A1", "HSD3B1", "INSL3", "NR5A1"),
  `Vascular/perivascular` = c("PECAM1", "CDH5", "VWF", "CLDN5", "KDR", "FLT1", "EMCN", "RGS5", "PDGFRB", "CSPG4"),
  `Immune/plasma` = c("PTPRC", "CD53", "LYZ", "TYROBP", "FCER1G", "GPNMB", "ACP5", "S100A8", "S100A9", "PGLYRP1", "JCHAIN", "CD74"),
  `Epithelial/mesothelial` = c("KRT8", "KRT18", "KRT19", "KRT7", "EPCAM", "DSP", "MUC16", "MSLN", "PRG4"),
  Oocyte = c("DDX4", "DAZL", "FIGLA", "NOBOX", "PADI6", "OOEP", "NPM2", "ZP2", "ZP3", "ZP4")
)
cts <- tryCatch(LayerData(obj[["Spatial"]], layer = "counts"), error = function(e) GetAssayData(obj[["Spatial"]], slot = "counts"))
ids <- mp$cell_id; hit_names <- character()
for (nm in names(programs)) {
  g <- intersect(programs[[nm]], rownames(cts)); hn <- paste0(gsub("[^A-Za-z0-9]+", "_", nm), "_hits")
  mp[[hn]] <- if (length(g)) as.integer(Matrix::colSums(cts[g, ids, drop = FALSE] > 0)) else 0L
  hit_names[nm] <- hn
}
mp[, predicted_program_hits := 0L]
for (nm in names(hit_names)) mp[predicted_label_canonical == nm, predicted_program_hits := get(hit_names[[nm]])]
anti_cols <- unname(hit_names)
mp[, max_competing_hits := 0L]
for (i in seq_len(nrow(mp))) {
  own <- hit_names[[mp$predicted_label_canonical[i]]]
  cols <- setdiff(anti_cols, if (is.null(own)) character() else own)
  mp$max_competing_hits[i] <- if (length(cols)) max(as.integer(mp[i, ..cols])) else 0L
}
mp[, fullfeature_program_consistent := predicted_program_hits >= min_program_hits & predicted_program_hits >= max_competing_hits]

context <- rbindlist(lapply(names(hit_names), function(nm) {
  hn <- hit_names[[nm]]
  mp[, .(cluster_n = .N, cluster_program_any_fraction = mean(get(hn) >= 1L),
         cluster_program_strong_fraction = mean(get(hn) >= min_program_hits), cluster_program_mean_hits = mean(get(hn))), by = cluster][, predicted_label_canonical := nm]
}))
mp <- merge(mp, context, by = c("cluster", "predicted_label_canonical"), all.x = TRUE, sort = FALSE)
mp[, cluster_context_consistent := cluster_program_any_fraction >= min_cluster_any & cluster_program_strong_fraction >= min_cluster_strong]

md <- obj@meta.data; xy_cols <- if (all(c("x", "y") %in% names(md))) c("x", "y") else if (all(c("sdimx", "sdimy") %in% names(md))) c("sdimx", "sdimy") else stop("spatial coordinate columns absent")
all_xy <- as.matrix(md[, xy_cols]); rownames(all_xy) <- colnames(obj)
led[, cell_id := as.character(cell_id)]
closed_ok <- tolower(as.character(led$closed)) %in% c("true", "1")
ref <- led[!cell_id %in% ids & state %in% c("defined_fine", "defined_broad_only") & confidence %in% c("high", "medium") & closed_ok & nzchar(as.character(broad_label))]
ref <- ref[!duplicated(cell_id) & cell_id %in% rownames(all_xy)]
ref[, broad_canonical := canon(broad_label)]
if (nrow(ref) < k) stop("insufficient spatial reference observations")
qxy <- all_xy[ids, , drop = FALSE]; rxy <- all_xy[ref$cell_id, , drop = FALSE]
sample_ids <- sample(rownames(all_xy), min(30000L, nrow(all_xy)))
grid_step <- median(nn2(all_xy[sample_ids, , drop = FALSE], all_xy[sample_ids, , drop = FALSE], k = 2)$nn.dists[, 2], na.rm = TRUE)
radius <- grid_step * radius_multiplier
nn <- nn2(rxy, qxy, k = min(k, nrow(rxy)))
labels <- matrix(ref$broad_canonical[nn$nn.idx], nrow = nrow(nn$nn.idx)); within <- nn$nn.dists <= radius
mp[, `:=`(spatial_neighbors_within_radius = rowSums(within), spatial_same_label_n = 0L, spatial_same_label_fraction = 0, nearest_defined_distance = nn$nn.dists[, 1])]
for (i in seq_len(nrow(mp))) {
  same <- labels[i, ] == mp$predicted_label_canonical[i] & within[i, ]
  mp$spatial_same_label_n[i] <- sum(same)
  mp$spatial_same_label_fraction[i] <- sum(same) / max(1L, sum(within[i, ]))
}
mp[, observed_density_spatial_consistent := spatial_neighbors_within_radius >= 5L & spatial_same_label_n >= 3L & spatial_same_label_fraction >= spatial_fraction]
mp[, calibrated_call := mapping_tier %in% c("high", "moderate") & tolower(as.character(meets_moderate_or_higher)) %in% c("true", "1", "yes")]
mp[, validated_broad_return := calibrated_call & fullfeature_program_consistent & cluster_context_consistent & observed_density_spatial_consistent]
mp[, `:=`(
  final_broad_label = ifelse(validated_broad_return, predicted_label_canonical, NA_character_),
  final_state = ifelse(validated_broad_return, "defined_broad_only", "qc_holdout"),
  final_confidence = ifelse(validated_broad_return & mapping_tier == "high", "high", ifelse(validated_broad_return, "medium", "low")),
  final_action = ifelse(validated_broad_return, "terminal_residual_qc_atlas_broad_rescue", "retain_terminal_qc_reject"),
  fine_anchor_eligible = FALSE, x = qxy[cell_id, 1], y = qxy[cell_id, 2]
)]

fwrite(mp, file.path(outdir, "tables", "qc_atlas_adjudication.tsv.gz"), sep = "\t")
sumtab <- mp[, .N, by = .(predicted_label_canonical, calibrated_call, fullfeature_program_consistent, cluster_context_consistent, observed_density_spatial_consistent, validated_broad_return, final_state, final_broad_label)][order(-N)]
fwrite(sumtab, file.path(outdir, "tables", "qc_atlas_adjudication_summary.tsv"), sep = "\t")
fwrite(context, file.path(outdir, "tables", "cluster_program_context.tsv"), sep = "\t")
plot_dt <- mp[, .(x, y, plot_label = ifelse(validated_broad_return, final_broad_label, "QC holdout"))]
pal <- c(Stromal = "#4DAF4A", `Immune/plasma` = "#E41A1C", `QC holdout` = "#D0D0D0")
p <- ggplot(plot_dt, aes(x, y, colour = plot_label)) + geom_point(size = 0.08) + scale_colour_manual(values = pal, na.value = "#984EA3") + scale_y_reverse() + coord_equal() + theme_void() + labs(title = "Post-recluster calibrated QC atlas adjudication", colour = NULL)
ggsave(file.path(outdir, "figures", "qc_atlas_adjudication.png"), p, width = 9, height = 8, dpi = 360, bg = "white")
ggsave(file.path(outdir, "figures", "qc_atlas_adjudication.pdf"), p, width = 9, height = 8, device = cairo_pdf, bg = "white")
manifest <- list(status = "PASS", n_query = nrow(mp), n_validated_broad_return = sum(mp$validated_broad_return), n_qc_retained = sum(!mp$validated_broad_return),
  cluster_membership_artifact = normalizePath(a$clusters), minimum_program_hits = min_program_hits,
  minimum_cluster_any_fraction = min_cluster_any, minimum_cluster_strong_fraction = min_cluster_strong,
  spatial_neighbors = k, observed_density_radius = radius, observed_grid_step = grid_step, minimum_spatial_fraction = spatial_fraction,
  thresholds_derived_from_current_query = TRUE, full_feature_validation = TRUE, observed_density_spatial_prior = TRUE,
  fine_anchor_eligible = FALSE, warning = "Broad-only rescue; rejected observations remain QC holdout and no atlas return becomes a subtype anchor.")
write_json(manifest, file.path(outdir, "provenance", "adjudication_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
capture.output(sessionInfo(), file = file.path(outdir, "provenance", "sessionInfo.txt"))
writeLines(c("status\tPASS", paste0("completed_at\t", format(Sys.time(), tz = "UTC", usetz = TRUE))), file.path(outdir, "RUN_COMPLETE.tsv"))
cat(toJSON(manifest, pretty = TRUE, auto_unbox = TRUE), "\n")
