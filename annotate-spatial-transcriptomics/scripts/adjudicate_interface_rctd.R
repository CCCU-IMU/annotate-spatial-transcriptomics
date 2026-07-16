#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(data.table)
  library(jsonlite); library(RANN); library(ggplot2)
})
parse_args <- function(x) { o <- list(); i <- 1L; while (i <= length(x)) { k <- sub("^--", "", x[i]); if (i == length(x) || startsWith(x[i + 1L], "--")) { o[[k]] <- TRUE; i <- i + 1L } else { o[[k]] <- x[i + 1L]; i <- i + 2L } }; o }
read_any <- function(p) if (grepl("\\.gz$", p)) fread(cmd = paste("gzip -dc", shQuote(p))) else fread(p)
canon <- function(broad, fine = "") {
  x <- paste(as.character(broad), as.character(fine)); out <- rep("Other", length(x))
  out[grepl("stromal", x, ignore.case = TRUE)] <- "Stromal"
  out[grepl("mural|perivascular", x, ignore.case = TRUE)] <- "Stromal/vascular-associated"
  out[grepl("vascular|endothelial", x, ignore.case = TRUE)] <- "Vascular"
  out[grepl("granulosa", x, ignore.case = TRUE)] <- "Granulosa"
  out[grepl("theca", x, ignore.case = TRUE)] <- "Theca"
  out[grepl("immune|plasma|macroph|neutroph", x, ignore.case = TRUE)] <- "Immune"
  out[grepl("epithelial|mesothelial|surface epitheli", x, ignore.case = TRUE)] <- "Epithelial/mesothelial"
  out[grepl("oocyte", x, ignore.case = TRUE)] <- "Oocyte"
  out
}
a <- parse_args(commandArgs(trailingOnly = TRUE)); need <- c("rds", "rctd-results", "clusters", "cell-ledger", "out", "cell-id-col", "cluster-col")
miss <- need[!need %in% names(a)]; if (length(miss)) stop("Missing: ", paste(miss, collapse = ", "))
cc <- a$`cell-id-col`; gc <- a$`cluster-col`; outdir <- a$out
dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE); dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE); dir.create(file.path(outdir, "provenance"), recursive = TRUE, showWarnings = FALSE)
k <- as.integer(ifelse(is.null(a$`spatial-neighbors`), 20, a$`spatial-neighbors`)); radius_multiplier <- as.numeric(ifelse(is.null(a$`radius-multiplier`), 10, a$`radius-multiplier`)); min_spatial_fraction <- as.numeric(ifelse(is.null(a$`min-spatial-fraction`), 0.25, a$`min-spatial-fraction`)); set.seed(as.integer(ifelse(is.null(a$seed), 20260713, a$seed)))
obj <- readRDS(a$rds); d <- read_any(a$`rctd-results`); cl <- read_any(a$clusters); led <- read_any(a$`cell-ledger`)
d[, cell_id := as.character(cell_id)]; cl[[cc]] <- as.character(cl[[cc]])
if (uniqueN(d$cell_id) != nrow(d) || uniqueN(cl[[cc]]) != nrow(cl) || !all(d$cell_id %in% cl[[cc]]) || !all(d$cell_id %in% colnames(obj))) stop("query boundary mismatch")
if (!"rctd_confidence_tier" %in% names(d) || any(!d$rctd_confidence_tier %in% c("high", "moderate", "low"))) stop("tiered RCTD evidence is required")
d <- merge(d, cl[, .(cell_id = get(cc), cluster = as.character(get(gc)))], by = "cell_id", all.x = TRUE, sort = FALSE)
type_map <- c(Collagen_ECM_stromal_anchor = "Stromal", Collagen_stromal_anchor = "Stromal", Mural_perivascular_anchor = "Stromal/vascular-associated", Steroidogenic_theca_anchor = "Theca", Granulosa_anchor = "Granulosa", Blood_endothelial_anchor = "Vascular", Immune_myeloid_anchor = "Immune", Epithelial_mesothelial_anchor = "Epithelial/mesothelial", Strict_oocyte_anchor = "Oocyte", Strict_oocyte_seed_anchor = "Oocyte")
d[, proposed_canonical := unname(type_map[as.character(first_type)])]

programs <- list(
  Stromal = c("DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "COL6A1", "PDGFRA", "TCF21", "NR2F2", "OGN", "SPARC", "MGP", "CXCL14", "IGFBP5"),
  Mural = c("RGS5", "CSPG4", "MCAM", "PDGFRB", "ACTA2", "TAGLN", "MYL9", "MYH11", "CAVIN1"),
  Granulosa = c("FOXL2", "FSHR", "FST", "AMH", "CDH2", "INHA", "INHBA", "INHBB", "GJA1", "HSD17B1", "CYP19", "CYP19A1", "SERPINE2", "IHH"),
  Theca = c("LHCGR", "CYP17A1", "STAR", "CYP11A1", "HSD3B1", "INSL3", "NR5A1", "FDX1"),
  Oocyte_core = c("DDX4", "DAZL", "FIGLA", "NOBOX", "PADI6", "OOEP", "NPM2"),
  Oocyte_zona = c("ZP2", "ZP3", "ZP4"),
  Plasma = c("JCHAIN", "CD74", "XBP1", "MZB1", "TENT5C"),
  Immune = c("PTPRC", "CD53", "LYZ", "TYROBP", "FCER1G", "GPNMB", "ACP5", "S100A8", "S100A9", "PGLYRP1"),
  Epithelial = c("KRT8", "KRT18", "KRT19", "EPCAM", "MSLN", "PRG4"),
  Endothelial = c("PECAM1", "CDH5", "VWF", "CLDN5", "KDR", "FLT1", "EMCN", "RAMP2")
)
cts <- tryCatch(LayerData(obj[["Spatial"]], layer = "counts"), error = function(e) GetAssayData(obj[["Spatial"]], slot = "counts")); ids <- d$cell_id
for (nm in names(programs)) { g <- intersect(programs[[nm]], rownames(cts)); d[[paste0(nm, "_hits")]] <- if (length(g)) as.integer(Matrix::colSums(cts[g, ids, drop = FALSE] > 0)) else 0L }
d[, rctd_program_consistent := FALSE]
d[proposed_canonical == "Stromal" & first_type == "Collagen_ECM_stromal_anchor", rctd_program_consistent := Stromal_hits >= 3 & Stromal_hits >= pmax(Granulosa_hits, Theca_hits, Immune_hits, Epithelial_hits)]
d[proposed_canonical == "Stromal" & first_type == "Mural_perivascular_anchor", rctd_program_consistent := Mural_hits >= 3 & Mural_hits >= pmax(Endothelial_hits, Immune_hits, Epithelial_hits)]
d[proposed_canonical == "Granulosa", rctd_program_consistent := Granulosa_hits >= 4 & Granulosa_hits >= pmax(Theca_hits, Immune_hits, Epithelial_hits)]
d[proposed_canonical == "Theca", rctd_program_consistent := Theca_hits >= 3 & Theca_hits >= pmax(Immune_hits, Epithelial_hits)]
d[proposed_canonical == "Vascular", rctd_program_consistent := Endothelial_hits >= 2 & Endothelial_hits >= pmax(Stromal_hits, Mural_hits, Immune_hits, Epithelial_hits)]
d[proposed_canonical == "Immune", rctd_program_consistent := Immune_hits >= 2 & Immune_hits >= pmax(Epithelial_hits, Endothelial_hits)]
d[proposed_canonical == "Epithelial/mesothelial", rctd_program_consistent := Epithelial_hits >= 2 & Epithelial_hits >= pmax(Immune_hits, Endothelial_hits)]
d[, strict_oocyte_program := Oocyte_core_hits >= 5 & Oocyte_zona_hits >= 2 & (Oocyte_core_hits + Oocyte_zona_hits) >= 7 & pmax(Granulosa_hits, Stromal_hits, Immune_hits, Epithelial_hits) <= 1]
d[proposed_canonical == "Oocyte", rctd_program_consistent := strict_oocyte_program]
d[, plasma_review_candidate := cluster == "6" & Plasma_hits >= 3 & Plasma_hits > Immune_hits]

context <- rbindlist(list(
  d[, .(cluster_n = .N, program_any_fraction = mean(Stromal_hits >= 1), program_strong_fraction = mean(Stromal_hits >= 3)), by = cluster][, proposed_canonical := "Stromal"],
  d[, .(cluster_n = .N, program_any_fraction = mean(Granulosa_hits >= 1), program_strong_fraction = mean(Granulosa_hits >= 4)), by = cluster][, proposed_canonical := "Granulosa"],
  d[, .(cluster_n = .N, program_any_fraction = mean(Theca_hits >= 1), program_strong_fraction = mean(Theca_hits >= 3)), by = cluster][, proposed_canonical := "Theca"],
  d[, .(cluster_n = .N, program_any_fraction = mean((Oocyte_core_hits + Oocyte_zona_hits) >= 1), program_strong_fraction = mean(strict_oocyte_program)), by = cluster][, proposed_canonical := "Oocyte"],
  d[, .(cluster_n = .N, program_any_fraction = mean(Endothelial_hits >= 1), program_strong_fraction = mean(Endothelial_hits >= 2)), by = cluster][, proposed_canonical := "Vascular"],
  d[, .(cluster_n = .N, program_any_fraction = mean(Immune_hits >= 1), program_strong_fraction = mean(Immune_hits >= 2)), by = cluster][, proposed_canonical := "Immune"],
  d[, .(cluster_n = .N, program_any_fraction = mean(Epithelial_hits >= 1), program_strong_fraction = mean(Epithelial_hits >= 2)), by = cluster][, proposed_canonical := "Epithelial/mesothelial"]
))
d <- merge(d, context, by = c("cluster", "proposed_canonical"), all.x = TRUE, sort = FALSE)
d[, cluster_context_consistent := fifelse(proposed_canonical == "Oocyte", program_strong_fraction >= 0.005, program_any_fraction >= 0.15 & program_strong_fraction >= 0.03, na = FALSE)]

md <- obj@meta.data; xy_cols <- if (all(c("x", "y") %in% names(md))) c("x", "y") else if (all(c("sdimx", "sdimy") %in% names(md))) c("sdimx", "sdimy") else stop("spatial coordinates absent")
all_xy <- as.matrix(md[, xy_cols]); rownames(all_xy) <- colnames(obj); qxy <- all_xy[ids, , drop = FALSE]
led[, cell_id := as.character(cell_id)]; closed_ok <- tolower(as.character(led$closed)) %in% c("true", "1")
ref <- led[!cell_id %in% ids & state %in% c("defined_fine", "defined_broad_only") & confidence %in% c("moderate", "moderate") & closed_ok & cell_id %in% rownames(all_xy)]
ref <- ref[!duplicated(cell_id)]; ref[, spatial_canonical := canon(broad_label, fine_label)]; ref <- ref[spatial_canonical != "Other"]
rxy <- all_xy[ref$cell_id, , drop = FALSE]; sample_ids <- sample(rownames(all_xy), min(30000L, nrow(all_xy))); grid_step <- median(nn2(all_xy[sample_ids, , drop = FALSE], all_xy[sample_ids, , drop = FALSE], k = 2)$nn.dists[, 2], na.rm = TRUE); radius <- grid_step * radius_multiplier
nn <- nn2(rxy, qxy, k = min(k, nrow(rxy))); lab <- matrix(ref$spatial_canonical[nn$nn.idx], nrow = nrow(nn$nn.idx)); within <- nn$nn.dists <= radius
d[, `:=`(spatial_neighbors_within_radius = rowSums(within), spatial_same_label_n = 0L, spatial_same_label_fraction = 0, nearest_defined_distance = nn$nn.dists[, 1])]
for (i in seq_len(nrow(d))) { same <- lab[i, ] == d$proposed_canonical[i] & within[i, ]; d$spatial_same_label_n[i] <- sum(same); d$spatial_same_label_fraction[i] <- sum(same) / max(1L, sum(within[i, ])) }
d[, observed_density_spatial_consistent := spatial_neighbors_within_radius >= 5 & spatial_same_label_n >= 3 & spatial_same_label_fraction >= min_spatial_fraction]
d[proposed_canonical == "Oocyte", observed_density_spatial_consistent := observed_density_spatial_consistent & nearest_defined_distance <= radius / 2]
d[, post_rctd_validated := rctd_confidence_tier %in% c("high", "moderate") & rctd_program_consistent & cluster_context_consistent & observed_density_spatial_consistent & !plasma_review_candidate]
d[, final_action := "route_to_qc_holdout"]
d[plasma_review_candidate == TRUE, final_action := "targeted_plasma_recluster"]
d[post_rctd_validated & proposed_canonical != "Oocyte", final_action := "validated_broad_return"]
d[post_rctd_validated & proposed_canonical == "Oocyte", final_action := "targeted_oocyte_recluster"]
d[, final_broad_label := fifelse(final_action == "validated_broad_return" & first_type == "Mural_perivascular_anchor", "Stromal/vascular-associated",
  fifelse(final_action == "validated_broad_return" & proposed_canonical == "Granulosa", "Follicular somatic",
  fifelse(final_action == "validated_broad_return", proposed_canonical, NA_character_)))]
d[, final_state := fifelse(final_action == "validated_broad_return", "defined_broad_only", fifelse(final_action == "route_to_qc_holdout", "qc_holdout", "pending_review"))]
d[, `:=`(fine_anchor_eligible = FALSE, x = qxy[cell_id, 1], y = qxy[cell_id, 2])]
fwrite(d, file.path(outdir, "tables", "interface_rctd_adjudication.tsv.gz"), sep = "\t")
fwrite(d[final_action == "route_to_qc_holdout"], file.path(outdir, "tables", "rctd_qc_holdout_reroute_membership.tsv.gz"), sep = "\t")
sumtab <- d[, .N, by = .(cluster, first_type, validated_singlet, rctd_program_consistent, cluster_context_consistent, observed_density_spatial_consistent, strict_oocyte_program, plasma_review_candidate, post_rctd_validated, final_action, final_broad_label, final_state)][order(-N)]
fwrite(sumtab, file.path(outdir, "tables", "interface_rctd_adjudication_summary.tsv"), sep = "\t"); fwrite(context, file.path(outdir, "tables", "cluster_program_context.tsv"), sep = "\t")
plot_dt <- d[, .(x, y, plot_label = fifelse(final_action == "validated_broad_return", final_broad_label, fifelse(final_action == "targeted_oocyte_recluster", "Targeted Oocyte recluster", fifelse(final_action == "targeted_plasma_recluster", "Targeted plasma recluster", "Final QC holdout"))))]
pal <- c(Stromal = "#4DAF4A", `Stromal/vascular-associated` = "#377EB8", `Follicular somatic` = "#FF7F00", Theca = "#984EA3", Vascular = "#1F78B4", Immune = "#E31A1C", `Epithelial/mesothelial` = "#A65628", `Targeted Oocyte recluster` = "#FB9A99", `Targeted plasma recluster` = "#B15928", `Final QC holdout` = "#D0D0D0")
p <- ggplot(plot_dt, aes(x, y, colour = plot_label)) + geom_point(size = 0.08) + scale_colour_manual(values = pal) + scale_y_reverse() + coord_equal() + theme_void() + labs(title = "Calibrated interface RCTD post-adjudication", colour = NULL)
ggsave(file.path(outdir, "figures", "interface_rctd_post_adjudication.png"), p, width = 9, height = 8, dpi = 360, bg = "white"); ggsave(file.path(outdir, "figures", "interface_rctd_post_adjudication.pdf"), p, width = 9, height = 8, device = cairo_pdf, bg = "white")
manifest <- list(status = "PASS", n_query = nrow(d), rctd_high_n = sum(d$rctd_confidence_tier == "high"), rctd_moderate_n = sum(d$rctd_confidence_tier == "moderate"), rctd_low_n = sum(d$rctd_confidence_tier == "low"), rctd_fine_return_n = 0, rctd_broad_return_n = sum(d$final_action == "validated_broad_return"), n_qc_holdout_rerouted = sum(d$final_action == "route_to_qc_holdout"), n_targeted_oocyte_recluster = sum(d$final_action == "targeted_oocyte_recluster"), n_targeted_plasma_recluster = sum(d$final_action == "targeted_plasma_recluster"), oocyte_policy = "RCTD cannot directly define Oocyte; the candidate enters the independent contamination-safe targeted cohort", observed_grid_step = grid_step, observed_density_radius = radius, minimum_spatial_fraction = min_spatial_fraction, independent_fine_evidence = FALSE, fine_anchor_eligible = FALSE, warning = "RCTD is low-priority assistance: high plus independent evidence may support fine, moderate may return broad-only, and every low or otherwise unresolved observation enters the final frozen QC holdout for terminal calibrated Atlas review without QC reclustering.")
write_json(manifest, file.path(outdir, "provenance", "adjudication_manifest.json"), pretty = TRUE, auto_unbox = TRUE); capture.output(sessionInfo(), file = file.path(outdir, "provenance", "sessionInfo.txt")); writeLines(c("status\tPASS", paste0("completed_at\t", format(Sys.time(), tz = "UTC", usetz = TRUE))), file.path(outdir, "RUN_COMPLETE.tsv")); cat(toJSON(manifest, pretty = TRUE, auto_unbox = TRUE), "\n")
