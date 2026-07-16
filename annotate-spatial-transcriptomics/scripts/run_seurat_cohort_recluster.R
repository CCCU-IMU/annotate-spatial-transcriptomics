#!/usr/bin/env Rscript
# Current public entry point for cohort reclustering.
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
if (length(file_arg) != 1L) stop("cannot resolve cohort runner location")
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
source(file.path(script_dir, "run_seurat_cohort_recluster_impl.R"), chdir = TRUE)
