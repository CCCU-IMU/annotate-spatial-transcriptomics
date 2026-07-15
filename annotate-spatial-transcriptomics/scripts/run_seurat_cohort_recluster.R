#!/usr/bin/env Rscript
# Current public entry point. The implementation filename is retained for
# byte-for-byte compatibility with projects created before v1.6.0.
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
if (length(file_arg) != 1L) stop("cannot resolve cohort runner location")
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
source(file.path(script_dir, "run_seurat_pool_recluster.R"), chdir = TRUE)
