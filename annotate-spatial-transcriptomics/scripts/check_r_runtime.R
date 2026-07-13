#!/usr/bin/env Rscript
pkgs<-c("Matrix","data.table","ggplot2","patchwork","scattermore","Seurat","SeuratObject","SingleCellExperiment","SummarizedExperiment","jsonlite","spacexr","igraph","uwot")
ok<-vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)
cat("R=",R.version.string,"\n",sep="");for(i in seq_along(pkgs))cat(pkgs[i],"\t",ok[i],"\n",sep="")
quit(status=ifelse(all(ok[c("Matrix","data.table","ggplot2")]),0,2))
