#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(Seurat);library(SeuratObject);library(data.table)})
parse_args<-function(x){o<-list();i<-1L;while(i<=length(x)){k<-sub("^--","",x[[i]]);if(i==length(x)||startsWith(x[[i+1L]],"--")){o[[k]]<-TRUE;i<-i+1L}else{o[[k]]<-x[[i+1L]];i<-i+2L}};o}
read_any<-function(p)if(grepl("\\.gz$",p,ignore.case=TRUE))fread(cmd=paste("gzip -dc",shQuote(p)))else fread(p)
first_col<-function(d,names,default=""){hit<-names[names%in%colnames(d)];if(length(hit))as.character(d[[hit[[1L]]]])else rep(default,nrow(d))}
blank_to_na<-function(x){x<-as.character(x);x[is.na(x)|!nzchar(x)]<-NA_character_;x}
a<-parse_args(commandArgs(trailingOnly=TRUE));required<-c("rds","membership","out","semantic-hash");missing<-required[!required%in%names(a)];if(length(missing))stop("Missing: ",paste(missing,collapse=", "))
obj<-readRDS(a$rds);if(!inherits(obj,"Seurat"))stop("Frozen annotation writeback requires a Seurat object")
membership<-read_any(a$membership);if(!all(c("cell_id","final_broad_label","final_cell_type")%in%names(membership)))stop("Membership lacks cell_id/final_broad_label/final_cell_type");membership[,cell_id:=as.character(cell_id)];if(anyDuplicated(membership$cell_id)||any(!nzchar(membership$cell_id)))stop("Membership has invalid cell IDs")
excluded<-data.table();if(!is.null(a$`excluded-membership`)){excluded<-read_any(a$`excluded-membership`);if(!"cell_id"%in%names(excluded))stop("excluded-membership lacks cell_id");excluded[,cell_id:=as.character(cell_id)];if(anyDuplicated(excluded$cell_id)||any(!nzchar(excluded$cell_id)))stop("excluded-membership has invalid cell IDs")}
cells<-colnames(obj);if(length(intersect(membership$cell_id,excluded$cell_id)))stop("Analysis and excluded memberships overlap");if(!setequal(cells,c(membership$cell_id,excluded$cell_id)))stop("Membership plus excluded-membership must exactly cover the Seurat observation universe")
meta<-data.table(cell_id=cells);meta<-merge(meta,membership,by="cell_id",all.x=TRUE,sort=FALSE);meta[,`.input_order`:=match(cell_id,cells)];setorder(meta,`.input_order`);meta[,`.input_order`:=NULL];if(!identical(meta$cell_id,cells))stop("Cell order restoration failed")
is_analysis<-meta$cell_id%in%membership$cell_id;broad<-blank_to_na(first_col(meta,c("final_broad_label")));fine<-blank_to_na(first_col(meta,c("final_fine_label")));cell_type<-blank_to_na(first_col(meta,c("final_cell_type")));state<-blank_to_na(first_col(meta,c("state_annotations","final_state_annotation")))
annotation_state<-first_col(meta,c("final_state","final_annotation_state"));annotation_state[!is_analysis]<-"excluded_initial_qc"
obj$spanno_v2_2_broad<-broad
obj$spanno_v2_2_broad_display<-ifelse(is_analysis,ifelse(is.na(broad),"QC/Unknown",broad),"Excluded initial QC")
obj$spanno_v2_2_cell_type<-cell_type
obj$spanno_v2_2_cell_type_display<-ifelse(is_analysis,ifelse(is.na(cell_type),"QC/Unknown",cell_type),"Excluded initial QC")
obj$spanno_v2_2_broad_confidence<-blank_to_na(first_col(meta,c("final_broad_confidence","confidence")))
obj$spanno_v2_2_fine<-fine
obj$spanno_v2_2_fine_confidence<-blank_to_na(first_col(meta,c("final_fine_confidence")))
obj$spanno_v2_2_state<-state
obj$spanno_v2_2_annotation_state<-blank_to_na(annotation_state)
obj$spanno_v2_2_broad_source<-blank_to_na(first_col(meta,c("broad_freeze_source","assignment_origin")))
obj$spanno_v2_2_fine_source<-blank_to_na(first_col(meta,c("final_fine_assignment_source")))
obj$spanno_v2_2_semantic_hash<-rep(as.character(a$`semantic-hash`),length(cells))
if(!is.null(a$`sample-id`))obj$spanno_v2_2_sample_id<-rep(as.character(a$`sample-id`),length(cells))
dir.create(dirname(a$out),recursive=TRUE,showWarnings=FALSE);saveRDS(obj,a$out,compress=TRUE)
cat("Wrote frozen v2.2 annotations for ",sum(is_analysis)," analysis observations and ",sum(!is_analysis)," excluded observations\n",sep="")
