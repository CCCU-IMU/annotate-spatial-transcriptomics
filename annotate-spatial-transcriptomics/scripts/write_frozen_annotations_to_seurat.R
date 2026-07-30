#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(Seurat);library(SeuratObject);library(data.table);library(jsonlite)})
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
dir.create(dirname(a$out),recursive=TRUE,showWarnings=FALSE)
threads<-as.integer(ifelse(is.null(a$threads),Sys.getenv("LSB_DJOB_NUMPROC",unset="1"),a$threads));if(!is.finite(threads)||threads<1L)threads<-1L
raw_tmp<-tempfile(pattern=paste0(basename(a$out),"."),tmpdir=dirname(a$out),fileext=".raw.tmp")
compressed_tmp<-paste0(a$out,".compressed.tmp")
pigz_error_tmp<-paste0(a$out,".pigz.stderr.tmp")
on.exit(unlink(c(raw_tmp,compressed_tmp,pigz_error_tmp),force=TRUE),add=TRUE)
saveRDS(obj,raw_tmp,compress=FALSE)
pigz<-Sys.which("pigz");compression_method<-"R_gzip";used_threads<-1L
if(nzchar(pigz)&&threads>1L){
  status<-system2(pigz,c("-p",as.character(threads),"-c",shQuote(raw_tmp)),stdout=compressed_tmp,stderr=pigz_error_tmp)
  if(!identical(as.integer(status),0L)){
    detail<-if(file.exists(pigz_error_tmp))paste(readLines(pigz_error_tmp,warn=FALSE),collapse="; ")else"no stderr"
    stop("pigz failed while writing the annotated RDS: ",detail)
  }
  compression_method<-"pigz";used_threads<-threads
}else{
  unlink(raw_tmp,force=TRUE);saveRDS(obj,compressed_tmp,compress="gzip")
}
if(!file.exists(compressed_tmp)||file.info(compressed_tmp)$size<=0)stop("annotated RDS temporary output is empty")
if(file.exists(a$out))unlink(a$out,force=TRUE)
if(!file.rename(compressed_tmp,a$out))stop("atomic annotated RDS rename failed")
sha256<-function(path){x<-system2("sha256sum",shQuote(path),stdout=TRUE);strsplit(x[[1L]],"\\s+")[[1L]][[1L]]}
manifest_path<-ifelse(is.null(a$`writer-manifest`),paste0(a$out,".writer_manifest.json"),a$`writer-manifest`)
dir.create(dirname(manifest_path),recursive=TRUE,showWarnings=FALSE)
write_json(list(status="PASS",stage="canonical_final_rds_writer",annotated_rds=list(path=normalizePath(a$out),sha256=sha256(a$out),n_bytes=unname(file.info(a$out)$size)),semantic_hash=as.character(a$`semantic-hash`),n_analysis=sum(is_analysis),n_excluded=sum(!is_analysis),compression_method=compression_method,compression_threads=used_threads,atomic_rename=TRUE),manifest_path,pretty=TRUE,auto_unbox=TRUE)
cat("Wrote frozen v2.2 annotations for ",sum(is_analysis)," analysis observations and ",sum(!is_analysis)," excluded observations\n",sep="")
