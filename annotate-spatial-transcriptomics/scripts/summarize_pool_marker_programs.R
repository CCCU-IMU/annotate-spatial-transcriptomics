#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(Seurat);library(SeuratObject);library(Matrix);library(data.table)})
parse_args<-function(x){o<-list();i<-1L;while(i<=length(x)){k<-sub("^--","",x[i]);if(i==length(x)||startsWith(x[i+1],"--")){o[[k]]<-TRUE;i<-i+1L}else{o[[k]]<-x[i+1];i<-i+2L}};o}
read_any<-function(p)if(grepl("\\.gz$",p))fread(cmd=paste("gzip -dc",shQuote(p)))else fread(p)
a<-parse_args(commandArgs(trailingOnly=TRUE));need<-c("rds","clusters","markers","out","cell-id-col","cluster-col");miss<-need[!need%in%names(a)];if(length(miss))stop("Missing: ",paste(miss,collapse=", "))
obj<-readRDS(a$rds);cl<-read_any(a$clusters);mk<-fread(a$markers);cc<-a$`cell-id-col`;gc<-a$`cluster-col`;if(!all(c(cc,gc)%in%names(cl)))stop("cluster membership columns missing");if(!all(c("gene","program")%in%names(mk)))stop("markers require gene and program")
cl[[cc]]<-as.character(cl[[cc]]);if(uniqueN(cl[[cc]])!=nrow(cl)||any(!cl[[cc]]%in%colnames(obj)))stop("invalid cluster membership")
assay<-ifelse(is.null(a$assay),DefaultAssay(obj),a$assay);layer<-ifelse(is.null(a$layer),"data",a$layer);mat<-tryCatch(LayerData(obj[[assay]],layer=layer),error=function(e)GetAssayData(obj[[assay]],slot=layer));present<-intersect(unique(mk$gene),rownames(mat));if(!length(present))stop("none of the marker genes are present")
ids<-cl[[cc]];x<-mat[present,ids,drop=FALSE];groups<-as.character(cl[[gc]]);out<-rbindlist(lapply(sort(unique(groups)),function(g){idx<-which(groups==g);data.table(cluster=g,gene=present,avg_expression=Matrix::rowMeans(x[,idx,drop=FALSE]),pct_expressed=100*Matrix::rowMeans(x[,idx,drop=FALSE]>0),n_observations=length(idx))}),fill=TRUE)
out<-merge(out,unique(mk[,.(gene,program)]),by="gene",all.x=TRUE,allow.cartesian=TRUE);out[,gene_order:=match(gene,mk$gene)];setorder(out,program,gene_order,cluster);dir.create(dirname(a$out),recursive=TRUE,showWarnings=FALSE);fwrite(out,a$out,sep="\t");fwrite(data.table(gene=mk$gene,present=mk$gene%in%rownames(mat),program=mk$program),paste0(a$out,".feature_audit.tsv"),sep="\t");cat("PASS",nrow(cl),uniqueN(groups),length(present),"\n")
