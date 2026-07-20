#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(Seurat);library(SeuratObject);library(Matrix);library(data.table);library(ggplot2);library(scattermore);library(patchwork)})
parse_args<-function(x){o<-list();i<-1L;while(i<=length(x)){k<-sub("^--","",x[i]);if(i==length(x)||startsWith(x[i+1],"--")){o[[k]]<-TRUE;i<-i+1L}else{o[[k]]<-x[i+1];i<-i+2L}};o};a<-parse_args(commandArgs(trailingOnly=TRUE))
need<-c("rds","coordinates","markers","out","cell-id-col");m<-need[!need%in%names(a)];if(length(m))stop("Missing: ",paste(m,collapse=", "))
dir.create(a$out,recursive=TRUE,showWarnings=FALSE);gene_dir<-file.path(a$out,"spatial_genes");dir.create(gene_dir,showWarnings=FALSE);dir.create(file.path(a$out,"tables"),showWarnings=FALSE);scope<-ifelse(is.null(a$scope),"all_input_observations",a$scope)
read_any<-function(p)if(grepl("\\.gz$",p))fread(cmd=paste("gzip -dc",shQuote(p)))else fread(p);safe<-function(x)substr(gsub("_+","_",gsub("[^A-Za-z0-9_.-]+","_",x)),1,150)
save_both<-function(p,stem,width=7.5,height=7){ggsave(paste0(stem,".png"),p,width=width,height=height,dpi=400,bg="white");ggsave(paste0(stem,".pdf"),p,width=width,height=height,device=cairo_pdf,bg="white")}
obj<-readRDS(a$rds);co<-read_any(a$coordinates);cc<-a$`cell-id-col`;co[[cc]]<-as.character(co[[cc]]);cells<-intersect(colnames(obj),co[[cc]]);co<-co[match(cells,get(cc))];xx<-intersect(c("sdimx","x","spatial_x"),names(co))[1];yy<-intersect(c("sdimy","y","spatial_y"),names(co))[1]
if(inherits(obj,"Seurat")){assay<-ifelse(is.null(a$assay),DefaultAssay(obj),a$assay);layer<-ifelse(is.null(a$layer),"data",a$layer);mat<-tryCatch(LayerData(obj[[assay]],layer=layer),error=function(e)GetAssayData(obj[[assay]],slot=layer))[,cells,drop=FALSE]}else if(inherits(obj,"SingleCellExperiment")||inherits(obj,"SummarizedExperiment")){if(!requireNamespace("SummarizedExperiment",quietly=TRUE))stop("SummarizedExperiment required");available<-SummarizedExperiment::assayNames(obj);data_name<-if(!is.null(a$`data-assay`))a$`data-assay` else if("normcounts"%in%available)"normcounts" else if("logcounts"%in%available)"logcounts" else stop("Specify --data-assay");mat<-SummarizedExperiment::assay(obj,data_name)[,cells,drop=FALSE]}else stop("Unsupported object class")
markers<-read_any(a$markers);stopifnot(all(c("gene","marker_group")%in%names(markers)));markers<-unique(markers[,.(gene=as.character(gene),marker_group=as.character(marker_group))])
rows<-list();group_plots<-list();z<-1L
for(i in seq_len(nrow(markers))){
  g<-markers$gene[i];group<-markers$marker_group[i]
  if(!g%in%rownames(mat)){
    rows[[z]]<-data.table(marker_group=group,gene=g,availability="missing_from_expression_matrix",scope=scope,n_observations=length(cells),png="",pdf="");z<-z+1L;next
  }
  d<-data.table(x=co[[xx]],y=co[[yy]],value=as.numeric(mat[g,]));cap<-quantile(d$value,.995,na.rm=TRUE);if(!is.finite(cap)||cap<=0)cap<-max(d$value,na.rm=TRUE)
  p<-ggplot(d,aes(x,y,colour=pmin(value,cap)))+scattermore::geom_scattermore(pointsize=.55,pixels=c(1800,1800))+scale_colour_viridis_c(option="magma",name="expression")+scale_y_reverse()+coord_equal()+theme_void()+labs(title=paste0(group," · ",g),subtitle=paste0("all observations; fixed point size; q99.5 clip=",signif(cap,3)))
  stem<-file.path(gene_dir,paste0(safe(group),"__",safe(g)));save_both(p,stem)
  rows[[z]]<-data.table(marker_group=group,gene=g,availability="available",scope=scope,n_observations=length(cells),png=paste0(stem,".png"),pdf=paste0(stem,".pdf"));z<-z+1L
  group_plots[[group]]<-c(group_plots[[group]],list(p))
}
asset_index<-rbindlist(rows,fill=TRUE);fwrite(asset_index,file.path(a$out,"tables","spatial_gene_asset_index.tsv"),sep="\t")
panel_dir<-file.path(a$out,"spatial_gene_panels");dir.create(panel_dir,showWarnings=FALSE);panel_rows<-list();z<-1L
for(group in unique(markers$marker_group)){
  plots<-group_plots[[group]];requested<-markers[marker_group==group,gene];available<-asset_index[marker_group==group&availability=="available",gene]
  stem<-file.path(panel_dir,safe(group))
  if(length(plots)){
    ncol<-min(3L,length(plots));panel<-patchwork::wrap_plots(plots,ncol=ncol)+patchwork::plot_annotation(title=paste0(group," canonical marker expression across all observations"));save_both(panel,stem,width=5.2*ncol,height=4.8*ceiling(length(plots)/ncol))
    png<-paste0(stem,".png");pdf<-paste0(stem,".pdf")
  }else{png<-"";pdf<-""}
  panel_rows[[z]]<-data.table(marker_group=group,requested_genes=paste(requested,collapse=";"),available_genes=paste(available,collapse=";"),n_requested=length(requested),n_available=length(available),scope=scope,n_observations=length(cells),point_size=.55,expression_scale="input normalized expression; per-gene q99.5 clip",png=png,pdf=pdf);z<-z+1L
}
fwrite(rbindlist(panel_rows,fill=TRUE),file.path(a$out,"tables","spatial_gene_group_asset_index.tsv"),sep="\t")
