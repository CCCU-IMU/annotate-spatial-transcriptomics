#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(Seurat);library(SeuratObject);library(data.table);library(ggplot2);library(scattermore)})
parse_args<-function(x){o<-list();i<-1L;while(i<=length(x)){k<-sub("^--","",x[i]);if(i==length(x)||startsWith(x[i+1],"--")){o[[k]]<-TRUE;i<-i+1L}else{o[[k]]<-x[i+1];i<-i+2L}};o}
a<-parse_args(commandArgs(trailingOnly=TRUE));need<-c("rds","metadata","out","cell-id-col");miss<-need[!need%in%names(a)];if(length(miss))stop("Missing: ",paste(miss,collapse=", "))
dir.create(file.path(a$out,"figures"),recursive=TRUE,showWarnings=FALSE)
read_any<-function(p)if(grepl("\\.gz$",p))fread(cmd=paste("gzip -dc",shQuote(p)))else fread(p)
save_both<-function(p,stem){ggsave(paste0(stem,".png"),p,width=9,height=7,dpi=400,bg="white",limitsize=FALSE);ggsave(paste0(stem,".pdf"),p,width=9,height=7,device=cairo_pdf,bg="white",limitsize=FALSE)}
obj<-readRDS(a$rds);meta<-read_any(a$metadata);cc<-a$`cell-id-col`;meta[[cc]]<-as.character(meta[[cc]]);if(uniqueN(meta[[cc]])!=nrow(meta))stop("duplicate metadata IDs")
cells<-intersect(colnames(obj),meta[[cc]]);meta<-meta[match(cells,get(cc))];if("analysis_scope"%in%names(meta)){keep<-meta$analysis_scope=="analysis_set";cells<-cells[keep];meta<-meta[keep]}
if(!is.null(a$umap)){um<-read_any(a$umap);um[[cc]]<-as.character(um[[cc]]);um<-um[match(cells,get(cc))];ux<-intersect(c("UMAP_1","umap_1","UMAP1"),names(um))[1];uy<-intersect(c("UMAP_2","umap_2","UMAP2"),names(um))[1];ud<-data.table(cell_id=cells,x=um[[ux]],y=um[[uy]])}else{red<-ifelse(is.null(a$reduction),"umap",a$reduction);em<-Embeddings(obj,red)[cells,1:2,drop=FALSE];ud<-data.table(cell_id=cells,x=em[,1],y=em[,2])}
sd<-NULL;if(!is.null(a$coordinates)){co<-read_any(a$coordinates);co[[cc]]<-as.character(co[[cc]]);co<-co[match(cells,get(cc))];xx<-intersect(c("sdimx","x","spatial_x"),names(co))[1];yy<-intersect(c("sdimy","y","spatial_y"),names(co))[1];sd<-data.table(cell_id=cells,x=co[[xx]],y=co[[yy]])}else if(all(c("x","y")%in%colnames(obj[[]]))){md<-obj[[]];sd<-data.table(cell_id=cells,x=as.numeric(md[cells,"x"]),y=as.numeric(md[cells,"y"]))}
plot_one<-function(base,labels,title,stem,spatial=FALSE){d<-copy(base);d[,label:=as.character(labels)];d<-d[!is.na(label)&nzchar(label)];if(!nrow(d))return(FALSE);lev<-sort(unique(d$label));pal<-setNames(hcl.colors(length(lev),"Dynamic"),lev);p<-ggplot(d,aes(x,y,colour=label))+scattermore::geom_scattermore(pointsize=ifelse(spatial,.5,.65),pixels=c(2200,2200))+scale_colour_manual(values=pal)+theme_classic(base_size=8)+labs(title=title,colour=NULL);if(spatial)p<-p+scale_y_reverse()+coord_equal()+theme_void()+labs(title=title);save_both(p,stem);TRUE}
rows<-list();z<-1L
for(view in c("strict","inclusive","display")){
  bcol<-paste0(view,"_broad_label");fcol<-paste0(view,"_fine_label");scol<-paste0(view,"_state");if(!all(c(bcol,fcol,scol)%in%names(meta)))stop("missing ",view," view columns")
  fine<-as.character(meta[[fcol]]);fine[is.na(fine)|!nzchar(fine)]<-as.character(meta[[bcol]])[is.na(fine)|!nzchar(fine)]
  for(level in c("broad","subtype")){labels<-if(level=="broad")meta[[bcol]]else fine;stem<-file.path(a$out,"figures",paste0(view,"_",level,"_UMAP"));plot_one(ud,labels,paste(view,level,"UMAP"),stem,FALSE);if(!is.null(sd)){stem2<-file.path(a$out,"figures",paste0(view,"_",level,"_spatial"));plot_one(sd,labels,paste(view,level,"spatial"),stem2,TRUE)};rows[[z]]<-data.table(view=view,level=level,n_labeled=sum(!is.na(labels)&nzchar(labels)),umap_png=paste0(stem,".png"),umap_pdf=paste0(stem,".pdf"),spatial_png=if(!is.null(sd))paste0(stem2,".png")else "",spatial_pdf=if(!is.null(sd))paste0(stem2,".pdf")else "");z<-z+1L}
}
fwrite(rbindlist(rows),file.path(a$out,"tables","annotation_view_overview_asset_index.tsv"),sep="\t")
