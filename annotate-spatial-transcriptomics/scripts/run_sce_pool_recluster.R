#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(Seurat)
  library(SeuratObject); library(Matrix); library(data.table); library(ggplot2)
  library(scattermore)
})

parse_args <- function(x) { o <- list(); i <- 1L; while (i <= length(x)) { k <- sub("^--", "", x[i]); if (i == length(x) || startsWith(x[i+1], "--")) { o[[k]] <- TRUE; i <- i+1L } else { o[[k]] <- x[i+1]; i <- i+2L } }; o }
a <- parse_args(commandArgs(trailingOnly=TRUE)); need <- c("rds","membership","out","cell-id-col","resolutions"); miss <- need[!need %in% names(a)]; if(length(miss)) stop("Missing: ", paste(miss, collapse=", "))
dir.create(a$out,recursive=TRUE,showWarnings=FALSE); dir.create(file.path(a$out,"tables"),showWarnings=FALSE); dir.create(file.path(a$out,"figures"),showWarnings=FALSE)
read_any <- function(p) if(grepl("\\.gz$",p)) fread(cmd=paste("gzip -dc",shQuote(p))) else fread(p)
save_both <- function(p,stem,w=9,h=7){ggsave(paste0(stem,".png"),p,width=w,height=h,dpi=360,bg="white",limitsize=FALSE);ggsave(paste0(stem,".pdf"),p,width=w,height=h,device=cairo_pdf,bg="white",limitsize=FALSE)}
tag <- function(x) gsub("\\.","p",as.character(x)); seed <- as.integer(ifelse(is.null(a$seed),20260713,a$seed)); dims_n <- as.integer(ifelse(is.null(a$dims),30,a$dims)); nfeatures <- as.integer(ifelse(is.null(a$nfeatures),3000,a$nfeatures)); resolutions <- as.numeric(strsplit(a$resolutions,",",fixed=TRUE)[[1]]); method <- ifelse(is.null(a$normalization),"SCT",toupper(a$normalization)); set.seed(seed)

sce <- readRDS(a$rds); if(!inherits(sce,"SummarizedExperiment")) stop("Expected SingleCellExperiment/SummarizedExperiment")
mem <- read_any(a$membership); cc <- a$`cell-id-col`; mem[[cc]] <- as.character(mem[[cc]]); if(uniqueN(mem[[cc]]) != nrow(mem)) stop("Duplicate membership IDs")
if("query_or_anchor"%in%names(mem)&&any(mem$query_or_anchor=="anchor"))stop("SCE runner is query-only and cannot claim anchor assistance; use the full-feature Seurat or standard H5AD anchor runner")
cells <- intersect(colnames(sce),mem[[cc]]); if(length(cells) != nrow(mem)) stop("Membership IDs do not exactly match object")
count_assay <- ifelse(is.null(a$`count-assay`),"counts",a$`count-assay`); if(!count_assay %in% assayNames(sce)) stop("Count assay absent: ",count_assay)
cts <- assay(sce,count_assay)[,cells,drop=FALSE]; umi <- Matrix::colSums(cts); zero_ids <- names(umi)[!is.finite(umi)|umi<=0]; fwrite(data.table(cell_id=zero_ids,route="qc_holdout",reason="zero_count_in_selected_assay"),file.path(a$out,"tables","zero_count_observations.tsv"),sep="\t"); if(length(zero_ids)==length(umi))stop("All pool observations have zero counts"); cts <- cts[,setdiff(colnames(cts),zero_ids),drop=FALSE]; q <- CreateSeuratObject(counts=cts,assay="RNA",project="framework_pool",min.cells=0,min.features=0)
mem <- mem[match(colnames(q),get(cc))]; for (nm in setdiff(names(mem),cc)) q[[nm]] <- mem[[nm]]
if(method=="SCT") { q <- SCTransform(q,assay="RNA",new.assay.name="SCT",variable.features.n=nfeatures,verbose=FALSE,return.only.var.genes=FALSE); DefaultAssay(q) <- "SCT" } else { q <- NormalizeData(q,verbose=FALSE); q <- FindVariableFeatures(q,nfeatures=nfeatures,verbose=FALSE); q <- ScaleData(q,verbose=FALSE) }
max_pc <- min(max(dims_n,30),ncol(q)-1,nrow(q)-1); if(max_pc < 2) stop("Pool too small for graph reclustering")
q <- RunPCA(q,npcs=max_pc,verbose=FALSE); use_dims <- seq_len(min(dims_n,max_pc)); q <- FindNeighbors(q,dims=use_dims,verbose=FALSE); q <- RunUMAP(q,dims=use_dims,seed.use=seed,verbose=FALSE)
um <- Embeddings(q,"umap"); umdt <- data.table(cell_id=rownames(um),UMAP_1=um[,1],UMAP_2=um[,2]); coord_names <- if(all(c("sdimx","sdimy")%in%names(mem)))c("sdimx","sdimy") else if(all(c("x","y")%in%names(mem)))c("x","y") else NULL
for(res in resolutions){
  cn <- paste0("framework_res",tag(res)); q <- FindClusters(q,resolution=res,cluster.name=cn,random.seed=seed,verbose=FALSE); lab <- as.character(q[[cn,drop=TRUE]]); dt <- data.table(cell_id=colnames(q),cluster=lab,resolution=res); fwrite(dt,file.path(a$out,"tables",paste0(cn,"_clusters.tsv")),sep="\t")
  Idents(q) <- factor(lab,levels=sort(unique(lab))); deg <- if(length(unique(lab))>1) as.data.table(FindAllMarkers(q,assay=DefaultAssay(q),slot="data",only.pos=TRUE,test.use="wilcox",min.pct=.05,logfc.threshold=.1,return.thresh=1,max.cells.per.ident=2000,random.seed=seed,densify=FALSE,verbose=FALSE)) else data.table(p_val=numeric(),avg_log2FC=numeric(),pct.1=numeric(),pct.2=numeric(),p_val_adj=numeric(),cluster=character(),gene=character()); fwrite(deg,file.path(a$out,"tables",paste0(cn,"_DEG_all.tsv")),sep="\t")
  lfc <- intersect(c("avg_log2FC","avg_logFC"),names(deg))[1]; top <- deg[0]; if(nrow(deg) && length(lfc) && all(c("p_val_adj","cluster")%in%names(deg))){ sig <- deg[p_val_adj<.05 & get(lfc)>0]; top <- sig[order(cluster,p_val_adj,-get(lfc)),head(.SD,100),by=cluster] }; fwrite(top,file.path(a$out,"tables",paste0(cn,"_DEG_top100.tsv")),sep="\t")
  pal <- setNames(hcl.colors(length(unique(lab)),"Dynamic"),sort(unique(lab))); u <- merge(umdt,dt,by="cell_id",sort=FALSE); pu <- ggplot(u,aes(UMAP_1,UMAP_2,colour=cluster))+scattermore::geom_scattermore(pointsize=.8,pixels=c(1600,1600))+scale_colour_manual(values=pal)+theme_classic(base_size=8)+labs(title=paste("Pool UMAP",cn)); save_both(pu,file.path(a$out,"figures",paste0(cn,"_UMAP")))
  if(!is.null(coord_names)){ s <- cbind(mem[,c(cc,coord_names),with=FALSE],cluster=lab); setnames(s,c(cc,coord_names),c("cell_id","x","y")); ps <- ggplot(s,aes(x,y,colour=cluster))+scattermore::geom_scattermore(pointsize=.65,pixels=c(1800,1800))+scale_colour_manual(values=pal)+scale_y_reverse()+coord_equal()+theme_void()+labs(title=paste("Pool spatial",cn)); save_both(ps,file.path(a$out,"figures",paste0(cn,"_spatial"))); for(g in sort(unique(lab))){s[,selected:=cluster==g];ph<-ggplot(s,aes(x,y))+scattermore::geom_scattermore(data=s[selected==FALSE],colour="#DDDDDD",pointsize=.4,pixels=c(1600,1600))+scattermore::geom_scattermore(data=s[selected==TRUE],colour="#D62728",pointsize=.8,pixels=c(1600,1600))+scale_y_reverse()+coord_equal()+theme_void()+labs(title=paste0(cn," cluster ",g," (n=",sum(s$selected),")"));save_both(ph,file.path(a$out,"figures",paste0(cn,"_cluster_",g,"_highlight"))) } }
}
saveRDS(q,file.path(a$out,"pool_reclustered_seurat.rds"),compress=FALSE); fwrite(data.table(parameter=c("source_class","normalization","dims","nfeatures","resolutions","seed","anchor_assisted","query_only_graph_umap_deg","n_observations_input","n_observations_analyzed","n_zero_count_routed_qc"),value=c(class(sce)[1],method,length(use_dims),nfeatures,paste(resolutions,collapse=","),seed,FALSE,TRUE,length(cells),ncol(q),length(zero_ids))),file.path(a$out,"run_manifest.tsv"),sep="\t")
