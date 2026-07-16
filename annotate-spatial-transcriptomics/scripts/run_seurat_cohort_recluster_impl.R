#!/usr/bin/env Rscript
# Query-only Seurat cohort reclustering implementation.

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(data.table)
  library(ggplot2); library(patchwork); library(scattermore)
})

parse_args <- function(x) { out <- list(); i <- 1L; while (i <= length(x)) { k <- sub("^--", "", x[i]); if (i == length(x) || startsWith(x[i+1], "--")) { out[[k]] <- TRUE; i <- i+1L } else { out[[k]] <- x[i+1]; i <- i+2L } }; out }
a <- parse_args(commandArgs(trailingOnly = TRUE))
need <- c("rds","membership","out","cell-id-col","resolutions"); miss <- need[!need %in% names(a)]; if (length(miss)) stop("Missing: ", paste(miss, collapse=", "))
dir.create(a$out, recursive=TRUE, showWarnings=FALSE); dir.create(file.path(a$out,"tables"),showWarnings=FALSE); dir.create(file.path(a$out,"figures"),showWarnings=FALSE)
seed <- as.integer(ifelse(is.null(a$seed),20260710,a$seed)); set.seed(seed)
future_globals_max_gb <- as.numeric(ifelse(is.null(a$`future-globals-max-gb`),100,a$`future-globals-max-gb`))
options(future.globals.maxSize=future_globals_max_gb*1024^3)
if(requireNamespace("future",quietly=TRUE))future::plan("sequential")
dims_n <- as.integer(ifelse(is.null(a$dims),30,a$dims)); nfeatures <- as.integer(ifelse(is.null(a$nfeatures),3000,a$nfeatures))
pca_npcs <- as.integer(ifelse(is.null(a$`pca-npcs`),max(dims_n,30),a$`pca-npcs`))
sct_ncells_cap <- as.integer(ifelse(is.null(a$`sct-ncells`),50000,a$`sct-ncells`))
umap_min_dist <- as.numeric(ifelse(is.null(a$`umap-min-dist`),0.3,a$`umap-min-dist`))
resolutions <- as.numeric(strsplit(a$resolutions, ",", fixed=TRUE)[[1]])
if(!length(resolutions)||any(!is.finite(resolutions))||anyDuplicated(resolutions))stop("--resolutions must be a nonempty unique numeric grid")
minimum_resolution <- as.numeric(ifelse(is.null(a$`minimum-resolution`),0.1,a$`minimum-resolution`))
if(!is.finite(minimum_resolution)||minimum_resolution<0)stop("--minimum-resolution must be finite and nonnegative")
if(any(resolutions < minimum_resolution))stop("resolution below formal minimum ",minimum_resolution," is forbidden; repair the graph instead of lowering resolution")
resolution_contract <- tolower(ifelse(is.null(a$`resolution-contract`),"generic",a$`resolution-contract`))
if(!resolution_contract%in%c("generic","sheep_ovary"))stop("--resolution-contract must be generic or sheep_ovary")
if(resolution_contract=="sheep_ovary"&&!identical(resolutions,c(0.1,0.2,0.3,0.4,0.6)))stop("sheep_ovary cohort grid must equal 0.1,0.2,0.3,0.4,0.6")
detect_scheduler_cpus <- function(){for(name in c("LSB_DJOB_NUMPROC","SLURM_CPUS_PER_TASK","NSLOTS","AIP_CPUS")){value<-suppressWarnings(as.integer(Sys.getenv(name,unset="")));if(length(value)&&is.finite(value)&&value>0L)return(list(cpus=value,source=name))};list(cpus=NA_integer_,source="not_detected")}
scheduler_cpu<-detect_scheduler_cpus();resolution_workers_default<-if(is.finite(scheduler_cpu$cpus))min(length(resolutions),scheduler_cpu$cpus)else 1L
resolution_workers_requested <- as.integer(ifelse(is.null(a$`resolution-workers`),resolution_workers_default,a$`resolution-workers`))
resolution_future_plan_requested <- tolower(ifelse(is.null(a$`resolution-future-plan`),"auto",a$`resolution-future-plan`))
if(!is.finite(resolution_workers_requested)||resolution_workers_requested<1L)stop("--resolution-workers must be >= 1")
if(is.finite(scheduler_cpu$cpus)&&resolution_workers_requested>scheduler_cpu$cpus)stop("--resolution-workers exceeds scheduler allocation from ",scheduler_cpu$source)
if(!resolution_future_plan_requested%in%c("auto","multicore","multisession","sequential"))stop("--resolution-future-plan must be auto, multicore, multisession or sequential")
activate_resolution_parallelism <- function(requested_workers,n_tasks,requested_plan){
  workers<-max(1L,min(as.integer(requested_workers),as.integer(n_tasks)));plan_used<-"sequential"
  if(workers>1L&&requested_plan!="sequential"){
    if(!requireNamespace("future",quietly=TRUE))stop("future is required when --resolution-workers is greater than one")
    plan_used<-requested_plan;if(plan_used=="auto")plan_used<-if(future::supportsMulticore())"multicore"else"multisession"
    if(plan_used=="multicore"&&!future::supportsMulticore())stop("multicore future plan was requested but is unsupported on this runtime")
    if(plan_used=="multicore")future::plan(future::multicore,workers=workers)else future::plan(future::multisession,workers=workers)
  }else{workers<-1L;if(requireNamespace("future",quietly=TRUE))future::plan(future::sequential)}
  list(workers=workers,plan=plan_used)
}
method <- ifelse(is.null(a$normalization),"SCT",toupper(a$normalization))
role_col <- ifelse(is.null(a$`role-col`),"query_or_anchor",a$`role-col`); anchor_col <- ifelse(is.null(a$`anchor-label-col`),"anchor_label",a$`anchor-label-col`)
query_value <- ifelse(is.null(a$`query-value`),"query",a$`query-value`); anchor_value <- ifelse(is.null(a$`anchor-value`),"anchor",a$`anchor-value`)
read_any <- function(p) if (grepl("\\.gz$",p)) fread(cmd=paste("gzip -dc",shQuote(p))) else fread(p)
save_both <- function(p,stem,w,h){ggsave(paste0(stem,".png"),p,width=w,height=h,dpi=360,bg="white",limitsize=FALSE);ggsave(paste0(stem,".pdf"),p,width=w,height=h,device=cairo_pdf,bg="white",limitsize=FALSE)}
tag <- function(x) gsub("\\.","p",as.character(x))

obj <- readRDS(a$rds); stopifnot(inherits(obj,"Seurat")); mem <- read_any(a$membership); cc <- a$`cell-id-col`; mem[[cc]] <- as.character(mem[[cc]])
stopifnot(uniqueN(mem[[cc]])==nrow(mem)); missing_ids <- setdiff(mem[[cc]],colnames(obj)); if(length(missing_ids))stop(length(missing_ids)," membership IDs missing from object")
anchor_mode <- role_col %in% names(mem)
if(anchor_mode){
  roles <- as.character(mem[[role_col]]); if(any(!roles%in%c(query_value,anchor_value)))stop("membership role values must be query or anchor")
  query_ids <- mem[get(role_col)==query_value,get(cc)]; anchor_ids <- mem[get(role_col)==anchor_value,get(cc)]
  if(!length(query_ids)||!length(anchor_ids))stop("anchor-assisted mode requires nonempty query and anchor memberships")
  if(!anchor_col%in%names(mem)||any(!nzchar(as.character(mem[get(role_col)==anchor_value,get(anchor_col)]))))stop("anchors require anchor_label")
  if(uniqueN(as.character(mem[get(role_col)==anchor_value,get(anchor_col)]))<2L)stop("anchor-assisted interpretation requires at least two anchor labels")
}else{query_ids<-mem[[cc]];anchor_ids<-character()}
all_ids <- c(query_ids,anchor_ids); joint <- subset(obj,cells=all_ids); mem <- mem[match(colnames(joint),get(cc))]; stopifnot(all(mem[[cc]]==colnames(joint)))
assay <- ifelse(is.null(a$assay),DefaultAssay(joint),a$assay); DefaultAssay(joint)<-assay
count_layer <- ifelse(is.null(a$`count-layer`),"counts",a$`count-layer`); cts <- tryCatch(LayerData(joint[[assay]],layer=count_layer),error=function(e)GetAssayData(joint[[assay]],slot=count_layer)); umi <- Matrix::colSums(cts)
zero_ids <- names(umi)[!is.finite(umi)|umi<=0]; zero_query <- intersect(zero_ids,query_ids); zero_anchor <- intersect(zero_ids,anchor_ids)
fwrite(data.table(cell_id=zero_query,route="qc_holdout",reason="zero_count_query_in_selected_assay"),file.path(a$out,"tables","zero_count_observations.tsv"),sep="\t")
if(length(zero_anchor))warning(length(zero_anchor)," zero-count anchors were excluded")
query_ids<-setdiff(query_ids,zero_ids);anchor_ids<-setdiff(anchor_ids,zero_ids);all_ids<-c(query_ids,anchor_ids);if(!length(query_ids))stop("All query observations have zero counts");if(anchor_mode&&!length(anchor_ids))stop("All anchors have zero counts")
joint<-subset(joint,cells=all_ids);mem<-mem[match(colnames(joint),get(cc))]
if(length(query_ids)<3L)stop("At least three nonzero query observations are required")
k_param <- if(is.null(a$k)) min(30L,max(5L,floor(sqrt(length(query_ids))))) else as.integer(a$k)
k_param <- min(k_param,length(query_ids)-1L)

if (method=="SCT") {
  sct_method <- ifelse(is.null(a$`sct-method`),"glmGamPoi",a$`sct-method`)
  if(sct_method=="glmGamPoi"&&!requireNamespace("glmGamPoi",quietly=TRUE))stop("glmGamPoi is required by the frozen SCT profile; refusing a silent method fallback")
  joint <- NormalizeData(joint,assay=assay,normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)
  joint <- SCTransform(joint,assay=assay,new.assay.name="SCT",vst.flavor="v2",variable.features.n=nfeatures,ncells=min(sct_ncells_cap,ncol(joint)),conserve.memory=TRUE,verbose=FALSE,return.only.var.genes=TRUE,method=sct_method,seed.use=seed)
  DefaultAssay(joint)<-"SCT"
} else {
  joint <- NormalizeData(joint,assay=assay,verbose=FALSE); joint <- FindVariableFeatures(joint,assay=assay,nfeatures=nfeatures,verbose=FALSE); joint <- ScaleData(joint,assay=assay,verbose=FALSE)
}
pca_npcs_use <- min(pca_npcs,ncol(joint)-2L,length(VariableFeatures(joint))-1L)
if(pca_npcs_use<2L)stop("Too few PCs are available for cohort reclustering")
dims_use_n <- min(dims_n,pca_npcs_use)
joint <- RunPCA(joint,npcs=pca_npcs_use,features=VariableFeatures(joint),seed.use=seed,verbose=FALSE)
pca_joint <- Embeddings(joint,"pca")
if(anchor_mode){
  anchor_meta <- mem[match(anchor_ids,get(cc))]; labels <- as.character(anchor_meta[[anchor_col]]); centroids <- rowsum(pca_joint[anchor_ids,seq_len(dims_use_n),drop=FALSE],labels)/as.numeric(table(labels)[rownames(rowsum(pca_joint[anchor_ids,seq_len(dims_use_n),drop=FALSE],labels))])
  qemb <- pca_joint[query_ids,seq_len(dims_use_n),drop=FALSE]; dmat <- sapply(seq_len(nrow(centroids)),function(i)rowSums((qemb-matrix(centroids[i,],nrow=nrow(qemb),ncol=ncol(qemb),byrow=TRUE))^2)); if(is.null(dim(dmat)))dmat<-matrix(dmat,ncol=1)
  colnames(dmat)<-rownames(centroids);ord<-t(apply(dmat,1,order));top1<-colnames(dmat)[ord[,1]];top1d<-dmat[cbind(seq_len(nrow(dmat)),ord[,1])];top2d<-if(ncol(dmat)>1)dmat[cbind(seq_len(nrow(dmat)),ord[,2])]else rep(NA_real_,nrow(dmat))
  anchor_evidence<-data.table(cell_id=query_ids,nearest_anchor_label=top1,nearest_anchor_distance=top1d,anchor_distance_margin=top2d-top1d);fwrite(anchor_evidence,file.path(a$out,"tables","query_anchor_distance_evidence.tsv"),sep="\t")
}else anchor_evidence<-data.table(cell_id=query_ids)
q <- subset(joint,cells=query_ids)
q <- FindNeighbors(q,reduction="pca",dims=seq_len(dims_use_n),k.param=k_param,nn.method="annoy",n.trees=50,annoy.metric="cosine",graph.name=c("COHORT_nn","COHORT_snn"),verbose=FALSE)
resolution_parallel <- activate_resolution_parallelism(resolution_workers_requested,length(resolutions),resolution_future_plan_requested)
message("Running ",length(resolutions)," cohort Leiden resolutions with ",resolution_parallel$workers," future worker(s) using ",resolution_parallel$plan)
cluster_results <- FindClusters(object=q[["COHORT_snn"]],algorithm=4,resolution=resolutions,random.seed=seed,verbose=FALSE)
q <- RunUMAP(q,reduction="pca",dims=seq_len(dims_use_n),n.neighbors=k_param,min.dist=umap_min_dist,metric="cosine",seed.use=seed,verbose=FALSE)
if(requireNamespace("future",quietly=TRUE))future::plan(future::sequential)
qmem<-mem[match(query_ids,get(cc))];stopifnot(all(qmem[[cc]]==query_ids))
for(nm in setdiff(names(qmem),cc))q[[nm]]<-qmem[[nm]]
metadata_names <- colnames(q[[]])
if (!is.null(a$`x-col`) || !is.null(a$`y-col`)) {
  if (is.null(a$`x-col`) || is.null(a$`y-col`)) stop("--x-col and --y-col must be supplied together")
  spatial_pair <- c(a$`x-col`, a$`y-col`)
  if (!all(spatial_pair %in% metadata_names)) stop("requested spatial coordinate columns are missing")
} else {
  coordinate_candidates <- list(c("x", "y"), c("sdimx", "sdimy"), c("imagecol", "imagerow"), c("col", "row"))
  hits <- vapply(coordinate_candidates, function(pair) all(pair %in% metadata_names), logical(1))
  spatial_pair <- if (any(hits)) coordinate_candidates[[which(hits)[1]]] else character()
}
coords <- if (length(spatial_pair)) {
  md_coordinates <- q[[]]
  data.table(
    cell_id = rownames(md_coordinates),
    x = as.numeric(md_coordinates[[spatial_pair[[1]]]]),
    y = as.numeric(md_coordinates[[spatial_pair[[2]]]])
  )
} else NULL
um <- Embeddings(q,"umap"); umdt <- data.table(cell_id=rownames(um),UMAP_1=um[,1],UMAP_2=um[,2])
source_cols<-intersect(c("source_key","state_tags","spatial_tags","qc_tags","candidate_lineages"),names(qmem)); composition_all<-list();qc_all<-list();anchor_all<-list();zz<-1L
for (res in resolutions) {
  cn <- paste0("framework_res",tag(res)); result_col <- paste0("res.",res)
  if(!result_col%in%names(cluster_results))stop("missing parallel Leiden result: ",result_col)
  lab <- factor(as.character(cluster_results[[result_col]]),levels=sort(unique(as.character(cluster_results[[result_col]]))))
  q[[cn]] <- lab
  dt <- data.table(cell_id=colnames(q),cluster=lab,resolution=res); fwrite(dt,file.path(a$out,"tables",paste0(cn,"_clusters.tsv")),sep="\t")
  Idents(q)<-lab; deg <- if(length(unique(lab))>1)as.data.table(FindAllMarkers(q,assay=assay,slot="data",only.pos=TRUE,test.use="wilcox",min.pct=0.05,logfc.threshold=0.1,return.thresh=1,max.cells.per.ident=2000,random.seed=seed,densify=FALSE,verbose=FALSE))else data.table(p_val=numeric(),avg_log2FC=numeric(),pct.1=numeric(),pct.2=numeric(),p_val_adj=numeric(),cluster=character(),gene=character())
  if(!ncol(deg))deg<-data.table(p_val=numeric(),avg_log2FC=numeric(),pct.1=numeric(),pct.2=numeric(),p_val_adj=numeric(),cluster=character(),gene=character())
  fwrite(deg,file.path(a$out,"tables",paste0(cn,"_DEG_all.tsv")),sep="\t"); lfc_col<-intersect(c("avg_log2FC","avg_logFC"),names(deg))[1];top<-deg[0]
  if(nrow(deg)&&length(lfc_col)&&all(c("p_val_adj","cluster")%in%names(deg))){sig<-deg[p_val_adj<0.05 & get(lfc_col)>0];top<-sig[order(cluster,p_val_adj,-get(lfc_col)),head(.SD,100),by=cluster]};fwrite(top,file.path(a$out,"tables",paste0(cn,"_DEG_top100.tsv")),sep="\t")
  joined<-merge(dt,qmem,by.x="cell_id",by.y=cc,sort=FALSE);if(length(source_cols))for(sc in source_cols){by_cols<-c("cluster",as.character(sc));tmp<-joined[,.(n=.N),by=by_cols];tmp[,`:=`(resolution=res,field=sc,value=as.character(get(sc)))];tmp[,fraction:=n/sum(n),by=cluster];composition_all[[zz]]<-tmp[,.(resolution,cluster,field,value,n,fraction)];zz<-zz+1L}
  md<-q[[]];qc_candidates<-intersect(c("nCount_Spatial","nFeature_Spatial","nCount_RNA","nFeature_RNA","percent.mt"),names(md));if(length(qc_candidates)){qq<-data.table(cell_id=rownames(md),md[,qc_candidates,drop=FALSE]);qq<-merge(dt,qq,by="cell_id");qc_all[[length(qc_all)+1L]]<-qq[,lapply(.SD,median,na.rm=TRUE),by=.(resolution,cluster),.SDcols=qc_candidates]}
  if(anchor_mode){ae<-merge(dt,anchor_evidence,by="cell_id");anchor_all[[length(anchor_all)+1L]]<-ae[,.(n=.N,mean_distance=mean(nearest_anchor_distance),mean_margin=mean(anchor_distance_margin,na.rm=TRUE)),by=.(resolution,cluster,nearest_anchor_label)][,fraction:=n/sum(n),by=.(resolution,cluster)]}
  pal<-setNames(hcl.colors(nlevels(lab),"Dynamic"),levels(lab)); u<-merge(umdt,dt,by="cell_id",sort=FALSE)
  p<-ggplot(u,aes(UMAP_1,UMAP_2,colour=cluster))+scattermore::geom_scattermore(pointsize=.65,pixels=c(1800,1800))+scale_colour_manual(values=pal)+theme_classic(base_size=8)+labs(title=paste("Query-only cohort UMAP",cn));save_both(p,file.path(a$out,"figures",paste0(cn,"_UMAP")),9,7)
  if (!is.null(coords)) { s<-merge(coords,dt,by="cell_id",sort=FALSE); ps<-ggplot(s,aes(x,y,colour=cluster))+scattermore::geom_scattermore(pointsize=.5,pixels=c(1800,1800))+scale_colour_manual(values=pal)+scale_y_reverse()+coord_equal()+theme_void()+labs(title=paste("Query-only cohort spatial",cn)); save_both(ps,file.path(a$out,"figures",paste0(cn,"_spatial")),9,7); for(g in levels(lab)){s[,selected:=as.character(cluster)==g];ph<-ggplot(s,aes(x,y))+scattermore::geom_scattermore(data=s[selected==FALSE],colour="#DDDDDD",pointsize=.4,pixels=c(1600,1600))+scattermore::geom_scattermore(data=s[selected==TRUE],colour="#D62728",pointsize=.8,pixels=c(1600,1600))+scale_y_reverse()+coord_equal()+theme_void()+labs(title=paste0(cn," cluster ",g," (n=",sum(s$selected),")"));save_both(ph,file.path(a$out,"figures",paste0(cn,"_cluster_",g,"_highlight")),9,7)} }
}
if(length(composition_all))fwrite(rbindlist(composition_all,fill=TRUE),file.path(a$out,"tables","cluster_source_state_composition.tsv"),sep="\t")
if(length(qc_all))fwrite(rbindlist(qc_all,fill=TRUE),file.path(a$out,"tables","cluster_QC_summary.tsv"),sep="\t")
if(length(anchor_all))fwrite(rbindlist(anchor_all,fill=TRUE),file.path(a$out,"tables","cluster_anchor_distance_summary.tsv"),sep="\t")
fwrite(mem, file.path(a$out,"tables","analyzed_membership.tsv.gz"), sep="\t")
saveRDS(q,file.path(a$out,"cohort_reclustered_query_seurat.rds"),compress=FALSE)
if(anchor_mode)saveRDS(joint,file.path(a$out,"joint_query_anchor_pca_seurat.rds"),compress=FALSE)
fwrite(data.table(parameter=c("normalization","full_feature_deg_assay","full_feature_normalization","sct_vst_flavor","sct_method","sct_ncells","sct_conserve_memory","sct_return_only_var_genes","pca_npcs_requested","pca_npcs_used","dims_requested","dims_used","neighbor_k","neighbor_method","neighbor_trees","neighbor_metric","umap_min_dist","umap_metric","nfeatures","resolutions","minimum_resolution","resolution_contract","resolution_selection","seed","future_globals_max_gb","future_plan","scheduler_cpus_detected","scheduler_cpu_source","resolution_workers_requested","resolution_workers_used","umap_threads","anchor_assisted","query_only_graph_umap_deg","spatial_x_col","spatial_y_col","n_query_input","n_query_analyzed","n_anchors_analyzed","n_zero_query_qc","n_zero_anchor_excluded"),value=c(method,assay,"LogNormalize_scale_factor_10000",if(method=="SCT")"v2"else"not_applicable",if(method=="SCT")sct_method else "not_applicable",if(method=="SCT")min(sct_ncells_cap,ncol(joint))else"not_applicable",if(method=="SCT")TRUE else"not_applicable",if(method=="SCT")TRUE else "not_applicable",pca_npcs,pca_npcs_use,dims_n,dims_use_n,k_param,"annoy",50,"cosine",umap_min_dist,"cosine",nfeatures,paste(resolutions,collapse=","),minimum_resolution,resolution_contract,"adaptive_cohort_review_required",seed,future_globals_max_gb,resolution_parallel$plan,if(is.finite(scheduler_cpu$cpus))scheduler_cpu$cpus else "",scheduler_cpu$source,resolution_workers_requested,resolution_parallel$workers,resolution_parallel$workers,anchor_mode,TRUE,if(length(spatial_pair))spatial_pair[[1]]else"",if(length(spatial_pair))spatial_pair[[2]]else"",length(query_ids)+length(zero_query),ncol(q),length(anchor_ids),length(zero_query),length(zero_anchor))),file.path(a$out,"run_manifest.tsv"),sep="\t")
capture.output(sessionInfo(), file=file.path(a$out,"sessionInfo.txt"))
writeLines(c("status\tPASS",paste0("completed_at\t",format(Sys.time(),tz="UTC",usetz=TRUE))),file.path(a$out,"RUN_COMPLETE.tsv"))
quit(save="no",status=0,runLast=FALSE)
