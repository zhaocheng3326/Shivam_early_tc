#' ---
#' title: "SV_ETC(time_course) IT"
#' output:
#'  html_document:
#'    code_folding: hide
#' ---


# ##loading R library
rm(list=ls())
rewrite=FALSE

#' check whether in local computer
if (grepl("KI-",Sys.info()['nodename'])) {
  print("local computer")
  source("/Users/cheng.zhao/chzhao_bioinfo/PC/SnkM/SgCell.R")
  base_dir <- "/Users/cheng.zhao/Documents"
} else {
  print("On server")
  condaENV <- "/home/chenzh/miniconda3/envs/R4.3" 
  LBpath <- paste0(condaENV ,"/lib/R/library")
  #.libPaths(LBpath)
  base_dir="/home/chenzh"
}



suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(Matrix)
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(tibble)
  library(scran)
  library(batchelor)
  library(SeuratWrappers)
  #library(slingshot)
  #library(tradeSeq)
  #library(ComplexHeatmap)
})


set.seed(1)
DIR <- paste0(base_dir,"/My_project/Shivam_early_tc")
knitr::opts_knit$set(root.dir=DIR)
setwd(DIR)


## define the EM brief structure data
#' Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
#source("src/local.quick.fun.R")
#source("src/figures.setting.R")

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)


options(digits = 4)
options(future.globals.maxSize= 3001289600)
TD="Nov_2026"

rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter


load("tmp_data/gene.meta.Rdata",verbose=T)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


human.mk.list <- list(alpha=c("GCG","TTR","PCSK2"), beta=c("INS","ENO1","EDARADD"),delta=c("SST","HHEX","LEPR"),S4P=c("FOXJ1"),epsilon=c("GHRL","NPY1R","LINC00852","VTN"),gamma=c("PPY","GPC5-AS1","THSD7A"),co=c("IAPP","ADCYAP1","NPTX2"), ductal=c("KRT19","SPP1","HNF1B"),acinar=c("PRSS1","PNLIP"),endothelial=c("GNG11","SPARCL","SNHPG7","SPARCL1"),PSC=c("CALD1","VIM","LGALS1"), endocrine=c("CHGB","SCG5","CPE"), CC=c("TOP2A","CDK1"),CCI=c("CDKN1C"),alpha_prec=c("INS","GCG","ARX"),SCEC_p=c("CHGA","TPH1","PCSK1","NKX6.1","TPH1","DDC","SLC18A1","LMX1A","ADRΑ2A","FEV","TAC1","CXCL14"), SCEC_N=c("G6PC2","NPTX2","ISL1","PDX1"))
human.mk.list.short <- human.mk.list %>% lapply(function(x) {head(x,3)})


#' for RNA 
if (file.exists(paste0("tmp_data/",TD,"/data.ob.rds"))) {
  data.ob <- readRDS(paste0("tmp_data/",TD,"/data.ob.rds"))
  data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/data.ob.umap.rds"))
  data.ob@meta.data$EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"EML"]
  temp.DM.list <- readRDS(paste0("tmp_data/",TD,"/data.mk.rds"))
  
}else{
  #'
  counts.filter<- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
  meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) #%>% rows_update(readRDS(paste0("tmp_data/",TD,"/data.umap.ob.rds")) %>% select(cell,EML),by="cell") %>% rows_update(readRDS(paste0("tmp_data/",TD,"/D15.only.data.umap.ob.rds")) %>% select(cell,EML),by="cell")
  
  # #' extra counts
  # extra.meta.filter <- readRDS("~/My_project/Shivam_SunS_scRNA/tmp_data/Oct_2025/data.umap.ob.rds") %>% select(cell:mt.perc,EML) %>% filter(devTime=="scpan_unS") %>% mutate(EML=recode(EML,"Early_beta"="beta","Late_beta"="beta"),devTime="D42_unS") %>% select(colnames(meta.filter))
  # extra.counts.filter <- readRDS("~/My_project/Shivam_SunS_scRNA/tmp_data/Oct_2025/counts.filter.rds")[,extra.meta.filter$cell]
  
  #' confirm
  #table(rownames(extra.counts.filter)==rownames(counts.filter ))
  #table(duplicated(meta.filter %>% bind_rows(extra.meta.filter) %>% pull(cell)))
  
  # merge
  #meta.filter <- meta.filter %>% bind_rows(extra.meta.filter)
  #counts.filter <- counts.filter %>% cbind(extra.counts.filter)
  
  
  # if (file.exists(paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))) {
  #   lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))
  # }else{
  #   expG.set <- list()
  #   for (b in unique(meta.filter$devTime  %>% unique() %>% as.vector())) {
  #     temp.cell <- meta.filter %>% filter(devTime==b) %>% pull(cell)
  #     expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]
  #   }
  #   sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()
  #   length(sel.expG )
  #   
  #   sce.ob <- list()
  #   for (b in unique(meta.filter$devTime  %>% unique() %>% as.vector())) {
  #     print(b)
  #     temp.M <- meta.filter %>% filter(devTime==b)
  #     temp.sce <-  SingleCellExperiment(list(counts=as.matrix(counts.filter[sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors()
  #     sce.ob[[b]] <- temp.sce
  #   }
  #   
  #   mBN.sce.ob <- multiBatchNorm(sce.ob$D15,sce.ob$D18,sce.ob$D20,sce.ob$D22,sce.ob$D25)
  #   lognormExp.mBN<- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)
  #   saveRDS(lognormExp.mBN,paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))
  # }
  lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))
  
  sel.od <- c("D14","D15","D16","D17","D18")
  temp.M <- meta.filter
  temp.sel.expG <- rownames(lognormExp.mBN)
  
  
  #data.temp <- CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>% ScaleData(verbose=F)%>% RunPCA(verbose=F) %>% RunUMAP(dims=1:25,verbose=F) %>% FindNeighbors( dims = 1:25,verbose = FALSE) %>%  FindClusters(reso=1,verbose=F)#,vars.to.regress=c("nGene"),nn.method="annoy",annoy.metric="cosine"
  
  data.merge <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) #%>% CellCycleScoring( s.features = s.genes, g2m.features = g2m.genes)
  data.merge@assays$RNA$data <- as.matrix(lognormExp.mBN[temp.sel.expG,rownames(data.merge@meta.data)])
  data.spt <- SplitObject(data.merge, split.by = "devTime")%>% lapply(function(x){x=FindVariableFeatures(x,verbose=F,nfeatures=2000)})
  data.spt <- data.spt[sel.od]
  
  nGene=2500;pc=25;
  mnn.VGs <- SelectIntegrationFeatures(object.list = data.spt, nfeatures = nGene)
  for (b in names(data.spt)) {
    VariableFeatures(data.spt[[b]]) <- mnn.VGs
  }
  
  #' release memory
  rm(data.merge)
  rm(counts.filter)
  rm(sce.ob)
  rm(mBN.sce.ob)
  rm(lognormExp.mBN)
  
  data.spt <- data.spt[sel.od]
  set.seed(123)
  data.ob <- RunFastMNN(data.spt,verbose=F,features=length(mnn.VGs)) %>% RunUMAP( reduction = "mnn", dims = 1:pc,verbose=F) %>% FindNeighbors( reduction = "mnn", dims = 1:pc)#nn.method="annoy",annoy.metric="cosine")
  data.temp <- data.ob %>% FindClusters(reso=1,verbose=F)
  data.temp <- JoinLayers(data.temp)
  
  plot_grid(
    DimPlot(data.temp,label=T)+NoAxes()+NoLegend(),
    DimPlot(data.temp,group.by="EML",label=T)+NoAxes()+NoLegend(),
    DimPlot(data.temp,group.by="devTime",label=T)+NoAxes()+NoLegend(),
    FeaturePlot(data.temp,"nGene")+NoAxes()+NoLegend(),
    DimPlot(data.temp,group.by="batch",label=T)+NoAxes()+NoLegend(),
    FeaturePlot(data.temp,"mt.perc")+NoAxes()+NoLegend()
    #DimPlot(data.temp,label=T,group.by = "Phase")+NoAxes()+NoLegend()
  )
  FeaturePlot(data.temp,c("NEUROG3","FEV","NEUROD1","CHGA","INSM1","PDX1","NKX6.1","ARX","TPH1","INS","GCG","IRX1","IRX2","TTR","CHGB","PAX6","LMX1A"))
  FeaturePlot(data.temp,c("ANXA2","GSTA1","FZD5","ONECUT1","TPH1","SST","INS","GCG","CHGA"))
  
  DimPlot(data.temp,label=T,split.by = "pj",group.by = "EML")+NoLegend()
  DimPlot(data.temp,label=T,split.by = "devTime",group.by = "EML")+NoLegend()
  DimPlot(data.temp,label=T,split.by = "pj")+NoLegend()
  plot_grid(plotlist = FunFP_plot(data.temp,unlist(human.mk.list.short)))
  
  #' direct UMAP
  # temp.M <- meta.filter 
  # temp.sel.expG <- rownames(counts.filter)[rowSums(counts.filter[,c(temp.M$cell)] >=1) >2]
  # 
  # data.ob <- CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) %>% CellCycleScoring( s.features = s.genes, g2m.features = g2m.genes)  %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE)  %>% ScaleData(verbose=F)%>% RunPCA(verbose=F) %>% RunUMAP(dims=1:25,verbose=F) %>% FindNeighbors( dims = 1:25,verbose = FALSE)  #"mt.perc",#%>% CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes)  #vars.to.regress=c("mt.perc","S.Score", "G2M.Score")
  # #data.ob@assays$RNA@data <- as.matrix(lognormExp.mBN[temp.sel.expG,colnames(data.ob)])
  # 
  # data.temp <- data.ob %>% FindClusters(reso=0.4,verbose=F)
  # plot_grid(
  #   DimPlot(data.temp,label=T)+NoAxes()+NoLegend(),
  #   DimPlot(data.temp,group.by="EML",label=T)+NoAxes()+NoLegend(),
  #   DimPlot(data.temp,group.by="devTime",label=T)+NoAxes()+NoLegend(),
  #   FeaturePlot(data.temp,"nGene")+NoAxes()+NoLegend(),
  #   DimPlot(data.temp,group.by="batch",label=T)+NoAxes()+NoLegend(),
  #   FeaturePlot(data.temp,"mt.perc")+NoAxes()+NoLegend(),
  #   DimPlot(data.temp,label=T,group.by = "Phase")+NoAxes()+NoLegend()
  # )
  data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df()  %>% select(cell,SID:mt.perc,EML) %>% mutate(SC=paste0("C",as.vector(Idents(data.temp))))  %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$mnn@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:mnn_10),by="cell")
  
  
  data.ob.umap <- data.ob.umap %>% mutate(EML=recode(SC,'C0'="beta",'C1'="SCEC",'C2'="beta",'C3'="SCEC",'C4'="SCEC",'C5'="alpha",'C6'="beta",'C7'="SCEC_Beta_P","C8"="exo_C8","C9"="delta","C10"="exo_C10","C11"="exo_C10","C12"="exo_C12"))
  #data.ob.umap <- data.ob.umap %>% mutate(EML=recode(SC,'C0'="beta",'C1'="SCEC",'C2'="beta",'C3'="beta",'C4'="SCEC_Beta_P",'C5'="SCEC_Beta_P",'C6'="alpha",'C7'="alpha","C8"="delta","C9"="exo","C10"="SCEC_Beta_P","C11"="alpha","C12"="SCEC_Beta_P"))
  
  #' output 
  data.ob.umap %>% write.table(paste0("tmp_data/",TD,"/data.ob.umap.tsv"),quote=F,sep="\t",row.names=F)
  #saveRDS(data.ob,file=paste0("tmp_data/",TD,"/data.ob.rds"))
  #saveRDS(data.ob.umap,file=paste0("tmp_data/",TD,"/data.ob.umap.rds"))
  
  
  
  
  #####' marker ### need check
  data.temp@meta.data$EML <- (data.ob.umap %>% column_to_rownames("cell"))[rownames(data.temp@meta.data),"EML"]
  #Idents(data.temp) <- factor(data.temp@meta.data$EML)
  Idents(data.temp) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"EML"])
  
  DimPlot(data.temp,group.by="EML",label=T)
  
  temp.DM.list <- list()
  data.deg <- data.temp
  temp.DM.list$seperate <- FindAllMarkers(data.deg,only.pos=T)  %>% tbl_df() %>% filter(p_val_adj < 0.05) %>% group_by(gene) %>% top_n(1,-1*p_val_adj) %>% ungroup() %>% inner_join(log1p(AverageExpression(data.deg)$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df() %>% gather(cluster,ave_exp,-gene) %>% group_by(gene) %>% top_n(1,ave_exp) %>% mutate(cluster=gsub("-","_",cluster)),by = c("cluster", "gene"))
  
  
 # data.deg <- subset(data.temp,cell=(data.ob.umap %>% filter(!EML %in% c("beta_P","BSProG","SCEC_P","SCEC_Beta_P")) %>% filter(devTime!="D15")  %>% pull(cell)))
#temp.DM.list$exclude.D15.seperate <-FindAllMarkers(data.deg,only.pos=T)  %>% tbl_df() %>% filter(p_val_adj < 0.05) %>% group_by(gene) %>% top_n(1,-1*p_val_adj) %>% ungroup() %>% inner_join(log1p(AverageExpression(data.deg)$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df() %>% gather(cluster,ave_exp,-gene) %>% group_by(gene) %>% top_n(1,ave_exp) %>% mutate(cluster=gsub("-","_",cluster)),by = c("cluster", "gene"))
  
  saveRDS(temp.DM.list,paste0("tmp_data/",TD,"/data.mk.rds"))
  saveRDS(data.ob,paste0("tmp_data/",TD,"/data.ob.rds"))
  saveRDS(data.ob.umap,paste0("tmp_data/",TD,"/data.ob.umap.rds"))
}


FeaturePlot(data.temp,c("NEUROG3","FEV","NEUROD1","INSM1","PDX1","NKX6.1","ARX","TPH1","INS","GCG","IRX1","IRX2","TTR","CHGB","PAX6","LMX1A"))
FeaturePlot(data.temp,c("CHGA","TPH1","SST","INS","GCG","INS","GCG","ANXA2","GSTA1","FZD5","ONECUT1","LGALS3","PRSS1","COL2A1"))
FeaturePlot(data.temp,c("CHGA","TPH1","SST","INS","GCG","INS","GCG","ANXA2","GSTA1","FZD5","ONECUT1"))



data.temp <- data.ob %>% FindClusters(reso=1,verbose=F)
data.temp <- JoinLayers(data.temp)

FeaturePlot(data.ob,c("NEUROG3","PAX","ARX"))
temp.DM.list$seperate %>% filter(cluster=="SCEC_Beta_P") %>% pull(gene) %>% head(6) %>% FeaturePlot(data.temp,.)
temp.DM.list$seperate %>% filter(cluster=="SCEC_Beta_P") %>% pull(gene) %>% head(25) %>% FunFP_plot(data.temp,.) %>% cowplot::plot_grid(plotlist = )
temp.DM.list$seperate %>% filter(cluster=="11") %>% pull(gene) %>% head(6) %>% FeaturePlot(data.temp,.)

FeaturePlot(data.temp,c("IGF2","KLF2","NEUROG3"))




plot_grid(
  DimPlot(data.temp %>% FindClusters(reso=1.5,verbose=F),label=T)+NoAxes()+NoLegend(),
  DimPlot(data.temp,group.by="EML",label=T)+NoAxes()+NoLegend(),
  DimPlot(data.temp,group.by="devTime",label=T)+NoAxes()+NoLegend(),
  FeaturePlot(data.temp,"nGene")+NoAxes()+NoLegend(),
  DimPlot(data.temp,group.by="batch",label=T)+NoAxes()+NoLegend(),
  FeaturePlot(data.temp,"mt.perc")+NoAxes()+NoLegend(),
  FeaturePlot(data.ob,c("NEUROG3"))+NoAxes()+NoLegend(),
  FeaturePlot(data.ob,c("ARX"))+NoAxes()+NoLegend(),
  FeaturePlot(data.ob,c("PAX4"))+NoAxes()+NoLegend()
  #DimPlot(data.temp,label=T,group.by = "Phase")+NoAxes()+NoLegend()
)
c("NEUROG3", "NEUROD1", "NKX6.1", "PDX1") %>% setdiff(rownames(data.temp@assays$RNA))

FeaturePlot(data.temp,c("NEUROG3", "NEUROD1", "NKX6-1", "PDX1"))

FeaturePlot(data.temp,c("NEUROG3","IGF2","KLF2"))
