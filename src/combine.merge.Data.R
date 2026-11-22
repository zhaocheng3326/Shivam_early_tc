#' ---
#' title: "combine all data"
#' output:
#'  html_document:
#'    code_folding: hide
#' ---


# ##loading R library
#R4.3
rm(list=ls())
condaENV <- "/home/chenzh/miniconda3/envs/R4.0"
#LBpath <- paste0(condaENV ,"/lib/R/library")
#.libPaths(LBpath)

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(scran)
  #library(batchelor)
  library(Seurat)
  #library(SeuratWrappers)
  #library(scuttle)
  #library(SeuratDisk)
  #library(uwot)
})


# working directory
DIR <- "~/My_project/Shivam_early_tc/"
knitr::opts_knit$set(root.dir=DIR)
setwd(DIR)


#' Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
source("src/local.quick.fun.R")


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


#' loading gene meta based on mouse genome
gtf.anno <- read.delim( "/home/chenzh/Genome_new/Human/RefSeq/clean_chr_refdata/genes/gene.gtf.anno",stringsAsFactors = F,row.names = 5,head=F)
mt.gene <- gtf.anno$V6[gtf.anno$V1=="chrMT"]
ribo.gene <- c(gtf.anno$V6[grepl("^RPS",gtf.anno$V6)],gtf.anno$V6[grepl("^RPL",gtf.anno$V6)],gtf.anno$V6[grepl("^MRPL",gtf.anno$V6)],gtf.anno$V6[grepl("^MRPS",gtf.anno$V6)])


dup.gene <- names(table(gtf.anno$V6))[table(gtf.anno$V6 ) > 1]


counts <- list()
metas <- list()

pj="SV_ETC"
temp.list <- list()

sa <- c()
for (s in c("D14","D15","D16","D17","D18")) {
  temp.list[[s]] <- (Read10X(data.dir = paste0("tmp_data/X10/SV_ETC/outs/per_sample_outs/",s,"/count/sample_feature_bc_matrix/")))[["Gene Expression"]]
}




GI=rownames(temp.list[["D15"]] )
SI=names(temp.list)
temp <- NULL
for ( n in names(temp.list)) {
  colnames(temp.list[[n]]) <- paste(n,"10X",colnames(temp.list[[n]]),sep="_") 
  temp.list[[n]] <- temp.list[[n]] %>% as.data.frame() %>% tbl_df() #%>% tibble::rownames_to_column("Gene")
}
counts[[pj]] <- data.frame(Gene=GI) %>% tbl_df() %>% mutate(Gene=as.vector(Gene)) %>% bind_cols(do.call("bind_cols",temp.list)) ## each sample has the same rownames
rm(temp.list)

temp.list <- list()
CI=colnames(counts[[pj]] )
for ( n in SI) {
  temp.list[[n]] <- data.frame(SID=paste(n,"10X",sep="_"),batch="batch1",seqType="10X",cellType=n,devTime=n,EML=n,cell=CI[grepl(n,CI)],subCT="None") %>% tbl_df()
}

metas[[pj]] <-  bind_rows(do.call("bind_rows",temp.list))
metas[[pj]] <- metas[[pj]]  %>% mutate(libsize=colSums(counts[[pj]] %>% select(-Gene) )[metas[[pj]]$cell]) %>% mutate(nGene=colSums((counts[[pj]]  %>% select(-Gene)) >0)[metas[[pj]]$cell], pj=pj)  





GI <- unique(gtf.anno$V6)

for (n in c("SV_ETC" )) {
  counts[[n]] <- (counts[[n]] %>% filter(Gene %in% c(dup.gene,paste0(dup.gene,".1")))  %>% mutate(Gene=sub("\\.1","",Gene)) %>% gather(cell,counts,-Gene)  %>% group_by(Gene,cell) %>% summarise(counts=sum(counts)) %>% spread(cell,counts) %>% select(colnames(counts[[n]])) %>% bind_rows(counts[[n]] %>%  filter(!Gene %in% c(dup.gene,paste0(dup.gene,".1"))))  %>% tibble::column_to_rownames("Gene"))[GI,]
  metas[[n]] <- metas[[n]] %>% left_join(data.frame(cell=colnames(counts[[n]]),mt.perc=colSums(counts[[n]][mt.gene,])/colSums(counts[[n]])) %>% tbl_df() ,by="cell")
}




counts.all <-do.call("bind_cols",counts)
meta.all <- do.call("bind_rows",metas)

meta.all <- meta.all  %>% replace(.,is.na(.),0)

#' loading the mouse mapping reads prop
#mouse.prop <- read.delim("tmp_data/mapping_data/mouse_check/UMI.extract.mapped.mouse.stat.tsv",stringsAsFactors = F,head=F) %>% tbl_df() %>% setNames(c("CB","MCT"))

#meta.all$CB <-meta.all$cell %>% gsub("SPEN_KO_bulk_Naive_10X_","",.) %>% gsub("SPEN_KO_CL5_Primed_10X_","",.) %>% gsub("WT_Primed_10X_","",.) %>% gsub("WT_RB_Naive_10X_","",.) %>% strsplit("-") %>% lapply(function(x) {x[1]}) %>% unlist()
#meta.all <- meta.all %>% left_join(mouse.prop,by="CB") %>% mutate(MCT.prop =MCT/libsize)

#' saving files
save(counts.all,meta.all,file=paste0("tmp_data/",TD,"/all.counts.meta.Rdata"))
save(gtf.anno,mt.gene,ribo.gene,file="tmp_data/gene.meta.Rdata")


