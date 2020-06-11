## STAGE 3 
## CODE 1
#### Getting Patient specific DE/DM genes  ####
# Created By: Banabithi Bose
# Date Created: 4/27/2020

#### Use BRCA FPKM data (before deconvolution) to get median and fold change #####
##Compute median expression of a gene across all samples
.libPaths(c("/users/home/bbose/R/x86_64-redhat-linux-gnu-library/3.5","/users/home/bbose/R/x86_64-redhat-linux-gnu-library/3.6","/usr/lib64/R/library","/usr/share/R/library"))
#.libPaths()
gc(reset=T)
library('statmod')
library("stringr")
library(data.table)
library(sqldf)
library(reshape)
library(reshape2)
library(stringr)
library(janitor)
library(stringr)
library(plyr)
library(readr)
library(gsubfn)
library(pathfindR)
library(parallel)

load("~/DRUG_PIPELINE1/DATASETS/TCGA/OV.Rda")## TCGA data
load("~/PATIENT_CELL_PIPELINE1/DATAFRAMES/Ensemble_hgnc.Rda")
GeneName<-data.table(row.names(RnaSeq_fpkm_data))#56537
RnaSeq_fpkm_data<-data.table(RnaSeq_fpkm_data)
RnaSeq_fpkm_data$ensembl_gene_id<-GeneName$V1
RnaSeq_fpkm_hgnc_data <- merge(Ensemble_hgnc,RnaSeq_fpkm_data,by = "ensembl_gene_id")
RnaSeq_fpkm_hgnc_data<-RnaSeq_fpkm_hgnc_data[,-1]#56537

library(matrixStats)
RnaSeq_fpkm_hgnc_data$Median_expr <- rowMedians(as.matrix(RnaSeq_fpkm_hgnc_data[,-1]),na.rm=T)#37109 genes total in hgnc

n<-ncol(RnaSeq_fpkm_hgnc_data)
processMethFiles <- function(fileNames){

  #SampleName <- fileNames[1]
  SampleName <- fileNames
  gene_expr<-data.table(RnaSeq_fpkm_hgnc_data[,SampleName])
  gene_med<-RnaSeq_fpkm_hgnc_data[,c(1,n)]
  gene_expr_med<-cbind(gene_med,gene_expr)
  colnames(gene_expr_med)<-c("gene","median_expr","expr")
  gene_expr_med[gene_expr_med$median_expr==0,2]<-0.00000000001 ## to avoid divison by 0 for fold change
  gene_expr_med[gene_expr_med$expr==0,3]<-0.00000000001 ## to get non zero fold change, removing fold change = 0 for down gene
  gene_expr_med$fold_change<-gene_expr_med$expr/gene_expr_med$median_expr
  ### Since fold change >2 or <.5 giving 10,000 genes among 40,000, I am taking 4 and .25 threshold,
  # new threshold is giving around 5000 genes
  upgene<-gene_expr_med[gene_expr_med$fold_change>=4,]## 4968 genes
  downgene<-gene_expr_med[gene_expr_med$fold_change<=0.25,]## 7068 genes
  pat_de_genes <- data.table(rbind(upgene,downgene))[,1]


  save(upgene,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Patient_DE_genes/UP_Gene/",str_sub(SampleName,1,16),"_upgene.Rda"))
  save(downgene,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Patient_DE_genes/DOWN_Gene/",str_sub(SampleName,1,16),"_downgene.Rda"))
  save(pat_de_genes,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Patient_DE_genes/DE_Gene/",str_sub(SampleName,1,16),"_pat_de_genes.Rda"))
  
}

print("cluster making started")
cl <- makeCluster(mc <- getOption("cl.cores", 70),type = "FORK")
invisible(clusterEvalQ(cl, library('glmnet')))
invisible(clusterEvalQ(cl, library(janitor)))
invisible(clusterEvalQ(cl, library(stringr)))
invisible(clusterEvalQ(cl, library(plyr)))
invisible(clusterEvalQ(cl, library(data.table)))
invisible(clusterEvalQ(cl, library(readr)))
invisible(clusterEvalQ(cl, library(gsubfn)))
invisible(clusterEvalQ(cl, library(parallel)))
invisible(clusterEvalQ(cl, library(sqldf)))
invisible(clusterEvalQ(cl, library(reshape)))
invisible(clusterEvalQ(cl, library(reshape2)))
invisible(clusterEvalQ(cl, library(pbapply)))

print("cluster export started")
parallel::clusterExport(cl=cl,c("RnaSeq_fpkm_hgnc_data"),envir=environment())#envir needed to be correct, export from global aswell as local
print("cluster making done")


fileNames<-colnames(RnaSeq_fpkm_hgnc_data[,-c(1,n)])#373 patients, they may not unique

parLapplyLB(cl,fileNames,processMethFiles)

stopCluster(cl)
print("cluster stopped")




## Checked intersect of de genes between two patients, it was around 2000 genes
## After this we need to prepare deconvoluted expression data with these high/low genes

##### Use BRCA Methylation data (before deconvolution) to get median and fold change #####
#Compute median expression of a gene across all samples 

# .libPaths(c("/users/home/bbose/R/x86_64-redhat-linux-gnu-library/3.5","/users/home/bbose/R/x86_64-redhat-linux-gnu-library/3.6","/usr/lib64/R/library","/usr/share/R/library"))
# #.libPaths()
# gc(reset=T)
# library('statmod')
# library("stringr")
# library(data.table)
# library(sqldf)
# library(reshape)
# library(reshape2)
# library(stringr)
# library(janitor)
# library(stringr)
# library(plyr)
# library(readr)
# library(gsubfn)
# library(pathfindR)
# library(parallel)
# 
# load("~/MIRDRIVER/PIPELINE2_OV_LASSO/DATAFRAMES/OV_All_methylation.Rda")
# x<-data.table(row.names(OV_All_methylation))
# All_methylation<-cbind(x,OV_All_methylation)
# All_methylation[All_methylation==0]<-0.00000000001 ## to avoid the log transformation of 0
# load("~/PATIENT_CELL_PIPELINE1/DATAFRAMES/Ensemble_hgnc.Rda")
# 
# #Writing a function for transforming beta values to M values
# Mtranform <- function(x) {log(x[1]/(1-x[1]))}
# n<-ncol(All_methylation)-1
# All_methylation[, 2:n] <- as.data.frame(lapply(All_methylation[, 2:n], FUN = function(x) {sapply(x, FUN = Mtranform)}))
# 
# hgnc_methylation <- merge(Ensemble_hgnc,All_methylation,by.x = "ensembl_gene_id",by.y = "V1")
# hgnc_methylation<-hgnc_methylation[,-1]
# 
# library(matrixStats)
# hgnc_methylation$Median_meth <- rowMedians(as.matrix(hgnc_methylation[,-1]),na.rm=T)#14808 genes total in hgnc
# 
# 
# processMethFiles <- function(fileNames){
#   
#   #SampleName <- fileNames[1]
#   SampleName <- fileNames
#   gene_meth<-data.table(hgnc_methylation[,SampleName])
#   gene_med<-hgnc_methylation[,c(1,(n+2))]
#   gene_meth_med<-cbind(gene_med,gene_meth)
#   colnames(gene_meth_med)<-c("gene","median_meth","meth")
#   gene_meth_med[gene_meth_med$median_meth==0,2]<-0.00000000001 ## to avoid divison by 0 for fold change
#   gene_meth_med[gene_meth_med$meth==0,3]<-0.00000000001 ## to get non zero fold change, removing fold change = 0 for down gene
#   gene_meth_med$fold_change<-abs(gene_meth_med$meth/gene_meth_med$median_meth)
#   ### Since fold change >2 or <.5 giving 10,000 genes among 40,000, I am taking 4 and .25 threshold,
#   # new threshold is giving around 5000 genes
#   upgene<-gene_meth_med[gene_meth_med$fold_change>=4,]## 271 genes
#   downgene<-gene_meth_med[gene_meth_med$fold_change<=0.25,]## 449 genes
#   pat_dm_genes <- data.table(rbind(upgene,downgene))[,1]
#   
#   
#   save(upgene,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Patient_DM_genes/UP_Gene/",str_sub(SampleName,1,16),"_upgene.Rda"))
#   save(downgene,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Patient_DM_genes/DOWN_Gene/",str_sub(SampleName,1,16),"_downgene.Rda"))
#   save(pat_dm_genes,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Patient_DM_genes/DE_Gene/",str_sub(SampleName,1,16),"_pat_dm_genes.Rda"))
#   
# }
# 
# print("cluster making started")
# cl <- makeCluster(mc <- getOption("cl.cores", 70),type = "FORK")
# invisible(clusterEvalQ(cl, library('glmnet')))
# invisible(clusterEvalQ(cl, library(janitor)))
# invisible(clusterEvalQ(cl, library(stringr)))
# invisible(clusterEvalQ(cl, library(plyr)))
# invisible(clusterEvalQ(cl, library(data.table)))
# invisible(clusterEvalQ(cl, library(readr)))
# invisible(clusterEvalQ(cl, library(gsubfn)))
# invisible(clusterEvalQ(cl, library(parallel)))
# invisible(clusterEvalQ(cl, library(sqldf)))
# invisible(clusterEvalQ(cl, library(reshape)))
# invisible(clusterEvalQ(cl, library(reshape2)))
# invisible(clusterEvalQ(cl, library(pbapply)))
# 
# print("cluster export started")
# parallel::clusterExport(cl=cl,c("hgnc_methylation"),envir=environment())#envir needed to be correct, export from global aswell as local
# print("cluster making done")
# 
# 
# fileNames<-colnames(hgnc_methylation[,-c(1,(n+2))])
# 
# parLapplyLB(cl,fileNames,processMethFiles)
# 
# stopCluster(cl)
# print("cluster stopped")