#### STAGE 3 
## CODE 2
#### Getting Cell line specific DE/DM genes  (This code is same for all cancer types, so needs to run it once####
# Created By: Banabithi Bose
# Date Created: 4/28/2020
#### Use CCLE expression RPKM to get median and fold change #####

#Compute median expression of a gene across all samples
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
library(matrixStats)
load("~/DRUG_PIPELINE1/CCLE_DATA/ccle_expr_rpkm.Rda")#RPKM i.e. FPKM
length(unique(row.names(ccle_expr)))#54271## too many genes, so will take genes which are coomon to our TCGA cohort
r1<-data.table(row.names(ccle_expr))
ccle_expr<-data.table(cbind(r1,ccle_expr))

load("~/DRUG_PIPELINE1/DATASETS/TCGA_HGNC/BRCA.Rda")#FPKM#1103 patients
genes<-data.table(unique(intersect(ccle_expr$V1,RnaSeq_fpkm_hgnc_data$Group.1)))#30681
ccle_expr<-merge(genes,ccle_expr,by="V1")

c1<-data.table(colnames(ccle_expr[,-1]))
c2<-data.table(str_split_fixed(c1$V1,"_",2))
c3<-cbind(c1,c2[,2])
i=length(unique(c3$V2))

for (i in 1:26){
  #i=1
  cell<-c3[c3$V2==unique(c3$V2)[i],]$V1
  ccle_expr_rpkm<-cbind(ccle_expr[,1],ccle_expr[,..cell])
ccle_expr_rpkm$Median_expr <- rowMedians(as.matrix(ccle_expr_rpkm[,-1]),na.rm=T)


processMethFiles <- function(fileNames){

  #SampleName <- fileNames[3]
  SampleName <- fileNames
  gene_expr<-data.table(ccle_expr_rpkm[,..SampleName])
  gene_med<-ccle_expr_rpkm[,c("V1","Median_expr")]
  gene_expr_med<-cbind(gene_med,gene_expr)
  colnames(gene_expr_med)<-c("gene","median_expr","expr")
  gene_expr_med[gene_expr_med$median_expr==0,2]<-0.00000000001 ## to avoid divison by 0 for fold change
  gene_expr_med[gene_expr_med$expr==0,3]<-0.00000000001 ## to get non zero fold change, removing fold change = 0 for down gene
  gene_expr_med$fold_change<-gene_expr_med$expr/gene_expr_med$median_expr
  ### Since fold change >2 or <.5 giving 10,000 genes among 40,000, I am taking 4 and .25 threshold,
  # new threshold is giving around 5000 genes
  upgene<-unique(gene_expr_med[gene_expr_med$fold_change>=4,])## 2237 genes
  downgene<-unique(gene_expr_med[gene_expr_med$fold_change<=0.25,])## 2582 genes
  cell_de_genes <- data.table(rbind(upgene,downgene))[,1]


  save(upgene,file=paste0("/users/home/bbose/PATIENT_CELL_PIPELINE3/CCLE_DE_genes/UP_Gene/",SampleName,"_upgene.Rda"))
  save(downgene,file=paste0("/users/home/bbose/PATIENT_CELL_PIPELINE3/CCLE_DE_genes/DOWN_Gene/",SampleName,"_downgene.Rda"))
  save(cell_de_genes,file=paste0("/users/home/bbose/PATIENT_CELL_PIPELINE3/CCLE_DE_genes/DE_Gene/",SampleName,"_cell_de_genes.Rda"))

}

print("cluster making started")
cl <- makeCluster(mc <- getOption("cl.cores", 60),type = "FORK")
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
parallel::clusterExport(cl=cl,c("ccle_expr_rpkm"),envir=environment())#envir needed to be correct, export from global aswell as local
print("cluster making done")


fileNames<-colnames(ccle_expr_rpkm[,-c("V1","Median_expr")])

parLapplyLB(cl,fileNames,processMethFiles)

stopCluster(cl)
print("cluster stopped")

}
# After this we need to find pathways for each cell line


################## ##### Use CCLE Methylation to get median and fold change #####
#Compute median expression of a gene across all samples 
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
library(matrixStats)
load("~/DRUG_PIPELINE1/CCLE_DATA/imputed_ccle_meth.Rda")
ccle_meth[ccle_meth==0]<-0.00000000001 ## to avoid zero in M transformation
#Writing a function for transforming beta values to M values
Mtranform <- function(x) {log(x[1]/(1-x[1]))}

ccle_meth[, 1:831] <- as.data.frame(lapply(ccle_meth[, 1:831], FUN = function(x) {sapply(x, FUN = Mtranform)}))


ccle_methylation<-na.omit(ccle_meth)#5000
gene<-data.table(row.names(ccle_methylation))
ccle_methylation1<-data.table(cbind(gene,ccle_methylation))
# we will take common genes with TCGA cohort for methylation
load("~/Ajani_PIPELINE2_BREAST_LASSO/PIPELINE2_DATAFRAMES/All_methylation.Rda")
x<-data.table(row.names(All_methylation))
load("~/PATIENT_CELL_PIPELINE1/DATAFRAMES/Ensemble_hgnc.Rda")
All_methylation<-cbind(x,All_methylation)
hgnc_methylation <- merge(Ensemble_hgnc,All_methylation,by.x = "ensembl_gene_id",by.y = "V1")
hgnc_methylation<-hgnc_methylation[,-1]
genes<-data.table(intersect(hgnc_methylation$hgnc_symbol,gene$V1))#10251
ccle_methylation1<-merge(genes,ccle_methylation1,by="V1")
c1<-data.table(colnames(ccle_meth))
c2<-data.table(str_split_fixed(c1$V1,"_",2))
c3<-cbind(c1,c2[,2])
i=length(unique(c3$V2))

for (i in 1:24){
  #i=1
  cell<-c3[c3$V2==unique(c3$V2)[i],]$V1
  ccle_methylation<-cbind(ccle_methylation1[,1],ccle_methylation1[,..cell])
  ccle_methylation$Median_expr <- rowMedians(as.matrix(ccle_methylation[,-1]),na.rm=T)


  processMethFiles <- function(fileNames){

    #SampleName <- fileNames[3]
    SampleName <- fileNames
    gene_meth<-data.table(ccle_methylation[,..SampleName])
    gene_med<-ccle_methylation[,c("V1","Median_expr")]
    gene_meth_med<-cbind(gene_med,gene_meth)
    colnames(gene_meth_med)<-c("gene","median_meth","meth")
    gene_meth_med[gene_meth_med$median_meth==0,2]<-0.00000000001 ## to avoid divison by 0 for fold change
    gene_meth_med[gene_meth_med$meth==0,3]<-0.00000000001 ## to get non zero fold change, removing fold change = 0 for down gene
    gene_meth_med$fold_change<-gene_meth_med$meth/gene_meth_med$median_meth
    ### Since fold change >2 or <.5 giving 10,000 genes among 40,000, I am taking 4 and .25 threshold,
    # new threshold is giving around 5000 genes
    upgene<-gene_meth_med[gene_meth_med$fold_change>=4,]## 643 genes
    downgene<-gene_meth_med[gene_meth_med$fold_change<=0.25,]## 1925 genes
    cell_dm_genes <- data.table(rbind(upgene,downgene))[,1]


    save(upgene,file=paste0("/users/home/bbose/PATIENT_CELL_PIPELINE3/CCLE_DM_genes/UP_Gene/",SampleName,"_upgene.Rda"))
    save(downgene,file=paste0("/users/home/bbose/PATIENT_CELL_PIPELINE3/CCLE_DM_genes/DOWN_Gene/",SampleName,"_downgene.Rda"))
    save(cell_dm_genes,file=paste0("/users/home/bbose/PATIENT_CELL_PIPELINE3/CCLE_DM_genes/DE_Gene/",SampleName,"_cell_dm_genes.Rda"))

  }

  print("cluster making started")
  cl <- makeCluster(mc <- getOption("cl.cores", 60),type = "FORK")
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
  parallel::clusterExport(cl=cl,c("ccle_methylation"),envir=environment())#envir needed to be correct, export from global aswell as local
  print("cluster making done")


  fileNames<-colnames(ccle_methylation[,-c("V1","Median_expr")])

  parLapplyLB(cl,fileNames,processMethFiles)

  stopCluster(cl)
  print("cluster stopped")

}
## After this we need to find pathways for each cell line

