##### This is for pathways with cancer genes from cell lines specific DE/DM genes  #####
# needs to run once for all cancer types
## STAGE 3 
## CODE 4

# Created By: Banabithi Bose
# Date Created: 4/28/2020

.libPaths(c("/users/home/bbose/R/x86_64-redhat-linux-gnu-library/3.5","/users/home/bbose/R/x86_64-redhat-linux-gnu-library/3.6","/usr/lib64/R/library","/usr/share/R/library"))

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
################# Mutated cancer genes ###################
load("/users/home/bbose/PATIENT_CELL_PIPELINE1/DATAFRAMES/expressionCancerGeneList.Rda")#703
#KEGG_universe <- data.table(unique(pathfindR::kegg_descriptions))
#universe=327

processMethFiles <- function(fileNames){
  
  #SampleName <- fileNames[1]
  
  SampleName <- fileNames
  load(paste0("~/PATIENT_CELL_PIPELINE3/CCLE_DE_genes/DE_Gene/",SampleName))
  common_gene_mut_de <-merge(cell_de_genes,CancerGeneList,by.x = "gene",by.y = "V1")
  
  #### Kegg Pathways with cell lines genes
  ekegg<-enrichment(
    input_genes=common_gene_mut_de$gene,
    genes_by_term = pathfindR::kegg_genes,
    term_descriptions = pathfindR::kegg_descriptions,
    adj_method = "bonferroni",
    enrichment_threshold = 0.05,
    sig_genes_vec=common_gene_mut_de$gene,
    background_genes=unlist(kegg_genes))
  
  reactome<-enrichment(
    input_genes=common_gene_mut_de$gene,
    genes_by_term = pathfindR::reactome_genes,
    term_descriptions = pathfindR::reactome_descriptions,
    adj_method = "bonferroni",
    enrichment_threshold = 0.05,
    sig_genes_vec=common_gene_mut_de$gene,
    background_genes=unlist(reactome_genes))
  
    if (is.null(nrow(ekegg))== TRUE){
      print("no cell line kegg")
      
    }else{
      d1<-data.frame(ekegg$Term_Description)
      d2<-data.frame(ekegg$adj_p)
      d3<-data.frame(ekegg$ID)
      cell_kegg<-cbind(d3,d1,d2)
      colnames(cell_kegg)<-c("ID","kegg_pathway","p.adjust.cell")
      save(cell_kegg,file=paste0("/users/home/bbose/PATIENT_CELL_PIPELINE3/CCLE_Pathway/Cell_kegg/",str_sub(SampleName,1,-19),".Rda"))
    }
      
      if (is.null(nrow(reactome))== TRUE){
        print("no cell line reactome")
        
      }else{
        d1<-data.frame(reactome$Term_Description)
        d2<-data.frame(reactome$adj_p)
        d3<-data.frame(reactome$ID)
        cell_reactome<-cbind(d3,d1,d2)
        colnames(cell_reactome)<-c("ID","reactome_pathway","p.adjust.cell")
        save(cell_reactome,file=paste0("/users/home/bbose/PATIENT_CELL_PIPELINE3/CCLE_Pathway/Cell_reactome/",str_sub(SampleName,1,-19),".Rda"))
        
      }
  
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
parallel::clusterExport(cl=cl,c("CancerGeneList"),envir=environment())#envir needed to be correct, export from global aswell as local
print("cluster making done")

fileNames<-list.files("/users/home/bbose/PATIENT_CELL_PIPELINE3/CCLE_DE_genes/DE_Gene/", pattern="*.Rda")#1019

parLapplyLB(cl,fileNames,processMethFiles)

stopCluster(cl)
print("cluster stopped")




