##### This is for pathways with cancer genes from patient specific DE genes  #####
## THIS CODE IS FOR FINDING PATHWAYS FOR PATIENTS DE/DM MUTATED GENES
## STAGE 3 
## CODE 3
# Created By: Banabithi Bose
# Date Created: 4/28/2020
##### This is for pathways with cancer genes from patient specific DE genes  #####
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
################# Mutated cancer genes ###################
load("/users/home/bbose/PATIENT_CELL_PIPELINE1/DATAFRAMES/expressionCancerGeneList.Rda")#703
#KEGG_universe <- data.table(unique(pathfindR::kegg_descriptions))
#universe=327

processMethFiles <- function(fileNames){

  #SampleName <- fileNames[1]

  SampleName <- fileNames
  load(paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Patient_DE_genes/DE_Gene/",SampleName))
  common_gene_mut_de <-merge(pat_de_genes,CancerGeneList,by.x = "gene",by.y = "V1")

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
      print("no patient kegg")

    }else{
      d1<-data.frame(ekegg$Term_Description)
      d2<-data.frame(ekegg$adj_p)
      d3<-data.frame(ekegg$ID)
      pat_kegg<-cbind(d3,d1,d2)
      colnames(pat_kegg)<-c("ID","kegg_pathway","p.adjust.pat")
      save(pat_kegg,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Patients_Pathway/Pat_kegg_expr/",str_sub(SampleName,1,28),".Rda"))
    }

      if (is.null(nrow(reactome))== TRUE){
        print("no patient reactome")

      }else{
        d1<-data.frame(reactome$Term_Description)
        d2<-data.frame(reactome$adj_p)
        d3<-data.frame(reactome$ID)
        pat_reactome<-cbind(d3,d1,d2)
        colnames(pat_reactome)<-c("ID","reactome_pathway","p.adjust.pat")
        save(pat_reactome,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Patients_Pathway/Pat_reactome_expr/",str_sub(SampleName,1,28),".Rda"))

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

fileNames<-list.files("/users/home/bbose/PATIENT_CELL_PIPELINE3/OVARIAN/Patient_DE_genes/DE_Gene/", pattern="*.Rda")#1102

parLapplyLB(cl,fileNames,processMethFiles)

stopCluster(cl)
print("cluster stopped")





## THIS CODE IS FOR FINDING PATHWAYS FOR PATIENTS DE/DM MUTATED GENES
## STAGE 3 
## CODE 3
# Created By: Banabithi Bose
# Date Created: 4/28/2020
##### This is for pathways with cancer genes from patient specific DM genes  #####

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
################# Mutated cancer genes ###################
load("/users/home/bbose/PATIENT_CELL_PIPELINE1/DATAFRAMES/expressionCancerGeneList.Rda")#703
#KEGG_universe <- data.table(unique(pathfindR::kegg_descriptions))
#universe=327

processMethFiles <- function(fileNames){
  
  #SampleName <- fileNames[1]
  
  SampleName <- fileNames
  load(paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Patient_DM_genes/DE_Gene/",SampleName))
  common_gene_mut_de <-merge(pat_dm_genes,CancerGeneList,by.x = "gene",by.y = "V1")
  
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
    print("no patient kegg")
    
  }else{
    d1<-data.frame(ekegg$Term_Description)
    d2<-data.frame(ekegg$adj_p)
    d3<-data.frame(ekegg$ID)
    pat_kegg<-cbind(d3,d1,d2)
    colnames(pat_kegg)<-c("ID","kegg_pathway","p.adjust.pat")
    save(pat_kegg,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Patients_Pathway/Pat_kegg_meth/",str_sub(SampleName,1,28),".Rda"))
  }
  
  if (is.null(nrow(reactome))== TRUE){
    print("no patient reactome")
    
  }else{
    d1<-data.frame(reactome$Term_Description)
    d2<-data.frame(reactome$adj_p)
    d3<-data.frame(reactome$ID)
    pat_reactome<-cbind(d3,d1,d2)
    colnames(pat_reactome)<-c("ID","reactome_pathway","p.adjust.pat")
    save(pat_reactome,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Patients_Pathway/Pat_reactome_meth/",str_sub(SampleName,1,28),".Rda"))
    
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


fileNames<-list.files("~/PATIENT_CELL_PIPELINE3/OVARIAN/Patient_DM_genes/DE_Gene/", pattern="*.Rda")#613


parLapplyLB(cl,fileNames,processMethFiles)

stopCluster(cl)
print("cluster stopped")





