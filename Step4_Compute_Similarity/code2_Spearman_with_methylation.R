
## THIS CODE IS FOR FINDING PATHWAYS ACTIVITY BASED CORRELATION FOR PATIENTS DM GENES IN PATHWAY GENES
## STAGE 4
## CODE 2
# Created By: Banabithi Bose
# Date Created: 4/28/2020

### FOR UNION 
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
library(Matrix)
library(CCA)
library(PMA)
library(matrixcalc)

PC_Correlation<-data.frame("","","","","","","","","","","","","")
colnames(PC_Correlation)<-c("patient","cellLine","reactome_pathway","description","eucledian","spearman","pearson","cca","pat_dm","cell_dm","union_dm","reactome_genes","dm_in_reactome")

processCorr<- function(cell_reactome_files,pat_reactome,pat_dm_genes,SampleName){
  
  tryCatch( 
    {
      load(paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Sample_methylation/",SampleName,".Rda", sep = ""))# patient's methylation
      load("~/DRUG_PIPELINE1/CCLE_DATA/imputed_ccle_meth.Rda")
      #cellLine_id<-cell_reactome_files[1]
      cellLine_id<-cell_reactome_files
      print(cellLine_id)
      #cellLine_id<-"22RV1_PROSTATE"
      load(paste0("~/PATIENT_CELL_PIPELINE3/CCLE_Pathway/Cell_reactome/",cellLine_id,".Rda"))## cell line's pathways
      union_reactome<-union(pat_reactome$ID,cell_reactome$ID)## union pathways bet patient-cell line
      load(paste0("~/PATIENT_CELL_PIPELINE2/CCLE_DM_genes/DE_Gene/",cellLine_id,"_cell_dm_genes.Rda"))## cell line's DM genes
      union_dm_genes<-data.table(union(pat_dm_genes$gene,cell_dm_genes$gene))# union DM genes bet patient cell line
      ## need to check how many union_dm_genes are
      #m<-data.table(row.names(methylation))
      #union_dm_genes1<-merge(m,union_dm_genes,by="V1")
      
      
      row_pat<-row.names(methylation)
      row_cell<-row.names(ccle_meth)
      row_pat_cell<-intersect(row_pat,row_cell)
      patients_DM<-intersect(pat_dm_genes$gene,row_pat_cell)
      cellLines_DM<-intersect(cell_dm_genes$gene,row_pat_cell)
      union_dm_genes1<-union(pat_dm_genes$gene,cell_dm_genes$gene)# union DE genes bet patient cell line
      union_dm_genes<-intersect(union_dm_genes1,row_pat_cell)
      
      Lx<-methylation[union_dm_genes,]## patient's expression with common DE genes
      r1<-data.table(row.names(Lx))
      Lx<-cbind(r1,Lx)
      ccle_meth<-as.matrix(ccle_meth)
      colnames(ccle_meth)<-gsub("^X","",colnames(ccle_meth))
      Ly<-data.frame(ccle_meth[union_dm_genes,cellLine_id,drop=F])## cell line's methylation with common DE genes
      r2<-data.table(row.names(Ly))
      Ly<-cbind(r2,Ly)
      ### For reactome
      reactome_pathNames<-union_reactome
      print(reactome_pathNames)
      if(length(reactome_pathNames)==0){
        reactome_pathNames<-"NA"
        reactome_descriptions<-"NA"
      }
      reactome_genes_by_term = data.frame(cbind(pathfindR::reactome_genes))
      reactome_term_descriptions = data.frame(pathfindR::reactome_descriptions)
      reactome_descriptions <- as.character(reactome_term_descriptions[reactome_pathNames,])
      reactome_pathway_genes<-data.table(unique(unlist(reactome_genes_by_term[reactome_pathNames,])))## taking union of genes in all intersected pathways
      Common_DM_genes_pathway <-unique(intersect(reactome_pathway_genes$V1,union_dm_genes))
      if (nrow(reactome_pathway_genes)==0){
        reactomeLx<-matrix(1,nrow = 0,ncol=1)
      }else{
        
        colnames(reactome_pathway_genes)<-"V1"
        reactomeLx <- merge(reactome_pathway_genes,Lx,by="V1")
        reactomeLy <- merge(reactome_pathway_genes,Ly,by="V1")
      }
      if (nrow(reactomeLx)==0||nrow(reactomeLx)==1){
        
        print("corr not possible")
        reactome_euc_Score<-"NA"
        reactome_pearScore<-"NA"
        reactome_spearScore<-"NA"
        reactome_ccaScore<-"NA"
        reactome_meth_sim_scores<-cbind(SampleName,cellLine_id,reactome_pathNames,reactome_descriptions,reactome_euc_Score,reactome_spearScore,reactome_pearScore,reactome_ccaScore,length(patients_DM),length(cellLines_DM),
                                    length(unique(union_dm_genes)),length(unique(reactome_pathway_genes$V1)),length(unique(Common_DM_genes_pathway)))
        reactome_meth_sim_scores<-data.table(reactome_meth_sim_scores)
        colnames(reactome_meth_sim_scores)<-c("patient","cellLine","reactome_pathway","description","eucledian","spearman","pearson","cca","pat_dm","cell_dm","union_dm","reactome_genes","dm_in_reactome")
        return(reactome_meth_sim_scores)
        
      }else{
        if(ncol(Lx)==5){
        ### Euclidean distance 
        reactome_euc1<-dist(rbind(unlist(reactomeLy[,2]),unlist(reactomeLx[,2])))
        reactome_euc2<-dist(rbind(unlist(reactomeLy[,2]),unlist(reactomeLx[,3])))
        reactome_euc3<-dist(rbind(unlist(reactomeLy[,2]),unlist(reactomeLx[,4])))
        reactome_euc4<-dist(rbind(unlist(reactomeLy[,2]),unlist(reactomeLx[,5])))
        reactome_euc<-data.frame(rbind(abs(reactome_euc1[1]),abs(reactome_euc2[1]),abs(reactome_euc3[1]),reactome_euc4[1]))
        colnames(reactome_euc)<-"V1"
        reactome_euc_Score<-mean(reactome_euc$V1,na.rm=TRUE)
        }else{
          reactome_euc1<-dist(rbind(unlist(reactomeLy[,2]),unlist(reactomeLx[,2])))
          reactome_euc2<-dist(rbind(unlist(reactomeLy[,2]),unlist(reactomeLx[,3])))
          reactome_euc3<-dist(rbind(unlist(reactomeLy[,2]),unlist(reactomeLx[,4])))
          #reactome_euc4<-dist(rbind(unlist(reactomeLy[,2]),unlist(reactomeLx[,5])))
          reactome_euc<-data.frame(rbind(abs(reactome_euc1[1]),abs(reactome_euc2[1]),abs(reactome_euc3[1])))
          colnames(reactome_euc)<-"V1"
          reactome_euc_Score<-mean(reactome_euc$V1,na.rm=TRUE)
        }
        ### Pearson correlation and Spearman
        if (nrow(reactomeLx)==1||nrow(na.omit(reactomeLy))==0||nrow(na.omit(reactomeLx))==0){
          reactome_euc_Score<-"NA"
          reactome_pearScore<-"NA"
          reactome_spearScore<-"NA"
          reactome_ccaScore<-"NA"
        }else{
          if(ncol(Lx)==5){
          reactome_pearCorr1<-cor(reactomeLy[,2],reactomeLx[,2],use ="complete",method=c("pearson"))
          reactome_pearCorr2<-cor(reactomeLy[,2],reactomeLx[,3],use ="complete",method=c("pearson"))
          reactome_pearCorr3<-cor(reactomeLy[,2],reactomeLx[,4],use ="complete",method=c("pearson"))
          reactome_pearCorr4<-cor(reactomeLy[,2],reactomeLx[,5],use ="complete",method=c("pearson"))
          x<-data.frame(rbind(reactome_pearCorr1[1],reactome_pearCorr2[1],reactome_pearCorr3[1],reactome_pearCorr4[1]))
          colnames(x)<-"V1"
          reactome_spearCorr1<-cor(reactomeLy[,2],reactomeLx[,2],use ="complete",method=c("spearman"))
          reactome_spearCorr2<-cor(reactomeLy[,2],reactomeLx[,3],use ="complete",method=c("spearman"))
          reactome_spearCorr3<-cor(reactomeLy[,2],reactomeLx[,4],use ="complete",method=c("spearman"))
          reactome_spearCorr4<-cor(reactomeLy[,2],reactomeLx[,5],use ="complete",method=c("spearman"))
          x1<-data.frame(rbind(reactome_spearCorr1[1],reactome_spearCorr2[1],reactome_spearCorr3[1],reactome_spearCorr4[1]))
          colnames(x1)<-"V1"
          reactome_ccaScore<-"NA"
          }else{
            reactome_pearCorr1<-cor(reactomeLy[,2],reactomeLx[,2],use ="complete",method=c("pearson"))
            reactome_pearCorr2<-cor(reactomeLy[,2],reactomeLx[,3],use ="complete",method=c("pearson"))
            reactome_pearCorr3<-cor(reactomeLy[,2],reactomeLx[,4],use ="complete",method=c("pearson"))
            #reactome_pearCorr4<-cor(reactomeLy[,2],reactomeLx[,5],use ="complete",method=c("pearson"))
            x<-data.frame(rbind(reactome_pearCorr1[1],reactome_pearCorr2[1],reactome_pearCorr3[1]))
            colnames(x)<-"V1"
            reactome_spearCorr1<-cor(reactomeLy[,2],reactomeLx[,2],use ="complete",method=c("spearman"))
            reactome_spearCorr2<-cor(reactomeLy[,2],reactomeLx[,3],use ="complete",method=c("spearman"))
            reactome_spearCorr3<-cor(reactomeLy[,2],reactomeLx[,4],use ="complete",method=c("spearman"))
            #reactome_spearCorr4<-cor(reactomeLy[,2],reactomeLx[,5],use ="complete",method=c("spearman"))
            x1<-data.frame(rbind(reactome_spearCorr1[1],reactome_spearCorr2[1],reactome_spearCorr3[1]))
            colnames(x1)<-"V1"
            reactome_ccaScore<-"NA"
          }
          if(abs(mean(x$V1,na.rm=TRUE))=="NaN"){
            
            reactome_pearScore<-"NA"
            reactome_spearScore<-"NA"
            reactome_ccaScore="NA"
          }else{
            reactome_pearScore<-abs(mean(x$V1,na.rm=TRUE))
            reactome_spearScore<-abs(mean(x1$V1,na.rm=TRUE))
            
            tryCatch(
              {correlcca <- cc(reactomeLx[,-1],reactomeLy[,-1,drop=F])
              reactome_ccaScore<-abs(correlcca$cor)
              },error = function(error_condition) {
                
                reactome_ccaScore<-"NA"
                #cat(paste0(SampleName,cellLine_id,intersect_reactome," 2 "),file="/users/home/bbose/PATIENT_CELL_PIPELINE3/ERROR_FILES/reactome_Cell_CellTypeCorrmeth.txt",sep="\n",append=TRUE)
                
              },finally={
                
              }
            )
          }
        }
        reactome_meth_sim_scores<-cbind(SampleName,cellLine_id,reactome_pathNames,reactome_descriptions,reactome_euc_Score,reactome_spearScore,reactome_pearScore,reactome_ccaScore,length(patients_DM),length(cellLines_DM),
                                    length(unique(union_dm_genes)),length(unique(reactome_pathway_genes$V1)),length(unique(Common_DM_genes_pathway)))
        reactome_meth_sim_scores<-data.table(reactome_meth_sim_scores)
        colnames(reactome_meth_sim_scores)<-c("patient","cellLine","reactome_pathway","description","eucledian","spearman","pearson","cca","pat_dm","cell_dm","union_dm","reactome_genes","dm_in_reactome")
        
      }
      
    },
    error = function(error_condition) {
      cat(paste0(error_condition),file="/users/home/bbose/PATIENT_CELL_PIPELINE3/ERROR_FILES/Sim1Reacmeth",sep="\n",append=TRUE)
      
    }, 
    
    finally={
      return(reactome_meth_sim_scores)
    }
    
  )
  return(reactome_meth_sim_scores) 
  
}

processMethFiles <- function(fileNames,processCorr){
  
  #SampleName <- str_sub(fileNames[1],end=-5)
  SampleName <- str_sub(fileNames,end=-5)
  print(SampleName)
  load(paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Patients_Pathway/Pat_reactome_meth/",SampleName,"_pat_dm_gene.Rda"))# patient's pathways
  load(paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Patient_DM_genes/DE_Gene/",SampleName,"_pat_dm_genes.Rda")) # patient's DM genes
  cell_files<-list.files("~/PATIENT_CELL_PIPELINE2/CCLE_DM_genes/DE_Gene/", pattern="*.Rda")#831 cell lines with DE genes
  cell_reactome_files1<-list.files("~/PATIENT_CELL_PIPELINE3/CCLE_Pathway/Cell_reactome/", pattern="*.Rda")#1018 cell lines with reactome pathways
  cells_meth1<-str_sub(cell_reactome_files1,1,-5)
  cells_meth2<-str_sub(cell_files,1,-19)
  common_cells<-intersect(cells_meth1,cells_meth2)#823 (cell lines we got from expression based pathways)
  common_cells<-gsub("^X","",common_cells)
  load("~/DRUG_PIPELINE1/CCLE_DATA/imputed_ccle_meth.Rda")
  ccle_meth<-as.matrix(ccle_meth)
  colnames(ccle_meth)<-gsub("^X","",colnames(ccle_meth))
  ccle_meth<-ccle_meth[,common_cells]
  #cell_reactome_files<-cell_reactome_files[1:10]
  cell_reactome_files<-colnames(ccle_meth)
  Corr_2<-lapply(cell_reactome_files,processCorr,pat_reactome,pat_dm_genes,SampleName)
  PC_Correlation1<-rbindlist(Corr_2)
  PC_Correlation<-rbind(PC_Correlation,PC_Correlation1)[-1,]
  #save(PC_Correlation,file=paste0("/users/home/bbose/PATIENT_CELL_PIPELINE3/SIMILARITIES/reactome_Meth_Sim/",SampleName,".Rda"))
  
  return(PC_Correlation)
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
invisible(clusterEvalQ(cl, library(matrixcalc)))
print("cluster export started")
parallel::clusterExport(cl=cl,c("PC_Correlation","processCorr"),envir=environment())#envir needed to be correct, export from global aswell as local
print("cluster making done")

fileNames1<-gsub("(_pat_dm_gene)","",unique(list.files("~/PATIENT_CELL_PIPELINE3/OVARIAN/Patients_Pathway/Pat_reactome_meth/", pattern="*.Rda")))#594#We will use the pathways we got from expression in this case also
fileNames2<-unique(list.files("~/PATIENT_CELL_PIPELINE3/OVARIAN/Sample_methylation/", pattern="*.Rda"))#623 patients
fileNames<-intersect(fileNames1,fileNames2)[201:594]#594, patients have sample methylation as well as reactome pathways

Patient_Cell_Reactome_Meth_Sim2<-rbindlist(parLapplyLB(cl,fileNames,processMethFiles,processCorr))
save(Patient_Cell_Reactome_Meth_Sim2,file="~/PATIENT_CELL_PIPELINE3/OVARIAN/SIMILARITIES1/Patient_Cell_Reactome_Meth_Sim2.Rda")
stopCluster(cl)
print("cluster stopped")
