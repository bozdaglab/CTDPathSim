## THIS CODE IS FOR FINDING PATHWAYS ACTIVITY BASED CORRELATION FOR PATIENT and CELL LINE
## STAGE 4
## CODE 1
# Created By: Banabithi Bose
# Date Created: 4/28/2020

##### We will combine all the pathways for a pair and take union of DE or DM genes for correlation.
## THIS CODE CAN BE USED FOR EXPRESSION AND METHYLATION BASED SCORES

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
# library(Matrix)
# library(CCA)
# library(PMA)
# library(matrixcalc)


FindSim<- function(pat_de_genes,pat_reactome,expression,cell_de_genes,cell_reactome,ccle_expr){
    SampleName<-NULL
    cellLine_id<-NULL
    tryCatch(
    {
        #load(paste0("~/EXTENDED_CTDPathSim/",f,"/Sample_expression/",SampleName,".Rda", sep = ""))
        row_pat<-row.names(expression)
        #load("~/DRUG_PIPELINE1/CCLE_DATA/ccle_expr_rpkm.Rda")
        row_cell<-row.names(ccle_expr)
        row_pat_cell<-intersect(row_pat,row_cell)##common genes

        #cellLine_id<-str_sub(cell_reactome_files,1,-5)
        #print(cellLine_id)
        #cellLine_id<-"22RV1_PROSTATE"
        #load(paste0("~/PATIENT_CELL_PIPELINE3/CCLE_Pathway/Cell_reactome/",cellLine_id,".Rda"))## cell line's pathways
        union_reactome<-union(pat_reactome$ID,cell_reactome$ID)## union of pathways bet patient-cell line
        #load(paste0("~/PATIENT_CELL_PIPELINE3/CCLE_DE_genes/DE_Gene/",cellLine_id,"_cell_de_genes.Rda"))## cell line's DE genes
        patients_DE<-intersect(pat_de_genes$gene,row_pat_cell)
        cellLines_DE<-intersect(cell_de_genes$gene,row_pat_cell)
        union_de_genes1<-union(pat_de_genes$gene,cell_de_genes$gene)# union DE genes bet patient cell line
        union_de_genes<-intersect(union_de_genes1,row_pat_cell)

        Lx<-expression[union_de_genes,]## patient's expression with common DE genes
        r1<-data.table(row.names(Lx))
        Lx<-cbind(r1,Lx)
        #ccle_expr<-as.matrix(ccle_expr)
        #Ly<-data.frame(ccle_expr[union_de_genes,cellLine_id,drop=FALSE])## cell line's expression with union DE genes
        Ly<-as.data.frame(ccle_expr[union_de_genes,,drop=FALSE])## cell line's expression with union DE genes
        r2<-data.table(row.names(Ly))
        Ly<-cbind(r2,Ly)
        ### For reactome
        reactome_pathNames<-union_reactome
        #print(reactome_pathNames)
        if(length(reactome_pathNames)==0){
            reactome_pathNames<-"NA"
            reactome_descriptions<-"NA"
        }else{
            reactome_genes_by_term = data.frame(cbind(pathfindR.data::reactome_genes))
            reactome_term_descriptions = data.frame(pathfindR.data::reactome_descriptions)
            reactome_descriptions <- as.character(reactome_term_descriptions[reactome_pathNames,])
            reactome_pathway_genes<-data.table(unique(unlist(reactome_genes_by_term[reactome_pathNames,])))## taking union of genes in all intersected pathways
            colnames(reactome_pathway_genes)<-"V1"
            Common_DE_genes_pathway <-unique(intersect(reactome_pathway_genes$V1,union_de_genes))
        }
        if (nrow(reactome_pathway_genes)==0){
            reactomeLx<-matrix(1,nrow = 0,ncol=1)
        }else{

            colnames(reactome_pathway_genes)<-"V1"
            reactomeLx <- merge(reactome_pathway_genes,Lx,by="V1")
            reactomeLy <- merge(reactome_pathway_genes,Ly,by="V1")
        }
        if (nrow(reactomeLx)==0||nrow(reactomeLx)==1){

            print("corr not possible")

            reactome_spearScore<-"NA"

            reactome_expr_sim_scores<-cbind(paste0(str_sub(SampleName,1,-5)),cellLine_id,reactome_pathNames,reactome_descriptions,reactome_spearScore,length(patients_DE),length(cellLines_DE),
                                            length(unique(union_de_genes)),length(unique(reactome_pathway_genes$V1)),length(unique(Common_DE_genes_pathway)))
            reactome_expr_sim_scores<-data.table(reactome_expr_sim_scores)
            colnames(reactome_expr_sim_scores)<-c("patient","cellLine","reactome_pathway","description","spearman","pat_de","cell_de","union_de","reactome_genes","de_in_reactome")
            reactome_expr_sim_scores<-unique(reactome_expr_sim_scores[,c(1,2,5:10)])
            #return(reactome_expr_sim_scores)

        }
        if (ncol(expression)==3){
            reactome_spearCorr1<-cor(reactomeLy[,2],reactomeLx[,2],use ="complete",method=c("spearman"))
            reactome_spearCorr2<-cor(reactomeLy[,2],reactomeLx[,3],use ="complete",method=c("spearman"))
            reactome_spearCorr3<-cor(reactomeLy[,2],reactomeLx[,4],use ="complete",method=c("spearman"))
            x1<-data.frame(rbind(reactome_spearCorr1[1],reactome_spearCorr2[1],reactome_spearCorr3[1]))
            colnames(x1)<-"V1"
        }
        if (ncol(expression)==4){
            reactome_spearCorr1<-cor(reactomeLy[,2],reactomeLx[,2],use ="complete",method=c("spearman"))
            reactome_spearCorr2<-cor(reactomeLy[,2],reactomeLx[,3],use ="complete",method=c("spearman"))
            reactome_spearCorr3<-cor(reactomeLy[,2],reactomeLx[,4],use ="complete",method=c("spearman"))
            reactome_spearCorr4<-cor(reactomeLy[,2],reactomeLx[,5],use ="complete",method=c("spearman"))
            x1<-data.frame(rbind(reactome_spearCorr1[1],reactome_spearCorr2[1],reactome_spearCorr3[1],
                                reactome_spearCorr4[1]))
            colnames(x1)<-"V1"
        }
        if (ncol(expression)==5){
            reactome_spearCorr1<-cor(reactomeLy[,2],reactomeLx[,2],use ="complete",method=c("spearman"))
            reactome_spearCorr2<-cor(reactomeLy[,2],reactomeLx[,3],use ="complete",method=c("spearman"))
            reactome_spearCorr3<-cor(reactomeLy[,2],reactomeLx[,4],use ="complete",method=c("spearman"))
            reactome_spearCorr4<-cor(reactomeLy[,2],reactomeLx[,5],use ="complete",method=c("spearman"))
            reactome_spearCorr5<-cor(reactomeLy[,2],reactomeLx[,6],use ="complete",method=c("spearman"))
            x1<-data.frame(rbind(reactome_spearCorr1[1],reactome_spearCorr2[1],reactome_spearCorr3[1]),reactome_spearCorr4[1],reactome_spearCorr5[1])
            colnames(x1)<-"V1"
        }

        if(abs(mean(x1$V1,na.rm=TRUE))=="NaN"){
            reactome_spearScore<-"NA"
        }else{
            reactome_spearScore<-abs(mean(x1$V1,na.rm=TRUE))
        }

        reactome_expr_sim_scores<-cbind(reactome_spearScore,length(patients_DE),length(cellLines_DE),
                                        length(unique(union_de_genes)),length(unique(reactome_pathway_genes$V1)),length(unique(Common_DE_genes_pathway)))
        reactome_expr_sim_scores<-data.table(reactome_expr_sim_scores)
        colnames(reactome_expr_sim_scores)<-c("spearman","pat_de","cell_de","union_de","reactome_genes","de_in_reactome")
    },
    error = function(error_condition) {
        cat(paste0(error_condition),file=paste0("~/EXTENDED_CTDPathSim/FindCorr.Error.txt",sep="\n",append=TRUE))
        return(reactome_expr_sim_scores)
    },finally={
        return(reactome_expr_sim_scores)
    })
}

# ##RUNNING THE FUNCTION
# load("~/CTDPathSim/inst/extdata/Example.patient1.de.genes.rda")
# load("~/CTDPathSim/inst/extdata/Example.patient1.pathways.rda")
# load("~/CTDPathSim/inst/extdata/Example.patient1.deconvoluted.expression.rda")
# load("~/CTDPathSim/inst/extdata/Example.cellLine1.de.genes.rda")
# load("~/CTDPathSim/inst/extdata/Example.cellLine1.pathways.rda")
# load("~/CTDPathSim/inst/extdata/Example.cellLine1.expression.rda")
#
# patient1.cellLine1.spearman<-FindCorr(pat_de_genes,pat_reactome,expression,cell_de_genes,cell_reactome,ccle_expr)
#
# ## RESULT IS SAVED IN CTDPathSim_Bioconductor/patient1.cellLine1.spearman.rda
