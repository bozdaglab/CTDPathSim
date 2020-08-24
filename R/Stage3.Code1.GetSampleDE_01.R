## STAGE 3
## CODE 1
#### Getting Patient specific DE/DM genes  ####
# Created By: Banabithi Bose
# Date Created: 4/27/2020
#
# #### Use FPKM data (before deconvolution) to get median and fold change #####
# ##Compute median expression of a gene across all samples
# .libPaths(c("/users/home/bbose/R/x86_64-redhat-linux-gnu-library/3.5","/users/home/bbose/R/x86_64-redhat-linux-gnu-library/3.6","/usr/lib64/R/library","/usr/share/R/library"))
# #.libPaths()
# gc(reset=TRUE)
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
# library(matrixStats)
##THIS FUCTION CAN BE USED FOR SAMPLE AND CELL LINE BOTH
GetSampleDE<-function(RnaSeq_fpkm_hgnc_data,parallel=c('TRUE','FALSE'),ncores=2){


    processFiles <- function(fileNames){

        ##MATT TAKE CARE OF THE FOLDER
        dir.create("~/Patient_DE_genes/")
        dir.create("~/Patient_DE_genes/UP_Gene/")
        dir.create("~/Patient_DE_genes/DOWN_Gene/")
        dir.create("~/Patient_DE_genes/DE_Gene/")
        #SampleName <- fileNames[1]
        SampleName <- fileNames
        RnaSeq_fpkm_hgnc_data<-as.data.frame(RnaSeq_fpkm_hgnc_data)
        gene_expr<-RnaSeq_fpkm_hgnc_data[, SampleName,drop=FALSE]
        gene_med<-RnaSeq_fpkm_hgnc_data[,c(1,n),drop=FALSE]
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



        save(upgene,file=paste0("~/Patient_DE_genes/UP_Gene/",SampleName,"_upgene.Rda"))
        save(downgene,file=paste0("~/Patient_DE_genes/DOWN_Gene/",SampleName,"_downgene.Rda"))
        save(pat_de_genes,file=paste0("~/Patient_DE_genes/DE_Gene/",SampleName,"_pat_de_genes.Rda"))

    }


    x <- data.table(rowMedians(RnaSeq_fpkm_hgnc_data,na.rm=TRUE))#37109 genes total in hgnc
    colnames(x)<-"Median_expr"
    RnaSeq_fpkm_hgnc_data<-cbind(RnaSeq_fpkm_hgnc_data,x)
    n<-ncol(RnaSeq_fpkm_hgnc_data)
    n<-as.integer(n)
    fileNames1<-colnames(RnaSeq_fpkm_hgnc_data)#373 patients, they may not unique
    #fileNames<-fileNames1[1:(n-1)]
    fileNames<-fileNames1[seq_len(n-1)]

    if (parallel==TRUE){
        print("cluster making started")
        cl <- makeCluster(mc <- getOption("cl.cores", ncores),type = "FORK")
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

        parLapplyLB(cl,fileNames,processFiles)

        stopCluster(cl)
        print("cluster stopped")
    }
    if (parallel==FALSE){
        lapply(fileNames,processFiles)
    }
}

# ##RUNNING THE FUNCTION
# load("~/CTDPathSim/inst/extdata/Example.Expression.rda")
#
# GetSampleDE(RnaSeq_fpkm_hgnc_data=RnaSeq_fpkm_hgnc_data,parallel=FALSE,ncores=5)
# ##EXAMPLE FILES ARE IN CTDPathSim_Bioconductor/Patient_DE_genes
