## Checked intersect of de genes between two patients, it was around 2000 genes
## After this we need to prepare deconvoluted expression data with these high/low genes

##### Use BRCA Methylation data (before deconvolution) to get median and fold change #####
#Compute median expression of a gene across all samples

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

GetSampleDM<-function(All_methylation,parallel=c('TRUE','FALSE'),ncores=2){

    processFiles <- function(fileNames){

        #SampleName <- fileNames[1]
        ##MATT TAKE CARE OF THE FOLDER
        dir.create("~/Patient_DM_genes/")
        dir.create("~/Patient_DM_genes/UP_Gene/")
        dir.create("~/Patient_DM_genes/DOWN_Gene/")
        dir.create("~/Patient_DM_genes/DE_Gene/")

        #######
        #SampleName <- fileNames[1]
        SampleName <- fileNames
        All_methylation<-as.data.frame(All_methylation)
        gene_meth<-All_methylation[,SampleName,drop=FALSE]
        gene_med<-All_methylation[,c(1,n),drop=FALSE]
        gene_meth_med<-cbind(gene_med,gene_meth)
        colnames(gene_meth_med)<-c("gene","median_meth","meth")
        gene_meth_med[gene_meth_med$median_meth==0,2]<-0.00000000001 ## to avoid divison by 0 for fold change
        gene_meth_med[gene_meth_med$meth==0,3]<-0.00000000001 ## to get non zero fold change, removing fold change = 0 for down gene
        gene_meth_med$fold_change<-abs(gene_meth_med$meth/gene_meth_med$median_meth)
        ### Since fold change >2 or <.5 giving 10,000 genes among 40,000, I am taking 4 and .25 threshold,
        # new threshold is giving around 5000 genes
        upgene<-gene_meth_med[gene_meth_med$fold_change>=4,]## 271 genes
        downgene<-gene_meth_med[gene_meth_med$fold_change<=0.25,]## 449 genes
        pat_dm_genes <- data.table(rbind(upgene,downgene))[,1]


        save(upgene,file=paste0("~/Patient_DM_genes/UP_Gene/",SampleName,"_upgene.Rda"))
        save(downgene,file=paste0("~/Patient_DM_genes/DOWN_Gene/",SampleName,"_downgene.Rda"))
        save(pat_dm_genes,file=paste0("~/Patient_DM_genes/DE_Gene/",SampleName,"_pat_dm_genes.Rda"))

    }


    All_methylation<-as.data.frame(All_methylation)
    All_methylation[All_methylation==0]<-0.00000000001 ## to avoid the log transformation of 0
    #Writing a function for transforming beta values to M values
    Mtranform <- function(x) {log(x[1]/(1-x[1]))}
    n<-ncol(All_methylation)
    ##All_methylation[, 1:n ] <- as.data.frame(lapply(All_methylation[, 1:n], FUN = function(x) {sapply(x, FUN = Mtranform)}))
    ## All_methylation[, seq_len(n) ] <- as.data.frame(lapply(All_methylation[, seq_len(n) ], FUN = function(x) {sapply(x, FUN = Mtranform)}))
    All_methylation[, seq_len(n) ] <- as.data.frame(lapply(All_methylation[, seq_len(n) ], FUN = function(x) {vapply(x, FUN = Mtranform, FUN.VALUE = 1.0)}))
    ##All_methylation[, seq_len(n) ] <- as.data.frame(lapply(All_methylation[, seq_len(n) ], FUN = function(x) { )}))


    ##library(matrixStats)
    #All_methylation$Median_meth <- rowMedians(All_methylation,na.rm=TRUE)#14808 genes total in hgnc
    x <- data.table(rowMedians(as.matrix(All_methylation),na.rm=TRUE))#37109 genes total in hgnc
    colnames(x)<-"Median_expr"
    All_methylation<-cbind(All_methylation,x)
    n<-ncol(All_methylation)
    n<-as.integer(n)
    fileNames1<-colnames(All_methylation)#373 patients, they may not unique
    ##fileNames<-fileNames1[1:(n-1)]
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
        parallel::clusterExport(cl=cl,c("All_methylation"),envir=environment())#envir needed to be correct, export from global aswell as local
        print("cluster making done")

        parLapplyLB(cl,fileNames,processFiles)

        stopCluster(cl)
        print("cluster stopped")
    }
    if (parallel==FALSE){
        lapply(fileNames,processFiles)
    }
}
#RUNNING THE FUNCTION
# load("~/CTDPathSim/inst/extdata/Example.gene.centric.methylation.50.rda")
#
# GetSampleDM(All_methylation=All_methylation,parallel=FALSE,ncores=5)
##EXAMPLE FILES ARE IN CTDPathSim_Bioconductor/Patient_DM_genes
