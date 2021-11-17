
GetSampleDGenesFun<-function(Expression.Sample, ncore = 1, DType, CTDDirectory="~"){
    
  
  if((is(Expression.Sample, "SummarizedExpirament"))){
    Expression.Sample<-assay(Expression.Sample)
  }
  
  if(DType == "DE"){
    
    if (!file.exists(file.path(CTDDirectory, "Patient_DE_genes"))) {
      dir.create(
        file.path(CTDDirectory, "Patient_DE_genes")
      )
    if (!file.exists(file.path(CTDDirectory, "Patient_DE_genes", "UP_Gene"))) {
        dir.create(
        file.path(CTDDirectory, "Patient_DE_genes", "UP_Gene")
      )
    
      if (!file.exists(file.path(CTDDirectory, "Patient_DE_genes", "DOWN_Gene"))) {
        dir.create(
          file.path(CTDDirectory, "Patient_DE_genes", "DOWN_Gene")
        )
        
        if (!file.exists(file.path(CTDDirectory, "Patient_DE_genes", "DE_Gene"))) {
          dir.create(
            file.path(CTDDirectory, "Patient_DE_genes", "DE_Gene")
          )
        }
      }
    }
    }
  }else if(DType == "DM"){
    if (!file.exists(file.path(CTDDirectory, "Patient_DM_genes"))) {
      dir.create(
        file.path(CTDDirectory, "Patient_DM_genes")
      )
      if (!file.exists(file.path(CTDDirectory, "Patient_DM_genes", "UP_Gene"))) {
        dir.create(
          file.path(CTDDirectory, "Patient_DM_genes", "UP_Gene")
        )
        
        if (!file.exists(file.path(CTDDirectory, "Patient_DM_genes", "DOWN_Gene"))) {
          dir.create(
            file.path(CTDDirectory, "Patient_DM_genes", "DOWN_Gene")
          )
          if (!file.exists(file.path(CTDDirectory, "Patient_DM_genes", "DM_Gene"))) {
            dir.create(
              file.path(CTDDirectory, "Patient_DM_genes", "DM_Gene")
            ) 
          }
        }
      }
  }
  
  
  
  processFiles <- function(fileNames){
    
    SampleName <- fileNames
    Expression.Sample<-as.data.frame(Expression.Sample)

    gene_expr<-Expression.Sample[, SampleName,drop=FALSE]
    gene_med<-Expression.Sample[,c(1,2),drop=FALSE]
    gene_expr_med<-cbind(gene_med,gene_expr)
    colnames(gene_expr_med)<-c("gene","median_expr","expr")
    gene_expr_med[gene_expr_med$median_expr==0,2]<-0.00000000001 ## to avoid divison by 0 for fold change
    gene_expr_med[gene_expr_med$expr==0,3]<-0.00000000001 ## to get non zero fold change, removing fold change = 0 for down gene
    gene_expr_med$fold_change<-gene_expr_med$expr/gene_expr_med$median_expr
    upgene<-gene_expr_med[gene_expr_med$fold_change>=4,]
    downgene<-gene_expr_med[gene_expr_med$fold_change<=0.25,]
    pat_de_genes <- data.table(rbind(upgene,downgene))[,1]
    
    if(DType == "DM"){
      
      save(upgene,file=file.path(CTDDirectory, "Patient_DM_genes", "UP_Gene", paste0(SampleName, "_upgene.Rda")))
      save(downgene,file=file.path(CTDDirectory, "Patient_DM_genes", "DOWN_Gene", paste0(SampleName, "_downgene.Rda")))
      save(pat_de_genes,file=file.path(CTDDirectory, "Patient_DM_genes", "DM_Gene", paste0(SampleName, "_pat_de_genes.Rda")))
      
    }else{
      save(upgene,file=file.path(CTDDirectory, "Patient_DE_genes", "UP_Gene", paste0(SampleName, "_upgene.Rda")))
      save(downgene,file=file.path(CTDDirectory, "Patient_DE_genes", "DOWN_Gene", paste0(SampleName, "_downgene.Rda")))
      save(pat_de_genes,file=file.path(CTDDirectory, "Patient_DE_genes", "DE_Gene", paste0(SampleName, "_pat_de_genes.Rda")))
    }
    
    
  }
  
  
  if(DType == "DM"){
    Expression.Sample<-as.data.frame(Expression.Sample)
    Expression.Sample[Expression.Sample==0]<-0.00000000001 ## to avoid the log transformation of 0
    #Writing a function for transforming beta values to M values
    Mtranform <- function(x) {log(x[1]/(1-x[1]))}
    n<-ncol(Expression.Sample)
    Expression.Sample[, seq_len(n) ] <- as.data.frame(lapply(Expression.Sample[, seq_len(n) ], FUN = function(x) {vapply(x, FUN = Mtranform, FUN.VALUE = 1.0)}))
    x <- data.table(rowMedians(as.matrix(Expression.Sample),na.rm=TRUE))#37109 genes total in hgnc
    
  }else{
    x <- data.table(rowMedians(Expression.Sample,na.rm=TRUE))#37109 genes total in hgnc
  }
  
  
  colnames(x)<-"Median_expr"
  row.names(x)<-c(rownames(Expression.Sample))

  Rnames<-rownames(Expression.Sample)
  Expression.Sample<-cbind(x,Expression.Sample)
  
  Expression.Sample<-as.data.frame(Expression.Sample)
  Expression.Sample<-cbind(Rnames, Expression.Sample)
  
  n<-ncol(Expression.Sample)
  n<-as.integer(n)
  fileNames1<-colnames(Expression.Sample)
  
  fileNames<-fileNames1[3:n]
  
  if (ncore > 1){
    print("cluster making started")
    cl <- makeCluster(mc <- getOption("cl.cores", ncores))
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
    parallel::clusterExport(cl=cl,c("Expression.Sample"),envir=environment())
    print("cluster making done")
    
    parLapplyLB(cl,fileNames,processFiles)
    
    stopCluster(cl)
    print("cluster stopped")
  }else{
    lapply(fileNames,processFiles)
  }
  return()
}}