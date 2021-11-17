# ## STAGE 3
# ## CODE 5
# #### Getting cell line specific DA genes  ####
# # Created By: Banabithi Bose
# # Date Created: 4/27/2020

GetCellLineDAFun <- function(CNV.Sample,
                             ncore = 1,
                             CTDDirectory = "~") {
  if ((is(CNV.Sample, "SummarizedExpirament"))) {
    CNV.Sample <- assay(CNV.Sample)
  }
  if (!file.exists(file.path(CTDDirectory, "CCLE_DA_genes"))) {
    dir.create(file.path(CTDDirectory, "CCLE_DA_genes"))
  }
  if (!file.exists(file.path(CTDDirectory, "CCLE_DA_genes", "UP_Gene"))) {
    dir.create(file.path(CTDDirectory, "CCLE_DA_genes", "UP_Gene"))
  }
  if (!file.exists(file.path(CTDDirectory, "CCLE_DA_genes", "DOWN_Gene"))) {
    dir.create(file.path(CTDDirectory, "CCLE_DA_genes", "DOWN_Gene"))
  }
  if (!file.exists(file.path(CTDDirectory, "CCLE_DA_genes", "DE_Gene"))) {
    dir.create(file.path(CTDDirectory, "CCLE_DA_genes", "DE_Gene"))
  }
  
  processFiles <- function(i) {
    #i=1
    d <-
      data.frame(CNV.Sample[, i] - rowMeans(CNV.Sample[,-i]))
    m <- data.frame(rowMeans(CNV.Sample[,]))
    df <- cbind(d, m)
    df <- na.omit(df)
    colnames(df) <- c("d", "m")
    df$z <- "ok"
    df[df$d >= 1 & df$m >= .5, 3] <- "upda"
    df[df$d <= (-1) & df$m >= .5, 3] <- "downda"
    SampleName <- colnames(CNV.Sample[, i, drop = FALSE])
    ampgene <- data.table(row.names(df[df$z == "upda",]))
    delgene <- data.table(row.names(df[df$z == "downda",]))
    cell_da_genes <- data.table(rbind(ampgene, delgene))
    
    save(ampgene,
         file = file.path(
           CTDDirectory,
           "CCLE_DA_genes",
           "UP_Gene",
           paste0(SampleName, "_upgene.Rda")
         ))
    save(delgene,
         file = file.path(
           CTDDirectory,
           "CCLE_DA_genes",
           "DOWN_Gene",
           paste0(SampleName, "_downgene.Rda")
         ))
    save(cell_da_genes,
         file = file.path(
           CTDDirectory,
           "CCLE_DA_genes",
           "DE_Gene",
           paste0(SampleName, "_pat_de_genes.Rda")
         ))
  }
  n <- ncol(CNV.Sample)
  n <- as.integer(n)
  i <- as.integer(1:n)
  if (ncore > 1) {
    print("cluster making started")
    cl <- makeCluster(mc <- getOption("cl.cores", ncores))
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
    parallel::clusterExport(cl = cl, NULL, envir = environment())
    parLapplyLB(cl, i, processFiles)
    stopCluster(cl)
    
  } else {
    lapply(i, processFiles)
  }
  return()
}
