

ProbeToGeneFun <-
  function(Methylation.Probe.Annotation,
           Deconv_methylation,
           Probe.Annotation) {
    if ((is(Methylation.Probe.Annotation, "SummarizedExpirament"))) {
      Methylation.Probe.Annotation <- assay(Methylation.Probe.Annotation)
    }
    
    if ((is(Deconv_methylation, "SummarizedExpirament"))) {
      Deconv_methylation <- assay(Deconv_methylation)
    }
    
    if (Probe.Annotation != "27K" & Probe.Annotation != "450K") {
      message("Incorrect gene format input. Input either 27K or 450K.")
    }
    
    
    if (Probe.Annotation == "27K") {
      A_27 <- Methylation.Probe.Annotation
      OV_Deconv_methylation <- NULL
      ###############Process Beta Values for 27K###############
      probes <- data.table(row.names(Deconv_methylation))
      colnames(probes) <- "ID"
    } else if (Probe.Annotation == "450K") {
      A_450 <- Methylation.Probe.Annotation
      
      An1 <-
        sqldf("Select * from A_450 where UCSC_RefGene_Group like '%TSS1500%' ")#59310
      An2 <-
        sqldf("Select * from A_450 where UCSC_RefGene_Group like '%TSS200%' ")#46466
      An3 <-
        sqldf("Select * from A_450 where UCSC_RefGene_Group like '%5UTR%' ")#38876
      Probe_Anno_450K <- rbind(An1, An2, An3)#144652
    }
    probes <- data.table(row.names(Deconv_methylation))
    colnames(probes) <- "ID"
    Deconv_methylation <- data.table(cbind(probes, Deconv_methylation))
    
    if (Probe.Annotation == "27K") {
      B <- merge(A_27, Deconv_methylation, by = "ID")
      B$gene <- paste(B$Symbol, ";", B$Synonym)
      B$D <-
        vapply(strsplit(as.character(B$gene), ";", fixed = TRUE), function(Methylation.Probe.Annotation)
          paste(unique(Methylation.Probe.Annotation), collapse = ";"), FUN.VALUE =
            character(1))
    } else if (Probe.Annotation == "450K") {
      B <- merge(Probe_Anno_450K, Deconv_methylation, by = "ID")#125139
      B$D <-
        vapply(strsplit(as.character(B$UCSC_REFGENE_NAME), ";", fixed = TRUE), function(Methylation.Probe.Annotation)
          paste(unique(Methylation.Probe.Annotation), collapse = ";"), FUN.VALUE =
            character(1))
    }
    
    D_Probe.Annotation <- separate_rows(B, D, convert = TRUE)
    D_Probe.Annotation <- data.frame(D_Probe.Annotation)
    n1 <- ncol(D_Probe.Annotation)
    if (Probe.Annotation == "27K") {
      n2 <- n1 - 2
    } else if (Probe.Annotation == "450K") {
      n2 <- n1 - 1
    }
    D_Probe.Annotation_1 <- D_Probe.Annotation[, c(n1, 4:n2)]
    Deconv_meth_gene <-
      aggregate(
        D_Probe.Annotation_1[,-1],
        by = list(D_Probe.Annotation_1$D),
        mean,
        na.rm = TRUE
      )
    
    row.names(Deconv_meth_gene) <- Deconv_meth_gene$Group.1
    Deconv_meth_gene <- Deconv_meth_gene[, -1]
    
    return(Deconv_meth_gene)
    
    
    
  }