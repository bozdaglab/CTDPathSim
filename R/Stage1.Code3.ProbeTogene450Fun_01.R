## STAGE 1
## CODE 3
### Processing methylation data from probes to genes  ####
# Created By: Banabithi Bose
# Date Created: 4/24/2020

ProbeTogene450Fun<-function(A_450,Deconv_methylation){
    ###############Process Beta Values for 450K###############
    An1<- sqldf("Select * from A_450 where UCSC_RefGene_Group like '%TSS1500%' ")#59310
    An2<- sqldf("Select * from A_450 where UCSC_RefGene_Group like '%TSS200%' ")#46466
    An3<- sqldf("Select * from A_450 where UCSC_RefGene_Group like '%5UTR%' ")#38876
    Probe_Anno_450K <- rbind(An1,An2,An3)#144652

    probes<-data.table(row.names(Deconv_methylation))
    colnames(probes)<-"ID"
    Deconv_methylation<-data.table(cbind(probes,Deconv_methylation))
    B<-merge(Probe_Anno_450K,Deconv_methylation,by="ID")#125139
    #B$D <- sapply(strsplit(as.character(B$UCSC_REFGENE_NAME), ";", fixed = TRUE), function(x)
    #paste(unique(x), collapse = ";"))
    ##B$D <- vapply(strsplit(as.character(B$UCSC_REFGENE_NAME),";"), function(x)
    ##paste(unique(x), FUN.VALUE = ";"))
    B$D <- vapply(strsplit(as.character(B$UCSC_REFGENE_NAME), ";", fixed = TRUE), function(x)
                    paste(unique(x), collapse = ";"),FUN.VALUE=character(1))

    ##library(tidyr)########## few lines of code missing
    D_450<-separate_rows(B,D, convert = TRUE)
    D_450<-data.frame(D_450)
    n1<-ncol(D_450)
    n2<-n1-1
    D_450_1<-D_450[,c(n1,4:n2)]
    Deconv_meth_gene <- aggregate(D_450_1[, -1], by = list(D_450_1$D), mean, na.rm = TRUE)

    row.names(Deconv_meth_gene)<-Deconv_meth_gene$Group.1
    Deconv_meth_gene<-Deconv_meth_gene[,-1]

    return(Deconv_meth_gene)
}

##RUNNING THE FUNCTION AND SAVING THE RESULT
# load("~/CTDPathSim/inst/extdata/Methylation.Probe.Annotation450.rda")
# load("~/CTDPathSim/inst/extdata/Example.Deconv.methylation.450K.rda")
#
# Deconv_meth_gene<-ProbeTogene450Fun(A_450=A_450,Deconv_methylation=Deconv_methylation)
# ## EXAMPLE Deconv_meth_gene  IN EXAMPLE DATASETS UNDER NAME Example.Deconv.meth.gene.450K.rda
