## STAGE 1
## CODE 4
### Processing methylation data from probes to genes  ####
# Created By: Banabithi Bose
# Date Created: 4/24/2020

ProbeTogene27Fun<-function(A_27,Deconv_methylation){
    OV_Deconv_methylation<-NULL
    ###############Process Beta Values for 27K###############
    probes<-data.table(row.names(Deconv_methylation))
    colnames(probes)<-"ID"
    Deconv_methylation<-data.table(cbind(probes,Deconv_methylation))
    B<-merge(A_27,OV_Deconv_methylation,by="ID")
    B$gene<-paste(B$Symbol,";",B$Synonym)
    #B$D <- sapply(strsplit(B$gene, ";", fixed = TRUE), function(x)
    #paste(unique(x), collapse = ";"))
    B$D <- vapply(strsplit(as.character(B$gene), ";", fixed = TRUE), function(x)
                    paste(unique(x), collapse = ";"),FUN.VALUE=character(1))

    D_27<-separate_rows(B,D, convert = TRUE)#88132
    D_27<-data.frame(D_27)
    n1<-ncol(D_27)
    n2<-n1-2
    D_27_1<-D_27[,c(n1,4:n2)]  ### THIS PART NEEDS TO BE GENERALIZED BASED ON CELL TYPES
    Deconv_meth_gene <- aggregate(D_27_1[, -1], by = list(D_27_1$D),
                                    mean, na.rm = TRUE) #43248

    row.names(Deconv_meth_gene)<-Deconv_meth_gene$Group.1
    Deconv_meth_gene<-Deconv_meth_gene[,-1]

    return(Deconv_meth_gene)
}

# ##RUNNING THE FUNCTION
# ## EXAMPLE A_27 FILE IS IN EXAMPLE DATASETS WITH NAME Methylation.Probe.Annotation27.rda
# Deconv_meth_gene<-ProbeTogene27Fun(A_27,Deconv_methylation)
