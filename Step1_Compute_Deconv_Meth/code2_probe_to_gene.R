## STAGE 1 
## CODE 2
### Processing methylation data from probes to genes  ####
# Created By: Banabithi Bose
# Date Created: 4/24/2020

load("~/PATIENT_CELL_PIPELINE3/OVARIAN/EDEC/OV_Deconv_expression450.Rda")

###############Process Beta Values for 450K###############
Annonation_450<- read.table("~/PATIENT_CELL_PIPELINE1/GPL18809-31921.txt", sep="\t",header = T)
A_450 <- Annonation_450[,c(1,25,27)]#242079
An1<- sqldf("Select * from A_450 where UCSC_RefGene_Group like '%TSS1500%' ")#59310
An2<- sqldf("Select * from A_450 where UCSC_RefGene_Group like '%TSS200%' ")#46466
An3<- sqldf("Select * from A_450 where UCSC_RefGene_Group like '%5UTR%' ")#38876
Probe_Anno_450K <- rbind(An1,An2,An3)#144652

probes<-data.table(row.names(OV_Deconv_methylation))
colnames(probes)<-"ID"
OV_Deconv_methylation<-data.table(cbind(probes,OV_Deconv_methylation))
B<-merge(Probe_Anno_450K,OV_Deconv_methylation,by="ID")#125139
B$D <- sapply(strsplit(as.character(B$UCSC_REFGENE_NAME), ";", fixed = TRUE), function(x) 
  paste(unique(x), collapse = ";"))


library(tidyr)########## few lines of code missing
OV_450<-separate_rows(B,D, convert = TRUE)
OV_450_1<-OV_450[,c(7,4:6)] #152210
OV_Deconv_meth_gene <- aggregate(OV_450_1[, -1], by = list(OV_450_1$D), 
                                   mean, na.rm = TRUE) #19395

row.names(OV_Deconv_meth_gene)<-OV_Deconv_meth_gene$Group.1
OV_Deconv_meth_gene<-OV_Deconv_meth_gene[,-1]
save(OV_Deconv_meth_gene,file="~/PATIENT_CELL_PIPELINE3/OVARIAN/EDEC/OV_Deconv_meth450_gene.Rda")

### Processing 27 k methylation data from probes to genes 
load("~/PATIENT_CELL_PIPELINE3/OVARIAN/EDEC/OV_Deconv_methylation27.Rda")

############### Process Beta Values for 27K ###############
Annonation_27<- fread("~/PATIENT_CELL_PIPELINE1/GPL8490-65.txt", sep="\t",header = T)
A_27 <- Annonation_27[,c(1,23,24)]
probes<-data.table(row.names(OV_Deconv_methylation))
colnames(probes)<-"ID"
OV_Deconv_methylation<-data.table(cbind(probes,OV_Deconv_methylation))
B<-merge(A_27,OV_Deconv_methylation,by="ID")#21675
B$gene<-paste(B$Symbol,";",B$Synonym)
B$D <- sapply(strsplit(B$gene, ";", fixed = TRUE), function(x) 
  paste(unique(x), collapse = ";"))
OV_27<-separate_rows(B,D, convert = TRUE)#88132
OV_27_1<-OV_27[,c(9,4:7)] #88132 ### THIS PART NEEDS TO BE GENERALIZED BASED ON CELL TYPES
OV_Deconv_meth_gene <- aggregate(OV_27_1[, -1], by = list(OV_27_1$D), 
                                   mean, na.rm = TRUE) #43248

row.names(OV_Deconv_meth_gene)<-OV_Deconv_meth_gene$Group.1
OV_Deconv_meth_gene<-OV_Deconv_meth_gene[,-1]
save(OV_Deconv_meth_gene,file="~/PATIENT_CELL_PIPELINE3/OVARIAN/EDEC/OV_Deconv_meth27_gene.Rda")
