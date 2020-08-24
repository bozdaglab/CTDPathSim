## STAGE 2
## CODE 2
#### Making Sample Specific Expression after Deconvolution ####
# Created By: Banabithi Bose
# Date Created: 4/24/2020

######## For 450k samples ########
library(data.table)
load("~/PATIENT_CELL_PIPELINE3/OVARIAN/EDEC/OV_Deconv_expression450.Rda")
OV_Deconv_expression<-abs(OV_Deconv_expression)
load("~/PATIENT_CELL_PIPELINE3/OVARIAN/EDEC/OV_Deconv_meth450_gene.Rda")
OV_Deconv_meth_gene<-abs(OV_Deconv_meth_gene)
load("~/PATIENT_CELL_PIPELINE3/OVARIAN/EDEC/OV_Deconv_proportions_450.Rda")
OV_Deconv_proportions<-abs(OV_Deconv_proportions)
processMethExp<- function(sample){
  tryCatch(
    {
      #sample=samples$V1[1]
      proportions<-OV_Deconv_proportions[unlist(sample),]
      methylation<-cbind(OV_Deconv_meth_gene[,1]*proportions[1],OV_Deconv_meth_gene[,2]*proportions[2],
                         OV_Deconv_meth_gene[,3]*proportions[3])
      methylation<-as.matrix(methylation)
      colnames(methylation)<-colnames(OV_Deconv_meth_gene)
      row.names(methylation)<-row.names(OV_Deconv_meth_gene)
      expression<-cbind(OV_Deconv_expression[,1]*proportions[1],OV_Deconv_expression[,2]*proportions[2],
                        OV_Deconv_expression[,3]*proportions[3])
      expression<-as.matrix(expression)
      colnames(expression)<-colnames(OV_Deconv_expression)
      row.names(expression)<-row.names(OV_Deconv_expression)
      #save(methylation,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Sample_methylation/",sample,".Rda"))
      save(methylation,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Sample_methylation/",str_sub(sample,1,16),".Rda"))
      
      #save(expression,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Sample_expression/",sample,".Rda"))
      save(expression,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Sample_expression/",str_sub(sample,1,16),".Rda"))
      
    }, error = function(error_condition) {
      cat(paste0(sample," : ",error_condition),file="~/PATIENT_CELL_PIPELINE3/ERROR_FILES/processMethExp.txt",sep="\n",append=TRUE)
      return()
    }, finally={
      
    }
  )
}

samples<-data.table(row.names(OV_Deconv_proportions))
lapply(samples$V1,processMethExp)

#fileNamesMeth<-list.files("~/PATIENT_CELL_PIPELINE3/OVARIAN/Sample_methylation/", pattern="*.Rda") #10
#fileNamesExpr<-list.files("~/PATIENT_CELL_PIPELINE3/OVARIAN/Sample_expression/", pattern="*.Rda") #10

######## For 27k samples ########

load("~/PATIENT_CELL_PIPELINE3/OVARIAN/EDEC/OV_Deconv_expression27.Rda")
OV_Deconv_expression<-abs(OV_Deconv_expression)
load("~/PATIENT_CELL_PIPELINE3/OVARIAN/EDEC/OV_Deconv_meth27_gene.Rda")
OV_Deconv_meth_gene<-abs(OV_Deconv_meth_gene)
load("~/PATIENT_CELL_PIPELINE3/OVARIAN/EDEC/OV_Deconv_proportions_27.Rda")
OV_Deconv_proportions<-abs(OV_Deconv_proportions)


processMethExp<- function(sample){
  tryCatch(
    {
      #sample=samples$V1[1]
      proportions<-OV_Deconv_proportions[unlist(sample),]
      methylation<-cbind(OV_Deconv_meth_gene[,1]*proportions[1],OV_Deconv_meth_gene[,2]*proportions[2],
                         OV_Deconv_meth_gene[,3]*proportions[3],OV_Deconv_meth_gene[,4]*proportions[4])
      methylation<-as.matrix(methylation)
      colnames(methylation)<-colnames(OV_Deconv_meth_gene)
      row.names(methylation)<-row.names(OV_Deconv_meth_gene)
      expression<-cbind(OV_Deconv_expression[,1]*proportions[1],OV_Deconv_expression[,2]*proportions[2],
                        OV_Deconv_expression[,3]*proportions[3],OV_Deconv_expression[,3]*proportions[3])
      expression<-as.matrix(expression)
      colnames(expression)<-colnames(OV_Deconv_expression)
      row.names(expression)<-row.names(OV_Deconv_expression)
      #save(methylation,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Sample_methylation/",sample,".Rda"))
      save(methylation,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Sample_methylation/",str_sub(sample,1,16),".Rda"))
      
      #save(expression,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Sample_expression/",sample,".Rda"))
      save(expression,file=paste0("~/PATIENT_CELL_PIPELINE3/OVARIAN/Sample_expression/",str_sub(sample,1,16),".Rda"))
      
    }, error = function(error_condition) {
      cat(paste0(sample," : ",error_condition),file="~/PATIENT_CELL_PIPELINE3/ERROR_FILES/processMethExp.txt",sep="\n",append=TRUE)
      return()
    }, finally={
      
    }
  )
}

samples<-data.table(row.names(OV_Deconv_proportions))#613
lapply(samples$V1,processMethExp)







