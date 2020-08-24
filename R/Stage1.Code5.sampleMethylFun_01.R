## STAGE 1
## CODE 5
#### Making Sample Specific Methylation after Deconvolution ####
# Created By: Banabithi Bose
# Date Created: 4/24/2020


# library(data.table)

sampleMethylFun<-function(Deconv_meth_gene,Deconv_proportions){
    Deconv_meth_gene<-abs(Deconv_meth_gene)
    Deconv_proportions<-abs(Deconv_proportions)
    dir.create("~/Sample_methylation")##MATT TAKE CARE OF THE DIRECTORY

    processMethExp<- function(sample){
        tryCatch(
        {
            #sample=samples$V1[5]
            proportions<-Deconv_proportions[unlist(sample),]
            n1<-ncol(Deconv_meth_gene)
            if(n1==3){
                methylation<-cbind(Deconv_meth_gene[,1]*proportions[1],Deconv_meth_gene[,2]*proportions[2],
                                    Deconv_meth_gene[,3]*proportions[3])
            }
            if(n1==4){
                methylation<-cbind(Deconv_meth_gene[,1]*proportions[1],Deconv_meth_gene[,2]*proportions[2],
                                    Deconv_meth_gene[,3]*proportions[3],Deconv_meth_gene[,4]*proportions[4])
            }
            if(n1==5){
                methylation<-cbind(Deconv_meth_gene[,1]*proportions[1],Deconv_meth_gene[,2]*proportions[2],
                                    Deconv_meth_gene[,3]*proportions[3],Deconv_meth_gene[,4]*proportions[4],
                                    Deconv_meth_gene[,5]*proportions[5])
            }

            methylation<-as.matrix(methylation)
            colnames(methylation)<-colnames(Deconv_meth_gene)
            row.names(methylation)<-row.names(Deconv_meth_gene)
            save(methylation,file=paste0("~/Sample_methylation/",sample,".Rda"))


        }, error = function(error_condition) {
            cat(paste0(sample," : ",error_condition),file="~/processMethExp.txt",sep="\n",append=TRUE)
            return()
        }, finally={
        }
    )}

    samples<-data.table(row.names(Deconv_proportions))
    lapply(samples$V1,processMethExp)
    return()
}


# ##RUNNING THE FUNCTION WITH EXAMPLE DATASETS
# load("~/CTDPathSim/inst/extdata/Example.Deconv.meth.gene.450K.rda")
# load("~/CTDPathSim/inst/extdata/Example.Deconv.proportions.450K.rda")
#
# ##SAMPLE SPECIFIC METHYLATION FILES ARE SAVED IN Sample_methylation FOLDER
# sampleMethylFun(Deconv_meth_gene,Deconv_proportions)
