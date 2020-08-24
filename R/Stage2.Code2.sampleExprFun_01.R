## STAGE 1
## CODE 5
#### Making Sample Specific Expression after Deconvolution ####
# Created By: Banabithi Bose
# Date Created: 4/24/2020

# library(data.table)

sampleExprFun<-function(Deconv_expression,Deconv_proportions){
    Deconv_expression<-abs(Deconv_expression)
    Deconv_proportions<-abs(Deconv_proportions)
    dir.create("~/Sample_expression")##MATT TAKE CARE OF THE DIRECTORY

    processMethExp<- function(sample){
        tryCatch(
        {
            #sample=samples$V1[5]
            proportions<-Deconv_proportions[unlist(sample),]
            n1<-ncol(Deconv_expression)
            if(n1==3){
                expression<-cbind(Deconv_expression[,1]*proportions[1],Deconv_expression[,2]*proportions[2],
                                    Deconv_expression[,3]*proportions[3])
            }
            if(n1==4){
                expression<-cbind(Deconv_expression[,1]*proportions[1],Deconv_expression[,2]*proportions[2],
                                    Deconv_expression[,3]*proportions[3],Deconv_expression[,4]*proportions[4])
            }
            if(n1==5){
                expression<-cbind(Deconv_expression[,1]*proportions[1],Deconv_expression[,2]*proportions[2],
                                    Deconv_expression[,3]*proportions[3],Deconv_expression[,4]*proportions[4],
                                    Deconv_expression[,5]*proportions[5])
            }

            expression<-as.matrix(expression)
            colnames(expression)<-colnames(Deconv_expression)
            row.names(expression)<-row.names(Deconv_expression)
            save(expression,file=paste0("~/Sample_expression/",sample,".Rda"))
        }, error = function(error_condition) {
            cat(paste0(sample," : ",error_condition),file="~/sampleExprFun.txt",sep="\n",append=TRUE)
            return()
        }, finally={
        }
    )}

    samples<-data.table(row.names(Deconv_proportions))
    lapply(samples$V1,processMethExp)
    return()
}

# ##RUNNING THE FUNCTION WITH EXAMPLE DATASETS
# load("~/CTDPathSim/inst/extdata/Example.Deconv.expression.rda")
# load("~/CTDPathSim/inst/extdata/Example.Deconv.proportions.450K.rda")
# ##SAMPLE SPECIFIC EXPRESSION FILES ARE SAVED IN Sample_expression FOLDER
# sampleExprFun(Deconv_expression,Deconv_proportions)
