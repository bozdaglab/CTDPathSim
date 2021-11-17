SampleDeconvProfileFun <- function(Deconv_expression,Deconv_proportions, ProfileType, CTDDirectory="~"){
  if((is(Deconv_expression, "SummarizedExpirament"))){
    Deconv_expression<-assay(Deconv_expression)
  }
  
  if((is(Deconv_proportions, "SummarizedExpirament"))){
    Deconv_proportions<-assay(Deconv_proportions)
  }
  
  Deconv_expression<-abs(Deconv_expression)
  Deconv_proportions<-abs(Deconv_proportions)
  
  if(ProfileType == "Expression"){
    
    if (!file.exists(file.path(CTDDirectory, "Sample_expression"))) {
      dir.create(
        file.path(CTDDirectory, "Sample_expression")
      )
    }
    
  }else{
    if (!file.exists(file.path(CTDDirectory, "Sample_methylation"))) {
      dir.create(
        file.path(CTDDirectory, "Sample_methylation")
      )
    }
  }
  

  processMethExp<- function(sample){
    tryCatch(
      {
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
        
        if(ProfileType == "Expression"){
          save(expression,file=file.path(CTDDirectory, "Sample_expression", paste0(sample, ".Rda")))
        }else{
          methylation<-expression
          save(expression,file=file.path(CTDDirectory, "Sample_methylation", paste0(sample, ".Rda")))
        }
        
        }, error = function(error_condition) {
          if(ProfileType == "Expression"){
            cat(paste0(sample," : ",error_condition),file=file.path("sampleExprFun.txt"),sep="\n",append=TRUE)
          }else{
            cat(paste0(sample," : ",error_condition),file=file.path("sampleMethylFun.txt"),sep="\n",append=TRUE)
          }
        return()
      }, finally={
      }
    )}
  
  samples<-data.table(row.names(Deconv_proportions))
  lapply(samples$V1,processMethExp)
  return()
  
}