## STAGE 2
## CODE 1
# Created By: Banabithi Bose
# Date Created: 4/23/2020

RunDeconvExprFun<-function(Expression.Sample,Deconv_proportions){
    
    if((is(Expression.Sample, "SummarizedExpirament"))){
        Expression.Sample<-assay(Expression.Sample)
    }
  
    if((is(Deconv_proportions, "SummarizedExpirament"))){
        Deconv_proportions<-assay(Deconv_proportions)
    }
  
    x1<-data.table(colnames(Expression.Sample))
    x2<-data.table(row.names(Deconv_proportions))
    common_samples<-unique(merge(x1,x2,by="V1"))
    stage2_result_tumors = run_edec_stage_2(
    gene_exp_bulk_samples = Expression.Sample[,common_samples$V1],
    cell_type_props = Deconv_proportions[common_samples$V1,])
    Deconv_expression <-stage2_result_tumors$means
    return(Deconv_expression)
}
