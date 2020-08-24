## STAGE 2
## CODE 1
# Created By: Banabithi Bose
# Date Created: 4/23/2020
### EDec for methylation ###
#devtools::install_github("BRL-BCM/EDec")

#Deconvolution of expression profiles

RunDeconvExpr<-function(RnaSeq_fpkm_hgnc_data,Deconv_proportions){
    x1<-data.table(colnames(RnaSeq_fpkm_hgnc_data))
    x2<-data.table(row.names(Deconv_proportions))
    common_samples<-unique(merge(x1,x2,by="V1"))
    stage2_result_tumors = EDec::run_edec_stage_2(
    gene_exp_bulk_samples = RnaSeq_fpkm_hgnc_data[,common_samples$V1],
    cell_type_props = Deconv_proportions[common_samples$V1,])
    Deconv_expression <-stage2_result_tumors$means
    return(Deconv_expression)
}

# # ##RUNNING THE FUNCTION WITH EXAMPLE DATASETS
# load("~/CTDPathSim/inst/extdata/Example.Expression.rda")
#  load("~/CTDPathSim_Bioconductor/Example.Deconv.proportions.450K.rda")
#  Deconv_expression<-RunDeconvExpr(RnaSeq_fpkm_hgnc_data,Deconv_proportions)
#  save(Deconv_expression,file="~/Example.Deconv.expression.rda")
#  ##RESULT HAS BEEN SAVED IN  ~/CTDPathSim_Bioconductor/ WITH NAME Example.Deconv.expression.rda
