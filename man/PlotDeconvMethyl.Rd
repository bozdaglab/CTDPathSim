\name{PlotDeconvMethyl}
\alias{PlotDeconvMethyl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to create a correlation plot for cell types with clusters predicted after deconvotion.
}
\description{
PlotDeconvMethyl creates a correlation plot for cell types with clusters predicted after deconvolution.
}
\usage{
PlotDeconvMethyl(Example.Methylation,markers_ovr,complete_ref_meth,stage1_result_ct,cn)
}
\arguments{
    \item{Example.Methylation}{
A R dataframe object with the expression values for each cell type
    }\item{markers_ovr}{
A R dataframe object that includes all the probes used in the study
    }\item{complete_ref_meth}{
A R dataframe object with the reference methylation values.
    }\item{stage1_result_ct}{
A list that contains the otuput of runDeconvMethyl
    }\item{cn}{
A R dataframe object with the genes and their class
    }
}
\details{
This function has no return value but creates a pdf of the correlation plot.
}
\value{
This function has no return value.
}
\author{
Banabithi Bose
}
\examples{
print("This function creates a correlation plot for cell types with clusters after deconvolution.")
#load(system.file("extdata", "Example.Methylation.450K.rda", package = "CTDPathSim"))
#load(system.file("extdata", "complete_ref_meth450.rda", package = "CTDPathSim"))
#load(system.file("extdata", "markers_ovr_450.rda", package = "CTDPathSim"))
#load(system.file("extdata", "stage1_result_ct_450K.rda", package = "CTDPathSim"))
# Create a vector of colors representing the class of each reference
#cn<-xlsx::read.xlsx2(file = "~/CTDPathSim/inst/extdata/ref_names.xlsx",header=TRUE,sheetIndex = 1)
#pdf("~/CTDPathSim/inst/extdata/_cor_matrix_stage1_450K.pdf")
#PlotDeconvMethyl(Example.Methylation=Example.Methylation.450K,markers_ovr=markers_ovr,
#                    complete_ref_meth=complete_ref_meth450,stage1_result_ct=stage1_result_ct,cn=cn)
#dev.off()
}
