\name{PlotDeconvMethyl}
\alias{PlotDeconvMethyl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to create a correlation plot for cell types with clusters predicted after deconvotion.
}
\description{
PlotDeconvMethyl creates a correlation plot for cell types with clusters predicted after deconvolution.
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
data("CTDData")
print(PlotDeconvMethylFun(Methylation.Sample=Methylation.Sample,
Reference.Methylation.Probes=Reference.Methylation.Probes,
Reference.Methylation.CellTypes=Reference.Methylation.CellTypes,
stage1_result_ct=stage1_result_ct,
Reference.CellTypes.Names=Reference.CellTypes.Names))
dev.off()
}
