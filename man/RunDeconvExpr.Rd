\name{RunDeconvExpr}
\alias{RunDeconvExpr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to find the deconvoluted methylation profile for each gene
}
\description{
RunDeconvExpr finds the deconvoluted methylation profile for each gene
}
\usage{
RunDeconvExpr(RnaSeq_fpkm_hgnc_data,Deconv_proportions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{RnaSeq_fpkm_hgnc_data}{
A R dataframe object with the expression values for each probe and cell type
    }\item{Deconv_proportions}{
A R dataframe object with the proportion of each cell type for each patient
    }
}
\value{
A R dataframe object with the expression values for each gene and cell type
}
\author{
Banabithi Bose
}
\examples{
load(system.file("extdata", "Example.Expression.50.rda", package = "CTDPathSim"))
load(system.file("extdata", "Example.Deconv.proportions.450K.rda", package = "CTDPathSim"))
Deconv_expression<-RunDeconvExpr(RnaSeq_fpkm_hgnc_data,Deconv_proportions)
save(Deconv_expression,file="~/Example.Deconv.expression.rda")
}
