\name{RunDeconvExpr}
\alias{RunDeconvExpr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to find the deconvoluted methylation profile for each gene
}
\description{
RunDeconvExpr finds the deconvoluted methylation profile for each gene
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
Deconv_expression<-RunDeconvExprFun(Expression.Sample                                                               ,stage1_result_ct$proportions)
}
