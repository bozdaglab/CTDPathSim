\name{sampleExprFun}
\alias{sampleExprFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to create an expression profile for each patient
}
\description{
sampleExprFun creates an expression profile for each patient
}
\usage{
sampleExprFun(Deconv_expression,Deconv_proportions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Deconv_expression}{
A R dataframe object with an expression value for each gene and cell type
    }\item{Deconv_proportions}{
A R dataframe object with proportions of each cell type for each patient
    }
}
\details{
Although there is no return value this function creates an expression profile for each patient
}
\value{
There is no return value
}
\author{
Banabithi Bose
}
\examples{
##RUNNING THE FUNCTION WITH EXAMPLE DATASETS
load(system.file("extdata", "Example.Deconv.expression.rda", package = "CTDPathSim"))
load(system.file("extdata", "Example.Deconv.proportions.450K.rda", package = "CTDPathSim"))
##SAMPLE SPECIFIC EXPRESSION FILES ARE SAVED IN Sample_expression FOLDER
sampleExprFun(Deconv_expression,Deconv_proportions)
}
