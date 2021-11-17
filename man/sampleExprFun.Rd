\name{SampleDeconvProfileFun}
\alias{SampleDeconvProfileFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to create an expression profile for each patient
}
\description{
SampleDeconvProfileFun creates an expression profile for each patient
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
Although there is no return value this function creates an expression 
profile for each patient
}
\value{
There is no return value
}
\author{
Banabithi Bose
}
\examples{
##SAMPLE SPECIFIC EXPRESSION FILES ARE SAVED IN Sample_expression FOLDER
load(system.file("CTDData.RData", package="CTDPathSim2"))
Deconv_proportions<-stage1_result_ct$proportions
SampleDeconvProfileFun(Deconv_expression,
                       Deconv_proportions, 
                       ProfileType="Expression", 
                       CTDDirectory="~/CTD123")
}
