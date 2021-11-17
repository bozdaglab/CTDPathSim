\name{GetSampleDGenes}
\alias{GetSampleDGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to find the DE genes for each patient
}
\description{
GetSampleDE fins the DE genes for each patient
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{RnaSeq_fpkm_hgnc_data}{
A R dataframe object with the expression values for each probe and cell type
    }\item{parallel}{
A Boolean value for whether or not the user wants to run this function in 
multi-core.
    }\item{ncores}{
An integer value for the number of cores the use if the code is run in parallel
    }
}
\details{
If one is running this function on windows the parallel must be FALSE. 
Although the function does not return anything this function 
makes a R dataframe object with the DE genes for each patient.
}
\value{
There is no return value
}
\author{
Banabithi Bose
}
\examples{
data("CTDData")
GetSampleDGenesFun(Expression.Sample=Expression.Sample,
                   parallel=FALSE, 
                   DType = "DE", 
                   ncores=5, 
                   CTDDirectory = "~/CTD123")
}
