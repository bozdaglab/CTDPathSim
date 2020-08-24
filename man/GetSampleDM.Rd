\name{GetSampleDM}
\alias{GetSampleDM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to find the DM genes for each patient
}
\description{
GetSampleDM fins the DM genes for each patient
}
\usage{
GetSampleDM(All_methylation,parallel=c('TRUE','FALSE'),ncores=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{All_methylation}{
A R dataframe object with the methylation value for each gene and patient
    }\item{parallel}{
A Boolean value for whether or not the user wants to run this function in multi-core. If the user is on windows parallel must be FALSE
    }\item{ncores}{
An integer value for the number of cores the use if the code is run in parallel
    }
}
\details{
If one is running this function on windows the parallel must be FALSE. Although the function does not return anything this function makes a R dataframe object with the DM genes for each patient.
}
\value{
There is no return value
}
\author{
Banabithi Bose
}
\examples{
load(system.file("extdata", "Example.gene.centric.methylation.50.rda", package = "CTDPathSim"))
GetSampleDM(All_methylation=All_methylation,parallel=FALSE,ncores=5)
}
