\name{sampleMethylFun}
\alias{sampleMethylFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to find the methylation profile for the patient
}
\description{
SampleMethylFun finds the methylation profile for each patient.
}
\usage{
sampleMethylFun(Deconv_meth_gene,Deconv_proportions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Deconv_meth_gene}{
A R dataframe object with the deconvoluted methylation values for each gene and cell type
    }\item{Deconv_proportions}{
A R dataframe object with the proportions of each cell type for each patient
    }
}
\details{
While there is no return value this function creates a methylation profile for each patient
}
\value{
There is no return value
}
\author{
Banabithi Bose
}
\examples{
##RUNNING THE FUNCTION WITH EXAMPLE DATASETS
load(system.file("extdata", "Example.Deconv.meth.gene.450K.rda", package = "CTDPathSim"))
load(system.file("extdata", "Example.Deconv.proportions.450K.rda", package = "CTDPathSim"))
##SAMPLE SPECIFIC METHYLATION FILES ARE SAVED IN Sample_methylation FOLDER
sampleMethylFun(Deconv_meth_gene,Deconv_proportions)
}
