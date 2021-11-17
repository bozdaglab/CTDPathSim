\name{RunDeconvMethyl}
\alias{RunDeconvMethyl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to make a deconvulted methylation profile for the genes.
}
\description{
This function is designed to make a deconvoluted methylation profile for the genes.
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Example.Methylation}{
A Dataframe for the example methylation values for the patients and the probes.
    } \item{markers_ovr}{
A vector of the probe markers.
    }
}
\value{
    \item{methylation}{
A R data frame object with the methylation values for each probe and cell type.
    }\item{proportions}{
A R data frame object with the proportion of each cell type for each probe.
    }\item{iterations}{
A double value with the iterations
    }\item{explained.variance}{
A double value with the explained variance
    }\item{res.sum.squares}{
A double value with residual sum of squares
    }\item{aic}{
A double value with the Akaike information criterion
    }\item{rss.per.iteration}{
A vector of doubles for the residual sum of squares per iteration
    }
}
\author{
Banabithi Bose
}
\examples{
stage1_result_c#t_450K.rda
load(system.file("CTDData.RData", package="CTDPathSim2"))
stage1_result_ct<-RunDeconvMethylFun(
                    Methylation.Sample=Methylation.Sample,                                           Reference.Methylation.Probes=Reference.Methylation.Markers)

}
