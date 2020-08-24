\name{RunDeconvMethyl}
\alias{RunDeconvMethyl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to make a deconvulted methylation profile for the genes.
}
\description{
This function is designed to make a deconvoluted methylation profile for the genes.
}
\usage{
RunDeconvMethyl(Example.Methylation, markers_ovr)
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
print("Function to make a deconvulted methylation profile for the genes.")
##RUNNING THE FUNCTION WITH EXAMPLE DATASETS
#load("~/CTDPathSim_Bioconductor/Example.Methylation.450K.50.rda")
#load("~/CTDPathSim_Bioconductor/markers_ovr_450.rda")
#stage1_result_ct<-RunDeconvMethyl(Example.Methylation=Example.Methylation.450K
#markers_ovr=markers_ovr)
##RESULT HAS BEEN SAVED IN  ~/CTDPathSim_Bioconductor/ WITH NAME stage1_result_c#t_450K.rda
}
