\name{GetPath}
\alias{GetPath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to find the reactome pathways for a patient
}
\description{
GetPath finds the reactome pathways for a patient
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{CancerGeneList}{
A R dataframe object with values of frequently mutated cancer driver genes
    }\item{pat_de_genes}{
A R dataframe object with the patients DE genes
    }
}
\value{
Returns a R dataframe object with the reactome pathways for a patient
}
\author{
Banabithi Bose
}
\examples{
data("CTDData")
x<-GetPathFun(CancerGeneList,pat_de_genes)
}
