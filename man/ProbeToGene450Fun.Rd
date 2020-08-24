\name{ProbeTogene450Fun}
\alias{ProbeTogene450Fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to convert the probes to genes.
}
\description{
This function converts dataframes with probes to dataframes with genes.
}
\usage{
ProbeTogene450Fun(A_450,Deconv_methylation)
}
\arguments{
    \item{A_450}{
A R dataframe object that includes the probes, the symbol and synonyms
    }\item{Deconv_methylation}{
A R dataframe with methylation values for each probe and cell type
    }
}

\value{
This function returns a R dataframe object with methylation values for each gene and cell types
}
\author{
Banabithi Bose
}
\note{
This function provides a similar functionality to the ProbeTogene450Fun but with the 450k probes.
}
\examples{
print("Function to convert the probes to genes.")
#load(system.file("extdata", "Methylation.Probe.Annotation450.rda", package = "CTDPathSim"))
#load(system.file("extdata", "Example.Deconv.methylation.450K.rda", package = "CTDPathSim"))
#Deconv_meth_gene<-ProbeTogene450Fun(A_450=A_450,Deconv_methylation=Deconv_methylation)
}

