\name{ProbeTogene27Fun}
\alias{ProbeTogene27Fun}
\title{
Function to convert the probes to genes.
}
\description{
This function converts dataframes with probes to dataframes with genes.
}
\usage{
ProbeTogene27Fun(A_27,Deconv_methylation)
}
\arguments{
    \item{A_27}{
A R dataframe object that includes the probes, the symbol and synonyms
    }\item{Deconv_methylation}{
A R dataframe with methylation values for each probe and cell type
    }
}
\value{
This function returns a R dataframe object with methylation values for each gene and cell type
}
\author{
Banabithi Bose
}
\note{
This function provides a similar functionality to the ProbeTogene450Fun but with the 27k probes.
}
\examples{
print("Function to convert the probes to genes.")
##Deconv_meth_gene<-ProbeTogene27Fun(A_27,Deconv_methylation)
}
