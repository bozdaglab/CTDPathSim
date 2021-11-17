\name{ProbeTogeneFun}
\alias{ProbeTogeneFun}
\title{
Function to convert the probes to genes.
}
\description{
This function converts dataframes with probes to dataframes with genes.
}
\arguments{
    \item{A_27}{
A R dataframe object that includes the probes, the symbol and synonyms
    }\item{Deconv_methylation}{
A R dataframe with methylation values for each probe and cell type
    }
}
\value{
This function returns a R dataframe object with methylation 
values for each gene and cell type
}
\author{
Banabithi Bose
}
\note{
This function provides a similar functionality to the 
ProbeTogene450Fun but with the 27k probes.
}
\examples{
load(system.file("CTDData.RData", package="CTDPathSim2"))
Deconv_methylation = stage1_result_ct$methylation
Deconv_meth_gene<-ProbeToGeneFun(Methylation.Probe.Annotation=Annotation450K,Deconv_methylation=Deconv_methylation, Probe.Annotation = "450K")
}
