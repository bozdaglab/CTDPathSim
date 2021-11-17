\name{pat_reactome}
\alias{pat_reactome}
\docType{data}
\title{
A R dataframe/SummarizedExperiment object of a sample's enriched pathways.
}
\description{
A R dataframe/SummarizedExperiment object of a sample's enriched pathways. This file is created in the pipeline}
\usage{data("pat_reactome")}
\format{
  A data frame with 34 observations on the following 3 variables.
\describe{
    \item{\code{ID}}{a numeric vector}
    \item{\code{reactome_pathway}}{a numeric vector}
    \item{\code{p.adjust.pat}}{a numeric vector}
  }
}
\examples{
data(pat_reactome)
## maybe str(pat_reactome) ; plot(pat_reactome) ...
}
\keyword{datasets}
