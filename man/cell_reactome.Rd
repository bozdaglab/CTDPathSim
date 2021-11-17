\name{cell_reactome}
\alias{cell_reactome}
\docType{data}
\title{
A R dataframe/SummarizedExperiment object of a cell lines's enriched pathways.}
\description{
A R dataframe/SummarizedExperiment object of a cell lines's enriched pathways.}
\usage{data("cell_reactome")}
\format{
  A data frame with 43 observations on the following 3 variables.
  \describe{
    \item{\code{ID}}{a character vector}
    \item{\code{reactome_pathway}}{a character vector}
    \item{\code{p.adjust.cell}}{a character vector}
  }
}
\source{
TCGA Database (https://portal.gdc.cancer.gov/)
}
\examples{
data(cell_reactome)
## maybe str(cell_reactome) ; plot(cell_reactome) ...
}
\keyword{datasets}
