\name{Reference.CellTypes.Names}
\alias{Reference.CellTypes.Names}
\docType{data}
\title{
 R dataframe/SummarizedExperiment object with the reference cell types used for deconvolution}
\description{
This R dataframe/SummarizedExperiment object has two columns containing the replicates and classes of the reference cell types used for deconvolution.}
\usage{data("Reference.CellTypes.Names")}
\format{
  A data frame with 38 observations on the following 2 variables.
  \describe{
    \item{\code{V1}}{a character vector}
    \item{\code{class}}{a character vector}
  }
}
\source{
TCGA Database (https://portal.gdc.cancer.gov/)
}
\examples{
data(Reference.CellTypes.Names)
## maybe str(Reference.CellTypes.Names) ; plot(Reference.CellTypes.Names) ...
}
\keyword{datasets}
