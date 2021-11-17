\name{cell_de_genes}
\alias{cell_de_genes}
\docType{data}
\title{
R dataframe/SummarizedExperiment object of a cell lines's DE genes
}
\description{
A R dataframe/SummarizedExperiment object of a cell lines's DE genes.
}
\usage{data("cell_de_genes")}
\format{
  A data frame with 5080 observations on the following variable.
  \describe{
    \item{\code{gene}}{a character vector}
  }
}
\source{
TCGA Database (https://portal.gdc.cancer.gov/)
}
\examples{
data(cell_de_genes)
## maybe str(cell_de_genes) ; plot(cell_de_genes) ...
}
\keyword{datasets}
