\name{ccle_expr}
\alias{ccle_expr}
\docType{data}
\title{
A R dataframe/SummarizedExperiment object with a cell lines's expression profile.
}
\description{
A R dataframe/SummarizedExperiment object with a cell lines's expression profile. Rows are the genes and column is the cell line.}
\usage{data("ccle_expr")}
\format{
  A data frame with 54271 observations on the following variable.
  \describe{
    \item{\samp{22RV1_PROSTATE}}{a numeric vector}
  }
}
\source{
TCGA Database (https://portal.gdc.cancer.gov/)
}
\examples{
data(ccle_expr)
## maybe str(ccle_expr) ; plot(ccle_expr) ...
}
\keyword{datasets}
