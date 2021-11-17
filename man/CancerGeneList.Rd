\name{CancerGeneList}
\alias{CancerGeneList}
\docType{data}
\title{
Dataframe with cancer-driver genes
}
\description{
This R dataframe object has a single column containing cancer-driver genes.
}
\usage{data("CancerGeneList")}
\format{
  A data frame with 703 observations on the following variable.
  \describe{
    \item{\code{V1}}{a character vector}
  }
}
\source{
TCGA Database (https://portal.gdc.cancer.gov/)
}
\examples{
data(CancerGeneList)
## maybe str(CancerGeneList) ; plot(CancerGeneList) ...
}
\keyword{datasets}
