\name{Annotation27K}
\alias{Annotation27K}
\docType{data}
\title{
R Dataframe or SummarizedExperiament Object with probes and gene names
}
\description{
This is a R dataframe/SummarizedExperiment object that includes methylation probe id's, and the corresponding reference gene name and synonyms. The columns are labeled to follow the convention of the annotation file in "Genome Browser" as "ID, Symbol, Synonym".
}
\usage{data("Annotation27K")}
\format{
  A data frame with 27578 observations on the following 3 variables.
  \describe{
    \item{\code{ID}}{a character vector}
    \item{\code{Symbol}}{a character vector}
    \item{\code{Synonym}}{a character vector}
  }
}
\source{
TCGA Database (https://portal.gdc.cancer.gov/)
}
\examples{
data(Annotation27K)
## maybe str(Annotation27K) ; plot(Annotation27K) ...
}
\keyword{datasets}
