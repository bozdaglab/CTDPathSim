\name{Annotation450K}
\alias{Annotation450K}
\docType{data}
\title{
R Dataframe or SummarizedExperiament Object with probes and gene names
}
\description{
This is a R dataframe/SummarizedExperiment object that includes methylation probe id's, and the corresponding reference gene name and reference group name. The columns in the dataframe/SummarizedExperiment have to be labeled to follow the convention of the annotation file in "Genome Browser" (https://genome.ucsc.edu/) as "ID, UCSC_REFGENE_NAME, UCSC_REFGENE_GROUP".
}
\usage{data("Annotation450K")}
\format{
  A data frame with 224 observations on the following 3 variables.
  \describe{
    \item{\code{ID}}{a character vector}
    \item{\code{UCSC_REFGENE_NAME}}{a character vector}
    \item{\code{UCSC_REFGENE_GROUP}}{a character vector}
  }
}
\source{
TCGA Database (https://portal.gdc.cancer.gov/)
}
\examples{
data(Annotation450K)
## maybe str(Annotation450K) ; plot(Annotation450K) ...
}
\keyword{datasets}
