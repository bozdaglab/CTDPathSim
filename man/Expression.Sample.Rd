\name{Expression.Sample}
\alias{Expression.Sample}
\docType{data}
\title{
R dataframe with expression data
}
\description{
This R dataframe object has gene expression data for patient samples. The columns are the samples and the genes are in the rows. The user can use RPKM or FPKM normalized data. The tool has not been tested on count data.}
\usage{data("Expression.Sample")}
\format{
  The format is:
 num [1:5, 1:79] 10.08 0.232 40.222 1.314 0.641 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:5] "TSPAN6" "TNMD" "DPM1" "SCYL3" ...
  ..$ : chr [1:79] "TCGA-OR-A5JT-01A" "TCGA-OR-A5K6-01A" "TCGA-OR-A5K3-01A" "TCGA-OR-A5JI-01A" ...
}
\source{
TCGA Database (https://portal.gdc.cancer.gov/)
}
\examples{
data(Expression.Sample)
## maybe str(Expression.Sample) ; plot(Expression.Sample) ...
}
\keyword{datasets}
