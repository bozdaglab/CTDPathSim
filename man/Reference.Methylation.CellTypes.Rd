\name{Reference.Methylation.CellTypes}
\alias{Reference.Methylation.CellTypes}
\docType{data}
\title{
R dataframe/SummarizedExperiment object includes reference methylation profile of reference cell types}
\description{
This R dataframe/SummarizedExperiment object includes reference methylation profile of reference cell types needed for deconvolution. The rows are methylation probes and columns are cell type names.}
\usage{data("Reference.Methylation.CellTypes")}
\format{
  The format is:
 num [1:500, 1:38] 0.872 0.834 0.888 0.86 0.778 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:500] "cg00323915" "cg21830221" "cg11348106" "cg12232731" ...
  ..$ : chr [1:38] "Neu" "NK" "Neu_1" "Bcell" ...
}
\source{
TCGA Database (https://portal.gdc.cancer.gov/)
}
\examples{
data(Reference.Methylation.CellTypes)
## maybe str(Reference.Methylation.CellTypes) ; plot(Reference.Methylation.CellTypes) ...
}
\keyword{datasets}
